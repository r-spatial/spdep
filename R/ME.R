# Copyright 2005-2008 by Roger Bivand and Pedro Peres-Neto (from Matlab)
#

ME <- function(formula, data, family = gaussian, weights, offset, listw, 
	alpha=0.05, nsim=99, verbose=NULL, stdev=FALSE) {
	MoraneI.boot <- function(var, i, ...) {
		var <- var[i]
		I <- (n/S0)*(crossprod(sW %*% var, var))/cpvar
		return(c(as(I, "matrix")))
	}

	MIR_a <- function(resids, sW, n, cpvar, S0, nsim, stdev=TRUE,
            par_boot_args=list()) {
		boot1 <- boot(resids, statistic=MoraneI.boot, R=nsim, 
			sim="permutation", sW=sW, n=n, S0=S0, cpvar=cpvar,
                        parallel=par_boot_args$parallel,
                        ncpus=par_boot_args$ncpus, cl=par_boot_args$cl)
		mi <- boot1$t0
		if (stdev) {
			zi <- (boot1$t0 - mean(boot1$t))/sqrt(var(boot1$t))
			pri <- pnorm(abs(zi), lower.tail=FALSE)
		} else {
			zi <- NA
			pri <- (sum(boot1$t >= mi)+1)/(nsim+1)
		}
		res <- list(estimate=mi, statistic=zi, p.value=pri)
		res
	}
        if (is.null(verbose)) verbose <- get("verbose", envir = .spdepOptions)
        stopifnot(is.logical(verbose))

	listw <- listw2U(listw) # make weights symmetric
	sW <- as(listw, "CsparseMatrix")
	
	Wmat <- listw2mat(listw) # convert to full matrix form
	n <- ncol(Wmat)
	S0 <- Szero(listw)

# argument handling copied from stats:::glm
#	call <- match.call()
    	if (is.character(family)) 
            family <- get(family, mode = "function", envir = parent.frame())
    	if (is.function(family)) family <- family()
    	if (is.null(family$family)) {
            print(family)
            stop("'family' not recognized")
    	}
     	if (missing(data)) data <- environment(formula)
    	mf <- match.call(expand.dots = FALSE)
    	m <- match(c("formula", "data", "weights", "offset"), names(mf), 0)
    	mf <- mf[c(1, m)]
    	mf$drop.unused.levels <- TRUE
    	mf[[1]] <- as.name("model.frame")
    	mf <- eval(mf, parent.frame())

	mt <- attr(mf, "terms")
	Y <- model.extract(mf, "response") # extract X and Y
	X <- model.matrix(mt, mf)

	weights <- model.weights(mf)
	if (!is.null(weights) && any(weights < 0)) 
            stop("negative weights not allowed")

	offset <- model.offset(mf)
	if (!is.null(offset) && length(offset) != NROW(Y))
	    stop("number of offsets should equal number of observations")

	glm_fit <- glm.fit(x=X, y=Y, weights=weights, offset=offset, 
	    family=family)
	glm_res <- glm_fit$y - glm_fit$fitted.values
	cpvar <- crossprod(glm_res)
        cores <- get.coresOption()
        if (is.null(cores)) {
        parallel <- "no"
        } else {
            parallel <- ifelse (get.mcOption(), "multicore", "snow")
        }
        ncpus <- ifelse(is.null(cores), 1L, cores)
        cl <- NULL
        if (parallel == "snow") {
            cl <- get.ClusterOption()
            if (is.null(cl)) {
                parallel <- "no"
                warning("no cluster in ClusterOption, parallel set to no")
            }
        }
        par_boot_args <- list(parallel=parallel, ncpus=ncpus, cl=cl)
	mRES <- MIR_a(glm_res, sW=sW, n=n, cpvar=cpvar, S0=S0, nsim=nsim,
		stdev=stdev, par_boot_args=par_boot_args)
	pIZ <- mRES$p.value
	tres <- c(NA, mRES$statistic, pIZ)
	if (pIZ > alpha) stop("base correlation larger than alpha")

	Cent <- diag(n) - (matrix(1/n, n, n))
        CWC <- Cent %*% Wmat %*% Cent
	rm(Cent, Wmat)
        CWC2 <- 0.5*(CWC+t(CWC))
        rm(CWC)
	eV <- eigen(CWC2)$vectors
        rm(CWC2)
	iZ <- numeric(n)
	for (i in 1:n) {
		iX <- cbind(X, eV[,i])
		i_glm <- glm.fit(x=iX, y=Y, weights=weights, offset=offset, 
	    		family=family)
		glm_res <- i_glm$y - i_glm$fitted.values
		cpvar <- crossprod(glm_res)
#		iZ[i] <- MIR_a(glm_res, sW=sW, n=n, cpvar=cpvar, S0=S0, 
#			nsim=nsim)$statistic
		iZ[i] <- c(as((n/S0)*(crossprod(sW %*% glm_res, glm_res)) /
			cpvar, "matrix"))
	}
	min_iZ <- which.min(abs(iZ))
	X <- cbind(X, eV[, min_iZ])
	glm_fit <- glm.fit(x=X, y=Y, weights=weights, offset=offset, 
	    family=family)
	glm_res <- glm_fit$y - glm_fit$fitted.values
	cpvar <- crossprod(glm_res)
	mRES <- MIR_a(glm_res, sW=sW, n=n, cpvar=cpvar, S0=S0, nsim=nsim,
		stdev=stdev, par_boot_args=par_boot_args)
	pIZ <- mRES$p.value
	used <- rep(FALSE, n)
	used[min_iZ] <- TRUE
	min_v <- min_iZ
	if (verbose) cat("eV[,", min_iZ, "], I: ", mRES$estimate, " ZI: ", 
		mRES$statistic, ", pr(ZI): ", pIZ, "\n", sep="")
	tres <- rbind(tres, c(min_iZ, mRES$statistic, pIZ))
	while (pIZ <= alpha) {
		for (i in 1:n) {
		    if (!used[i]) {
			iX <- cbind(X, eV[,i])
			i_glm <- glm.fit(x=iX, y=Y, weights=weights, 
				offset=offset, family=family)
			glm_res <- i_glm$y - i_glm$fitted.values
			cpvar <- crossprod(glm_res)
#			iZ[i] <- MIR_a(glm_res, sW=sW, n=n, cpvar=cpvar, S0=S0,
#				nsim=nsim)$statistic
			iZ[i] <- c(as((n/S0)*(crossprod(sW %*% glm_res, 
				glm_res))/cpvar, "matrix"))
		    } else iZ[i] <- NA
		}
		min_iZ <- which.min(abs(iZ))
		X <- cbind(X, eV[, min_iZ])
		glm_fit <- glm.fit(x=X, y=Y, weights=weights, offset=offset, 
	    		family=family)
		glm_res <- glm_fit$y - glm_fit$fitted.values
		cpvar <- crossprod(glm_res)
		mRES <- MIR_a(glm_res, sW=sW, n=n, cpvar=cpvar, S0=S0, 
			nsim=nsim, stdev=stdev, par_boot_args=par_boot_args)
		pIZ <- mRES$p.value
		used[min_iZ] <- TRUE
		min_v <- c(min_v, min_iZ)
		if (verbose) cat("eV[,", min_iZ, "], I: ", mRES$estimate, 
			" ZI: ", mRES$statistic, ", pr(ZI): ", pIZ, 
			"\n", sep="")
		tres <- rbind(tres, c(min_iZ, mRES$statistic, pIZ))
	}
	sv <- eV[,min_v, drop=FALSE]
	colnames(sv) <- paste("vec", min_v, sep="")
	colnames(tres) <- c("Eigenvector", "ZI", "pr(ZI)")
	rownames(tres) <- 0:(nrow(tres)-1)
	res <- list(selection=tres, vectors=sv)
	class(res) <- "ME_res"
	res
}

print.ME_res <- function(x, ...) {
	print(x$selection)
}

fitted.ME_res <- function(object, ...) {
	object$vectors
}


