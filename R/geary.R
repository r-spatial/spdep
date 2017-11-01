# Copyright 2001-9 by Roger Bivand 
#


geary <- function(x, listw, n, n1, S0, zero.policy=NULL) {
        if (is.null(zero.policy))
            zero.policy <- get("zeroPolicy", envir = .spdepOptions)
        stopifnot(is.logical(zero.policy))
        stopifnot(is.vector(x))
	z <- scale(x, scale=FALSE)
	zz <- sum(z^2)
	K <- (n*sum(z^4))/(zz^2)
	res <- geary.intern(x, listw, n, zero.policy, type="geary")
	C <- (n1 / (2*S0)) * (sum(res) / zz)
	res <- list(C=C, K=K)
	res
}

geary.intern <- function(x, listw, n, zero.policy=NULL, type="geary") {
        if (is.null(zero.policy))
            zero.policy <- get("zeroPolicy", envir = .spdepOptions)
        stopifnot(is.logical(zero.policy))
	cardnb <- card(listw$neighbours)
	if (type == "geary") ft <- TRUE
	else if (type == "sokal") ft <- FALSE
	else stop("type unknown")
	res <- .Call("gearyw", listw$neighbours, listw$weights,
		as.numeric(x), as.integer(cardnb),
		as.logical(zero.policy), as.logical(ft), PACKAGE="spdep")
	if (any(is.na(res))) warning("NAs in lagged values")
	res
}

geary.test <- function(x, listw, randomisation=TRUE, zero.policy=NULL,
    alternative="greater", spChk=NULL, adjust.n=TRUE) {
        if (is.null(zero.policy))
            zero.policy <- get("zeroPolicy", envir = .spdepOptions)
        stopifnot(is.logical(zero.policy))
	alternative <- match.arg(alternative, c("less", "greater", "two.sided"))
	if(!inherits(listw, "listw")) stop(paste(deparse(substitute(listw)),
		"is not a listw object"))
	if(!is.numeric(x)) stop(paste(deparse(substitute(x)),
		"is not a numeric vector"))
	if (any(is.na(x))) stop("NA in X")
	n <- length(listw$neighbours)
	if (n != length(x)) stop("objects of different length")
	if (is.null(spChk)) spChk <- get.spChkOption()
	if (spChk && !chkIDs(x, listw))
		stop("Check of data and weights ID integrity failed")
	
	wc <- spweights.constants(listw, zero.policy, adjust.n=adjust.n)
	S02 <- wc$S0*wc$S0
	res <- geary(x, listw, wc$n, wc$n1, wc$S0, zero.policy)
	C <- res$C
	if (is.na(C)) stop("NAs generated in geary - check zero.policy")
	K <- res$K
	EC <- 1
	if(randomisation) {
		VC <- (wc$n1*wc$S1*(wc$nn - 3*n + 3 - K*wc$n1))
		VC <- VC - ((1/4) * (wc$n1*wc$S2*(wc$nn + 3*n - 6 - 
			K*(wc$nn - n + 2))))
		VC <- VC + (S02*(wc$nn - 3 - K*(wc$n1^2)))
		VC <- VC / (n*wc$n2*wc$n3*S02)
	} else {
		VC <- ((2*wc$S1 + wc$S2)*wc$n1 - 4*S02) / (2*(n + 1)*S02)
	}
#	ZC <- (C - EC) / sqrt(VC)
# order changed 090609 RSB (C&O 1973, p. 21)
	ZC <- (EC - C) / sqrt(VC)
	statistic <- ZC
	names(statistic) <- "Geary C statistic standard deviate"
	PrC <- NA
	if (is.finite(ZC)) {
        	if (alternative == "two.sided") PrC <- 2 * pnorm(abs(ZC), 
			lower.tail=FALSE)
        	else if (alternative == "greater")
            	PrC <- pnorm(ZC, lower.tail=FALSE)
        	else PrC <- pnorm(ZC)
		if (!is.finite(PrC) || PrC < 0 || PrC > 1) 
		    warning("Out-of-range p-value: reconsider test arguments")
	}
	vec <- c(C, EC, VC)
	names(vec) <- c("Geary C statistic", "Expectation", "Variance")
	method <- paste("Geary C test under", ifelse(randomisation,
	    "randomisation", "normality"))
	data.name <- paste(deparse(substitute(x)), "\nweights:",
	    deparse(substitute(listw)), "\n")
	res <- list(statistic=statistic, p.value=PrC, estimate=vec, 
	    alternative=ifelse(alternative == "two.sided", alternative, 
	    paste("Expectation", alternative, "than statistic")), 
	    method=method, data.name=data.name)
	class(res) <- "htest"
	res
}

geary.mc <- function(x, listw, nsim, zero.policy=NULL,
	alternative="greater", spChk=NULL, adjust.n=TRUE, return_boot=FALSE) {
        if (is.null(zero.policy))
            zero.policy <- get("zeroPolicy", envir = .spdepOptions)
        stopifnot(is.logical(zero.policy))
        stopifnot(is.vector(x))
	alternative <- match.arg(alternative, c("less", "greater"))
	if(!inherits(listw, "listw")) stop(paste(deparse(substitute(listw)),
		"is not a listw object"))
	if(!is.numeric(x)) stop(paste(deparse(substitute(x)),
		"is not a numeric vector"))
	if(missing(nsim)) stop("nsim must be given")
	if (any(is.na(x))) stop("NA in X")
	n <- length(listw$neighbours)
	if (n != length(x)) stop("objects of different length")
	if (is.null(spChk)) spChk <- get.spChkOption()
	if (spChk && !chkIDs(x, listw))
		stop("Check of data and weights ID integrity failed")
        gamres <- suppressWarnings(nsim > gamma(n + 1))
        if (gamres) stop("nsim too large for this number of observations")
	if (nsim < 1) stop("non-positive nsim")
	wc <- spweights.constants(listw, zero.policy, adjust.n=adjust.n)
        if (return_boot) {
            geary_boot <- function(var, i, ...) {
                var <- var[i]
                return(geary(x=var, ...)$C)
            }
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
            res <- boot(x, statistic=geary_boot, R=nsim,
                sim="permutation", listw=listw, n=n, n1=wc$n1, S0=wc$S0, 
                zero.policy=zero.policy, parallel=parallel, ncpus=ncpus, cl=cl)
            return(res)
        }
	res <- numeric(length=nsim+1)
	for (i in 1:nsim) res[i] <- geary(sample(x), listw, n, wc$n1, wc$S0,
	    zero.policy)$C
	res[nsim+1] <- geary(x, listw, n, wc$n1, wc$S0, zero.policy)$C
	rankres <- rank(res)
	xrank <- rankres[length(res)]
	diff <- nsim - xrank
	diff <- ifelse(diff > 0, diff, 0)
# order changed 110411 RSB (C&O 1973, p. 21) Thanks to Daniel Garavito
        if (alternative == "greater") 
        	pval <- punif((diff + 1)/(nsim + 1), lower.tail=FALSE)
    	else if (alternative == "less") 
        	pval <- punif((diff + 1)/(nsim + 1))
	if (!is.finite(pval) || pval < 0 || pval > 1) 
		warning("Out-of-range p-value: reconsider test arguments")
	statistic <- res[nsim+1]
	names(statistic) <- "statistic"
	parameter <- xrank
	names(parameter) <- "observed rank"
	method <- "Monte-Carlo simulation of Geary C"
	data.name <- paste(deparse(substitute(x)), "\nweights:",
	    deparse(substitute(listw)), "\nnumber of simulations + 1:",
	    nsim+1, "\n")
	lres <- list(statistic=statistic, parameter=parameter,
	    p.value=pval, alternative=alternative, method=method, 
	    data.name=data.name, res=res)
	class(lres) <- c("htest", "mc.sim")
	lres
}

