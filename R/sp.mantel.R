# Copyright 2002-2010 by Roger Bivand 
#

sp.mantel.mc <- function(var, listw, nsim, type="moran", zero.policy=NULL,
	alternative="greater", spChk=NULL, return_boot=FALSE) {
        if (is.null(zero.policy))
            zero.policy <- get("zeroPolicy", envir = .spdepOptions)
        stopifnot(is.logical(zero.policy))
        stopifnot(is.vector(var))
	alternative <- match.arg(alternative, c("greater", "less"))
	if(!inherits(listw, "listw")) stop(paste(deparse(substitute(listw)),
		"is not a listw object"))
	if(!is.numeric(var)) stop(paste(deparse(substitute(var)),
		"is not a numeric vector"))
	if(missing(nsim)) stop("nsim must be given")
	n <- length(listw$neighbours)
	if (n < 1) stop("non-positive n")
        gamres <- suppressWarnings(nsim > gamma(n + 1))
        if (gamres) stop("nsim too large for this number of observations")
	if (nsim < 1) stop("non-positive nsim")

	if (any(is.na(var))) stop("NA in var")
	if (n != length(var)) stop("objects of different length")
	if (is.null(spChk)) spChk <- get.spChkOption()
	if (spChk && !chkIDs(var, listw))
		stop("Check of data and weights ID integrity failed")
	
	listw.U <- listw2U(listw)
	mantel.moran <- function(x, listwU, zero.policy) {
		res <- x * lag.listw(listwU, x, zero.policy=zero.policy)
		res <- sum(res)
		res
	}
	mantel.geary <- function(x, listwU, zero.policy) {
		res <- geary.intern(x, listwU, length(x), 
			zero.policy=zero.policy, type="geary")
		res <- sum(res)
		res
	} 
	mantel.sokal <- function(x, listwU, zero.policy) {
		res <- geary.intern(x, listwU, length(x), 
			zero.policy=zero.policy, type="sokal")
		res <- sum(res)
		res
	}
	if (type == "moran") f <- mantel.moran
	else if (type == "geary") f <- mantel.geary
	else if (type == "sokal") f <- mantel.sokal
	else stop("unknown type")
	xs <- scale(var)
        if (return_boot) {
            mantel_boot <- function(var, i, ...) {
                var <- var[i]
                return(f(x=var, ...))
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
            res <- boot(xs, statistic=mantel_boot, R=nsim,
                sim="permutation", listwU=listw.U, 
                zero.policy=zero.policy, parallel=parallel, ncpus=ncpus, cl=cl)
            return(res)
        }

	res <- numeric(length=nsim+1)
	for (i in 1:nsim) {
		y <- sample(xs)
		res[i] <- f(y, listw.U, zero.policy=zero.policy)
	}
	res[nsim+1] <- f(xs, listw.U, zero.policy=zero.policy)
	rankres <- rank(res)
	xrank <- rankres[length(res)]
	diff <- nsim - xrank
	diff <- ifelse(diff > 0, diff, 0)
        if (alternative == "less") 
        	pval <- punif((diff + 1)/(nsim + 1), lower.tail=FALSE)
    	else if (alternative == "greater") 
        	pval <- punif((diff + 1)/(nsim + 1))
	if (!is.finite(pval) || pval < 0 || pval > 1) 
		warning("Out-of-range p-value: reconsider test arguments")
	statistic <- res[nsim+1]
	names(statistic) <- "statistic"
	parameter <- xrank
	names(parameter) <- "observed rank"
	method <- paste("Mantel permutation test for", type, "measure")
	data.name <- paste(deparse(substitute(var)), "\nweights:",
	    deparse(substitute(listw)), "\nnumber of simulations + 1:",
	    nsim+1, "\n")
	est <- c(mean(res[1:nsim]), sd(res[1:nsim]))
	names(est) <- c("mean of permutations", "sd of permutations")
	lres <- list(statistic=statistic, parameter=parameter,
	    alternative=alternative, method=method, data.name=data.name, 
	    p.value=pval, res=res, estimate=est)
	class(lres) <- c("htest", "mc.sim")
	lres
}

plot.mc.sim <- function(x, xlim, xlab, main, sub, ..., ptype="density") {
	res <- x$res
	if (missing(xlim)) xlim <- range(res)*1.1
	n <- length(res)
	obs <- res[n]
	res <- res[-n]
        if (missing(xlab)) xlab <- strsplit(x$data.name, "\n")[[1]][1]
        if (missing(sub)) sub <- x$method
	if(ptype == "density") {
            if (missing(main)) main <- "Density plot of permutation outcomes"
            plot(density(res), xlim=xlim, xlab=xlab, main=main, sub=sub, ...)
        } else {
            if (missing(main)) main <- "Histogram of permutation outcomes"
            hist(res, xlim=xlim, xlab=xlab, main=main, sub=sub, ...)
        }
	abline(v=obs, col=1, lwd=2)
}
