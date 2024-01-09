# Copyright 2001-24 by Roger Bivand 
#


geary <- function(x, listw, n, n1, S0, zero.policy=attr(listw, "zero.policy")) {
        if (is.null(zero.policy))
            zero.policy <- get("zeroPolicy", envir = .spdepOptions)
        stopifnot(is.logical(zero.policy))
        stopifnot(is.vector(x))
        stopifnot(all(is.finite(x)))
#	z <- scale(x, scale=FALSE)
        z <- scale(x)
	zz <- sum(z^2)
	K <- (n*sum(z^4))/(zz^2)
#	res <- geary.intern(x, listw, n, zero.policy, type="geary")
	res <- geary.intern(z, listw, n, zero.policy, type="geary")
	C <- (n1 / (2*S0)) * (sum(res) / zz)
	res <- list(C=C, K=K)
	res
}

geary.intern <- function(x, listw, n, zero.policy=attr(listw, "zero.policy"), type="geary") {
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
#	if (any(is.na(res))) warning("NAs in lagged values")
	res
}

geary.test <- function(x, listw, randomisation=TRUE, zero.policy=attr(listw, "zero.policy"),
    alternative="greater", spChk=NULL, adjust.n=TRUE, na.action=na.fail) {
        if (is.null(zero.policy))
            zero.policy <- get("zeroPolicy", envir = .spdepOptions)
        stopifnot(is.logical(zero.policy))
	alternative <- match.arg(alternative, c("less", "greater", "two.sided"))
	wname <- deparse(substitute(listw))
	if (!inherits(listw, "listw")) stop(wname, "is not a listw object")
        xname <- deparse(substitute(x))
	if(!is.numeric(x)) stop(xname,
		" is not a numeric vector")
	if (deparse(substitute(na.action)) == "na.pass")
	    stop("na.pass not permitted")
	x <- na.action(x)
	na.act <- attr(x, "na.action")
	if (!is.null(na.act)) {
	    subset <- !(1:length(listw$neighbours) %in% na.act)
	    listw <- subset(listw, subset, zero.policy=zero.policy)
	}
	n <- length(listw$neighbours)
	if (n != length(x)) stop("objects of different length")
	if (is.null(spChk)) spChk <- get.spChkOption()
	if (spChk && !chkIDs(x, listw))
		stop("Check of data and weights ID integrity failed")
	
	wc <- spweights.constants(listw, zero.policy, adjust.n=adjust.n)
	S02 <- wc$S0*wc$S0
	res <- geary(as.vector(x), listw, wc$n, wc$n1, wc$S0, zero.policy)
	C <- res$C
	if (is.na(C)) stop("NAs generated in geary - check zero.policy and na.action")
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
	data.name <- paste(xname, "\nweights:", wname,
            ifelse(is.null(na.act), "", paste("\nomitted:", 
	    paste(na.act, collapse=", "))),
            ifelse(adjust.n && isTRUE(any(sum(card(listw$neighbours) == 0L))),
            "\nn reduced by no-neighbour observations", ""), "\n")
	res <- list(statistic=statistic, p.value=PrC, estimate=vec, 
	    alternative=ifelse(alternative == "two.sided", alternative, 
	    paste("Expectation", alternative, "than statistic")), 
	    method=method, data.name=data.name)
	if (!is.null(na.act)) attr(res, "na.action") <- na.act
	class(res) <- "htest"
	res
}

geary.mc <- function(x, listw, nsim, zero.policy=attr(listw, "zero.policy"),
	alternative="greater", spChk=NULL, adjust.n=TRUE, return_boot=FALSE,
        na.action=na.fail) {
        if (is.null(zero.policy))
            zero.policy <- get("zeroPolicy", envir = .spdepOptions)
        stopifnot(is.logical(zero.policy))
        stopifnot(is.vector(x))
	alternative <- match.arg(alternative, c("less", "greater", "two.sided"))
	wname <- deparse(substitute(listw))
	if (!inherits(listw, "listw")) stop(wname, "is not a listw object")
        xname <- deparse(substitute(x))
	if(!is.numeric(x)) stop(xname, " is not a numeric vector")
	if(missing(nsim)) stop("nsim must be given")
	if (deparse(substitute(na.action)) == "na.pass")
	    stop("na.pass not permitted")
	x <- na.action(x)
	na.act <- attr(x, "na.action")
	if (!is.null(na.act)) {
	    subset <- !(1:length(listw$neighbours) %in% na.act)
	    listw <- subset(listw, subset, zero.policy=zero.policy)
            if (return_boot) 
              message("NA observations omitted: ", paste(na.act, collapse=", "))
	}
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
            p_setup <- parallel_setup(NULL)
            parallel <- p_setup$parallel
            ncpus <- p_setup$ncpus
            cl <- p_setup$cl
            res <- boot(x, statistic=geary_boot, R=nsim,
                sim="permutation", listw=listw, n=n, n1=wc$n1, S0=wc$S0, 
                zero.policy=zero.policy, parallel=parallel, ncpus=ncpus, cl=cl)
            return(res)
        }
	res <- numeric(length=nsim+1)
	for (i in 1:nsim) res[i] <- geary(sample(x), listw, n, wc$n1, wc$S0,
	    zero.policy)$C
	res[nsim+1] <- geary(as.vector(x), listw, n, wc$n1, wc$S0, zero.policy)$C
	rankres <- rank(res)
	xrank <- rankres[length(res)]
	diff <- nsim - xrank
	diff <- ifelse(diff > 0, diff, 0)
# order changed 110411 RSB (C&O 1973, p. 21) Thanks to Daniel Garavito
        if (alternative == "greater") 
        	pval <- punif((diff + 1)/(nsim + 1), lower.tail=FALSE)
    	else if (alternative == "less") 
        	pval <- punif((diff + 1)/(nsim + 1))
        else pval <- punif(abs(xrank - (nsim+1)/2)/(nsim + 1), 0, 0.5,
                lower.tail=FALSE)
	if (!is.finite(pval) || pval < 0 || pval > 1) 
		warning("Out-of-range p-value: reconsider test arguments")
	statistic <- res[nsim+1]
	names(statistic) <- "statistic"
	parameter <- xrank
	names(parameter) <- "observed rank"
	method <- "Monte-Carlo simulation of Geary C"
	data.name <- paste(xname, "\nweights:", wname,
            ifelse(is.null(na.act), "", paste("\nomitted:", paste(na.act, 
                collapse=", "))), "\nnumber of simulations + 1:", nsim+1, "\n")
	lres <- list(statistic=statistic, parameter=parameter,
	    p.value=pval, alternative=alternative, method=method, 
	    data.name=data.name, res=res)
	if (!is.null(na.act)) attr(res, "na.action") <- na.act
	class(lres) <- c("htest", "mc.sim")
	lres
}

