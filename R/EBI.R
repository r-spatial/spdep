# Copyright 2002-2003 by Roger Bivand and Marilia Carvalho
# Addition of Martuzzi and Elliott Copyright 2006 Olaf Berke and Roger Bivand
#

EBImoran <- function (z, listw, nn, S0, zero.policy = NULL, subtract_mean_in_numerator=TRUE) 
{
#default subtract_mean_in_numerator=TRUE 160219 RA
        if (is.null(zero.policy))
            zero.policy <- get("zeroPolicy", envir = .spdepOptions)
        stopifnot(is.logical(zero.policy))
    stopifnot(is.vector(z))
    zm <- mean(z)
    zz <- sum((z - zm)^2)
    if (subtract_mean_in_numerator) z <- z - zm
    lz <- lag.listw(listw, z, zero.policy = zero.policy)
    EBI <- (nn/S0) * (sum(z * lz)/zz)
    res <- EBI
    res
}

EBImoran.mc <- function (n, x, listw, nsim, zero.policy = NULL,
 alternative = "greater", spChk = NULL, return_boot=FALSE, subtract_mean_in_numerator=TRUE) 
{
#default subtract_mean_in_numerator=TRUE 160219 RA
    message("The default for subtract_mean_in_numerator set TRUE from February 2016")
        if (is.null(zero.policy))
            zero.policy <- get("zeroPolicy", envir = .spdepOptions)
        stopifnot(is.logical(zero.policy))
    alternative <- match.arg(alternative, c("greater", "less"))
    if (!inherits(listw, "listw")) 
        stop(paste(deparse(substitute(listw)), "is not a listw object"))
    if (missing(nsim)) 
        stop("nsim must be given")
    m <- length(listw$neighbours)
    if (m != length(x)) 
        stop("objects of different length")
    if (m != length(n)) 
        stop("objects of different length")
    if (is.null(spChk)) 
        spChk <- get.spChkOption()
    if (spChk && !chkIDs(x, listw)) 
        stop("Check of data and weights ID integrity failed")
    if (spChk && !chkIDs(n, listw)) 
        stop("Check of data and weights ID integrity failed")
    gamres <- suppressWarnings(nsim > gamma(m + 1))
    if (gamres) stop("nsim too large for this number of observations")
    if (nsim < 1) stop ("nsim too small")
    S0 <- Szero(listw)
    EB <- EBest(n, x)
    p <- EB$raw
    b <- attr(EB, "parameters")$b
    a <- attr(EB, "parameters")$a
    v <- a + (b/x)
    v[v < 0] <- b/x
    z <- (p - b)/sqrt(v)
    if (return_boot) {
        EBI_boot <- function(var, i, ...) {
            var <- var[i]
            return(EBImoran(z=var, ...))
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
        res <- boot(z, statistic=EBI_boot, R=nsim,
            sim="permutation", listw=listw, nn=m, S0=S0,
            zero.policy=zero.policy,
            subtract_mean_in_numerator=subtract_mean_in_numerator,
            parallel=parallel, ncpus=ncpus, cl=cl)
        return(res)
    }
    res <- numeric(length = nsim + 1)
    for (i in 1:nsim) res[i] <- EBImoran(sample(z), listw, m, 
        S0, zero.policy, subtract_mean_in_numerator=subtract_mean_in_numerator)
    res[nsim + 1] <- EBImoran(z, listw, m, S0, zero.policy,
            subtract_mean_in_numerator=subtract_mean_in_numerator)
    rankres <- rank(res)
    zrank <- rankres[length(res)]
    diff <- nsim - zrank
    diff <- ifelse(diff > 0, diff, 0)
    if (alternative == "less") 
        pval <- punif((diff + 1)/(nsim + 1), lower.tail=FALSE)
    else if (alternative == "greater") 
        pval <- punif((diff + 1)/(nsim + 1))
    if (!is.finite(pval) || pval < 0 || pval > 1) 
	warning("Out-of-range p-value: reconsider test arguments")
    statistic <- res[nsim + 1]
    names(statistic) <- "statistic"
    parameter <- zrank
    names(parameter) <- "observed rank"
    method <- "Monte-Carlo simulation of Empirical Bayes Index"
    if (subtract_mean_in_numerator) method <- paste(method, "(mean subtracted)")
    data.name <- paste("cases: ", deparse(substitute(n)),
        ", risk population: ", deparse(substitute(x)), "\nweights: ",
        deparse(substitute(listw)), "\nnumber of simulations + 1: ",
        nsim + 1, "\n", sep = "")
    lres <- list(statistic = statistic, parameter = parameter, 
        p.value = pval, alternative = alternative, method = method, 
        data.name = data.name, res = res, z = z)
    class(lres) <- c("htest", "mc.sim")
    lres
}

probmap <- function(n, x, row.names=NULL, alternative="less") {
    alternative <- match.arg(alternative, c("greater", "less"))
    if (!is.numeric(x)) 
        stop(paste(deparse(substitute(x)), "is not a numeric vector"))
    if (!is.numeric(n)) 
        stop(paste(deparse(substitute(n)), "is not a numeric vector"))
    if (any(is.na(x))) 
        stop("NA in at risk population")
    if (any(is.na(n))) 
        stop("NA in cases")
    if (any(x < 0)) 
        stop("negative risk population")
    if (any(n < 0)) 
        stop("negative number of cases")
    p <- n/x
    nsum <- sum(n)
    xsum <- sum(x)
    b <- nsum/xsum
    expCount <- x*b
    relRisk <- 100*(n/expCount)
#    pmap <- ppois(n, expCount, lower.tail=(alternative=="less"))
    if (alternative == "less") {
        pmap <- ppois(n, expCount)
    } else {
        pmap <- 1 - ppois(n-1, expCount)
    }

    if (is.null(row.names)) 
	res <- data.frame(raw=p, expCount=expCount, relRisk=relRisk, 
	pmap=pmap)
    else
    	res <- data.frame(raw=p, expCount=expCount, relRisk=relRisk, 
	pmap=pmap, row.names=row.names)
    res
}


EBest <- function(n, x, family="poisson") {
    if (!is.numeric(x)) 
        stop(paste(deparse(substitute(x)), "is not a numeric vector"))
    if (!is.numeric(n)) 
        stop(paste(deparse(substitute(n)), "is not a numeric vector"))
    if (any(is.na(x))) 
        stop("NA in at risk population")
    if (any(is.na(n))) 
        stop("NA in cases")
    if (any(x < .Machine$double.eps)) 
        stop("non-positive risk population")
    if (any(n < 0)) 
        stop("negative number of cases")
    if (length(x) != length(n)) stop("vectors of different length")
    m <- length(n)
    p <- n/x
    nsum <- sum(n)
    xsum <- sum(x)
    b <- nsum/xsum
    s2 <- sum(x * (((p - b)^2)/xsum))
    if (family == "poisson") {
        a <- s2 - (b/(xsum/m))
        if (a < 0) a <- 0
        est <- b + (a*(p - b)) / (a + (b/x))
        res <- data.frame(raw=p, estmm=est)
        attr(res, "family") <- family
        attr(res, "parameters") <- list(a=a, b=b)
    } else if (family == "binomial") {
# contributed by Olaf Berke
	rho <- (x*s2 - (x/mean(x))*(b*(1-b))) / 
		((x-1)*s2 + ((mean(x)-x) /mean(x))*(b*(1-b)))
	est <- rho*p + (1-rho)*b
	res <- data.frame(raw = p, estmm = est)
        attr(res, "family") <- family
	attr(res, "parameters") <- list(a=s2, b=b)
    } else stop("family unknown")
    res
}

#EBestB <- function (n,x)
#{
#p <- n/x
#b <- sum(n)/sum(x)
#    nsum <- sum(n)
#    xsum <- sum(x)
#    b <- nsum/xsum
#s2 <- sum(x*(p-b)^2) / xsum
#roh <- (x*s2 - (x/mean(x))*(b*(1-b))) / 
#  ((x-1)*s2 + ((mean(x)-x) /mean(x))*(b*(1-b)))
#est <- roh*p + (1-roh)*b
#result <- data.frame(raw = p, ebest = est)
#attr(result, "parameters") <- list(a = b, b = s2)
#result
#}

EBlocal <- function(ri, ni, nb, zero.policy = NULL,
    spChk = NULL, geoda = FALSE) {
        if (is.null(zero.policy))
            zero.policy <- get("zeroPolicy", envir = .spdepOptions)
        stopifnot(is.logical(zero.policy))
# class to inherits Jari Oksanen 080603
    if (!inherits(nb, "nb")) 
        stop(paste(deparse(substitute(nb)), "is not an nb object"))
    lnb <- length(nb)
    if (lnb < 1) stop("zero length neighbour list")
    if (lnb != length(ri)) 
        stop("objects of different length")
    if (lnb != length(ni)) 
        stop("objects of different length")
    if (is.null(spChk)) 
        spChk <- get.spChkOption()
    if (spChk && !chkIDs(ni, nb)) 
        stop("Check of data and neighbour ID integrity failed")
    if (spChk && !chkIDs(ri, nb)) 
        stop("Check of data and neighbour ID integrity failed")
    if (!is.numeric(ni)) 
        stop(paste(deparse(substitute(ni)), "is not a numeric vector"))
    if (!is.numeric(ri)) 
        stop(paste(deparse(substitute(ri)), "is not a numeric vector"))
    if (any(is.na(ni))) 
        stop("NA in at risk population")
    if (any(ni == 0)) 
        stop("zero in at risk population")
    if (any(is.na(ri))) 
        stop("NA in cases")
    if (any(ni < 0)) 
        stop("negative risk population")
    if (any(ri < 0)) 
        stop("negative number of cases")
    lw <- nb2listw(include.self(nb), style="B", zero.policy=zero.policy)
    xi <- ri/ni
    r.i <- lag.listw(lw, ri, zero.policy = zero.policy)
    n.i <- lag.listw(lw, ni, zero.policy = zero.policy)
    nbar.i <- lag.listw(nb2listw(include.self(nb), style="W",
        zero.policy=zero.policy), ni, zero.policy = zero.policy)

    m.i <- r.i/n.i
    if (geoda) {
	work <- vector(mode="list", length=lnb)
	cnb <- card(lw$neighbours)
	for (i in 1:lnb) {
		work[[i]] <- rep(m.i[i], cnb[i])
		work[[i]] <- xi[lw$neighbours[[i]]] - work[[i]]
		work[[i]] <- work[[i]]^2
		work[[i]] <- ni[lw$neighbours[[i]]] * work[[i]]
	}
	C.i <- sapply(work, sum)
    } else {
        C.i <- lag.listw(lw, (ni * (xi - m.i)^2), zero.policy = zero.policy)
    }
    a.i <- (C.i/n.i) - (m.i/nbar.i)
    a.i[a.i < 0] <- 0
    est <- m.i + (xi - m.i) * (a.i / (a.i + (m.i/ni)))

    res <- data.frame(raw=xi, est=est)
    attr(res, "parameters") <- list(a=a.i, m=m.i)
    res
}
