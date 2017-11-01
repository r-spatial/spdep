# Copyright (c) 2007-2008 Markus Reder and Roger Bivand

localmoran.exact <- function(model, select, nb, glist = NULL, style = "W",
    zero.policy = NULL, alternative = "greater", spChk=NULL, 
    resfun=weighted.residuals, save.Vi = FALSE, useTP=FALSE, truncErr=1e-6,
    zeroTreat=0.1) {
# need to impose check on weights TODO!!
# class to inherits Jari Oksanen 080603
  if (!inherits(nb, "nb"))
        stop(paste(deparse(substitute(nb)), "not an nb object"))
        if (is.null(zero.policy))
            zero.policy <- get("zeroPolicy", envir = .spdepOptions)
        stopifnot(is.logical(zero.policy))
    n <- length(nb)
    u <- resfun(model)
    if (n != length(u)) 
        stop("objects of different length")
    if (is.null(spChk)) spChk <- get.spChkOption()
    if (spChk && !chkIDs(u, nb2listw(nb, zero.policy=zero.policy)))
	stop("Check of data and weights ID integrity failed")
    if (!(alternative %in% c("greater", "less", "two.sided")))
	stop("alternative must be one of: \"greater\", \"less\", or \"two.sided\"")
    if (missing(select)) select <- 1:n
    u <- as.vector(u)
    select <- unique(as.integer(select))
    if (length(select) < 1L) stop("select too short")
    if (any(select < 1 || select > n))
        stop("select out of range")
    utu <- c(crossprod(u))
    p <- model$rank
    p1 <- 1:p
    nacoefs <- which(is.na(coefficients(model)))
    m <- n - p - 2
    XtXinv <- chol2inv(model$qr$qr[p1, p1, drop = FALSE])
    X <- model.matrix(terms(model), model.frame(model))
# fixed after looking at TOWN dummy in Boston data
    if (length(nacoefs) > 0L) X <- X[,-nacoefs]
    if (!is.null(wts <- weights(model))) {
	X <- sqrt(diag(wts)) %*% X
    }
    B <- listw2U(nb2listw(nb, glist=glist, style="B",
	zero.policy=zero.policy))
    D <- NULL
    a <- NULL
    if (style == "W") {
        D <- 1/sapply(B$weights, sum)
    } else if (style == "S") {
        D <- 1 / sqrt(sapply(B$weights, function(x) sum(x^2)))
# correction by Danlin Yu, 25 March 2004
	a <- sum(sapply(B$weights, function(x) sqrt(sum(x^2))))
    } else if (style == "C") a <- sum(unlist(B$weights))
    res <- vector(mode="list", length=length(select))
    for (i in 1:length(select)) {
        Vi <- listw2star(B, select[i], style=style, n, D, a,
	    zero.policy=zero.policy)
        Viu <- lag.listw(Vi, u, zero.policy=TRUE)
	Ii <- c(crossprod(u, Viu) / utu)
        ViX <- lag.listw(Vi, X, zero.policy=TRUE)
        MViM <- t(X) %*% ViX %*% XtXinv
        t1 <- -sum(diag(MViM))
        sumsq.Vi <- function(x) {
            if (is.null(x)) NA
	    else sum(x^2)
	}
	trVi2 <- sum(sapply(Vi$weights, sumsq.Vi), na.rm=TRUE)
	t2a <- sum(diag(t(ViX) %*% ViX %*% XtXinv))
	t2b <- sum(diag(MViM %*% MViM))
	t2 <- trVi2 - 2*t2a + t2b
	e1 <- 0.5 * (t1 + sqrt(2*t2 - t1^2))
	en <- 0.5 * (t1 - sqrt(2*t2 - t1^2))
        gamma_1n <- c(c(en), c(e1))
        obj <- exactMoran(Ii, gamma_1n, alternative=alternative,
            type="Local", np2=n-(2+p), useTP=useTP, truncErr=truncErr,
            zeroTreat=zeroTreat)
        data.name <- paste("region:", select[i],
	    attr(nb, "region.id")[select[i]],
	    "\n", paste(strwrap(paste("model: ", gsub("[ ]+", " ", 
	    paste(deparse(model$call), sep="", collapse="")))),
	    collapse="\n"),
            "\nneighbours:", deparse(substitute(nb)),
	    "style:", style, "\n")
        obj$data.name <- data.name
        obj$df <- (n-p)
        obj$i <- paste(select[i], attr(nb, "region.id")[select[i]])
        obj$Vi <- if(save.Vi) Vi else NULL
	res[[i]] <- obj
    }
    class(res) <- "localmoranex"
    res
}


print.localmoranex <- function(x, ...) {
    extract <- function(x, i) {x[[i]]}
    regnames <- sapply(x, extract, 10)
    est <- sapply(x, extract, 3)
    sad <- sapply(x, extract, 1)
    pval <- sapply(x, extract, 2)
    oT <- sapply(x, extract, 7)
    res <- as.matrix(cbind(est, sad, pval))
    rownames(res) <- regnames
    colnames(res) <- c("Local Morans I", "Exact SD", "Pr. (exact)")
    print(res, ...)
    if (any(oT != "E")) warning(paste("Normal reported for:",
        paste(which(oT != "E"), collapse=", ")), call.=FALSE)
    invisible(res)
}

as.data.frame.localmoranex <- function(x, row.names=NULL, optional=FALSE, ...) {
    n <- length(x)
    if (n < 1) stop("x too short")
    res <- matrix(0, nrow=n, ncol=11)
    regnames <- NULL
    if (!is.null(row.names)) 
	if (length(row.names) == n) regnames <- row.names
    if (is.null(regnames))for (i in 1:n) regnames <- c(regnames, x[[i]]$i)
    for (i in 1:n) {
        tau <- x[[i]]$gamma
	df <- x[[i]]$df
        if (length(tau) == 2L) tau <- c(tau[1], rep(0, df-2), tau[2])
        max.I <- tau[1]
        min.I <- tau[length(tau)]
        E.I <- sum(tau)/df
        tau <- tau - E.I
        V.I <- (2*sum(tau^2)) / (df*(df+2))
        Z.I <- (x[[i]]$estimate - E.I) / sqrt(V.I)
	if (x[[i]]$alternative == "two.sided") 
	    P.I <- 2 * pnorm(abs(Z.I), lower.tail=FALSE)
        else if (x[[i]]$alternative == "greater")
            P.I <- pnorm(Z.I, lower.tail=FALSE)
        else P.I <- pnorm(Z.I)
        Sk.I <- ((8*sum(tau^3))/(df*(df+2)*(df+4))) / (V.I^(3/2))
        Kur.I <- ((48*sum(tau^4) + 12*(sum(tau^2))^2) /
            (df*(df+2)*(df+4)*(df+6))) / (V.I^2)
	res[i,] <- c(x[[i]]$estimate, Z.I, P.I, x[[i]]$statistic,
	    x[[i]]$p.value, E.I, V.I, Sk.I, Kur.I, min.I, max.I)
    }
    colnames(res) <- c("Local Morans I", "Stand. dev. (N)", "Pr. (N)",
        "Exact SD", "Pr. (exact)", "Expectation", "Variance",
        "Skewness", "Kurtosis", "Minimum", "Maximum")
    rownames(res) <- regnames
    res <- as.data.frame(res)
    extract <- function(x, i) {x[[i]]}
    res$oT <- sapply(x, extract, 7)
    res
}




