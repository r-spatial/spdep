# Copyright (c) 2007-2008 Markus Reder and Roger Bivand

localmoran.exact.alt <- function(model, select, nb, glist = NULL, style = "W",
    zero.policy = NULL, alternative = "greater", spChk=NULL, 
    resfun=weighted.residuals, Omega=NULL, save.Vi = FALSE, save.M=FALSE,
    useTP=FALSE, truncErr=1e-6, zeroTreat=0.1) {
# need to impose check on weights TODO!!
# class to inherits Jari Oksanen 080603
    if (!inherits(nb, "nb"))
        stop(paste(deparse(substitute(nb)), "not an nb object"))
        if (is.null(zero.policy))
            zero.policy <- get("zeroPolicy", envir = .spdepOptions)
        stopifnot(is.logical(zero.policy))
#    if (class(model) != "lm") 
#        stop(paste(deparse(substitute(model)), "not an lm object"))
    dmc <- deparse(model$call)
    n <- length(nb)
    if (!inherits(model, "lm"))
     	stop(paste(deparse(substitute(model)), "not an lm object"))
    if (is.null(Omega)) Omega <- diag(n)
    else {
        if (dim(Omega)[1] != n) stop("Omega of different size than data")
        Omega <- chol(Omega)
    }
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
    M <- diag(n) - X %*% tcrossprod(XtXinv, X)
    M1 <- Omega %*% M
    M2 <- M %*% t(Omega)
    B <- listw2U(nb2listw(nb, glist=glist, style="B",
	zero.policy=zero.policy))
    D <- NULL
    a <- NULL
    if (style == "W") {
        D <- 1/sapply(B$weights, sum)
    } else if (style == "S") {
        D <- 1 / sqrt(sapply(B$weights, function(x) sum(x^2)))
#        a <- sum(unlist(B$weights))
# correction by Danlin Yu, 25 March 2004
	a <- sum(sapply(B$weights, function(x) sqrt(sum(x^2))))
    } else if (style == "C") a <- sum(unlist(B$weights))
    res <- vector(mode="list", length=length(select))
    for (i in 1:length(select)) {
        Vi <- listw2star(B, select[i], style=style, n, D, a,
	    zero.policy=zero.policy)
        Viu <- lag.listw(Vi, u, zero.policy=TRUE)
	Ii <- c(crossprod(u, Viu) / utu)

        obj <- exactLocalMoranAlt(Ii=Ii, Vi=Vi, M1=M1, M2=M2, n=n,
            alternative=alternative, useTP=useTP, truncErr=truncErr,
                zeroTreat=zeroTreat)
        data.name <- paste("region:", select[i],
	    attr(nb, "region.id")[select[i]],
	    "\n", paste(strwrap(paste("model: ", gsub("[ ]+", " ", 
	    paste(dmc, sep="", collapse="")))),
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
    if (save.M) attr(res, "M") <- list(M1=M1, M2=M2)
    res
}

exactLocalMoranAlt <- function(Ii, Vi, M1, M2, n, alternative,
    type="Alternative", useTP=FALSE, truncErr=1e-6, zeroTreat=0.1) {
    ViI <- listw2mat(Vi) - Ii * diag(n)
    innerTerm <- M1 %*% ViI %*% M2
    evalue <- eigen(innerTerm, only.values=TRUE)$values
    gamma <- c(evalue)
    obj <- exactMoran(Ii, gamma, alternative=alternative,
        type=type, useTP=useTP, truncErr=truncErr, zeroTreat=zeroTreat)
    obj
}

