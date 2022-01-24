# Copyright (c) 2007-2022 Markus Reder and Roger Bivand

localmoran.exact.alt <- function(model, select, nb, glist = NULL, style = "W",
    zero.policy = NULL, alternative = "two.sided", spChk=NULL, 
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
    if (any(select < 1) || any(select > n))
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

    exactLocalMoranAlt_int <- function(i, B, select, style, n, D, a, 
        zero.policy, u, utu, M1, M2, useTP, truncErr, zeroTreat, alternative) {
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
	obj
    }
    if (parallel == "snow") {
      if (requireNamespace("parallel", quietly = TRUE)) {
        sI <- parallel::splitIndices(n, length(cl))
        env <- new.env()
        assign("select", select, envir=env)
        assign("B", B, envir=env)
        assign("style", style, envir=env)
        assign("n", n, envir=env)
        assign("D", D, envir=env)
        assign("a", a, envir=env)
        assign("zero.policy", zero.policy, envir=env)
        assign("alternative", alternative, envir=env)
        assign("u", u, envir=env)
        assign("utu", utu, envir=env)
        assign("M1", M2, envir=env)
        assign("M2", M2, envir=env)
        assign("useTP", useTP, envir=env)
        assign("truncErr", truncErr, envir=env)
        assign("zeroTreat", zeroTreat, envir=env)
        parallel::clusterExport(cl, varlist=c("select", "B", "style",
            "n", "D", "a", "zero.policy", "alternative", "u", "M1", "M2",
            "utu", "useTP", "truncErr", "zeroTreat"), envir=env)
        oo <- parallel::clusterApply(cl, x = sI, fun=lapply, function(i) {
            exactLocalMoranAlt_int(i, B, select, style, n, D, a, 
            zero.policy, u, utu, M1, M2, useTP, truncErr, zeroTreat,
            alternative)})
        res <- do.call("c", oo)
        rm(env)
      } else {
        stop("parallel not available")
      }
    } else if (parallel == "multicore") {
      if (requireNamespace("parallel", quietly = TRUE)) {
        sI <- parallel::splitIndices(n, ncpus)
        oo <- parallel::mclapply(sI, FUN=lapply, function(i) {
            exactLocalMoranAlt_int(i, B, select, style, n, D, a, 
            zero.policy, u, utu, M1, M2, useTP, truncErr, zeroTreat,
            alternative)}, mc.cores=ncpus)
        res <- do.call("c", oo)
      } else {
        stop("parallel not available")
      }
    } else {
        res <- lapply(1:n, function(i) exactLocalMoranAlt_int(i, B, select, 
            style, n, D, a, zero.policy, u, utu, M1, M2, useTP, truncErr, 
            zeroTreat, alternative))
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

