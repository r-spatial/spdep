# Copyright 2002-2008 by Roger Bivand and Michael Tiefelsdorf
#

lm.morantest.sad <- function (model, listw, zero.policy = NULL, 
    alternative = "greater", spChk=NULL, resfun=weighted.residuals, 
    tol = .Machine$double.eps^0.5, maxiter = 1000, tol.bounds=0.0001,
    zero.tol=1.0e-7, Omega=NULL, save.M=NULL, save.U=NULL) 
{
    if (!inherits(listw, "listw")) 
        stop(paste(deparse(substitute(listw)), "is not a listw object"))
    if (!inherits(model, "lm")) 
        stop(paste(deparse(substitute(model)), "not an lm object"))
        if (is.null(zero.policy))
            zero.policy <- get("zeroPolicy", envir = .spdepOptions)
        stopifnot(is.logical(zero.policy))
    N <- length(listw$neighbours)
    u <- resfun(model)
    if (N != length(u)) 
        stop("objects of different length")
    if (is.null(spChk)) spChk <- get.spChkOption()
    if (spChk && !chkIDs(u, listw))
	stop("Check of data and weights ID integrity failed")
    if (!(alternative %in% c("greater", "less", "two.sided")))
	stop("alternative must be one of: \"greater\", \"less\", or \"two.sided\"")
    u <- as.vector(u)
    listw.U <- listw2U(listw)
    S0 <- sum(unlist(listw.U$weights))
    lu <- lag.listw(listw.U, u, zero.policy = zero.policy)
    Nnn <- N
    if (zero.policy) Nnn <- length(which(card(listw$neighbours) > 0L))
    I <- (Nnn/S0) * ((t(u) %*% lu)/(t(u) %*% u))
    I_save <- I
    if (!isTRUE(all.equal((Nnn/S0), 1))) I <- I * (S0/Nnn)
    p <- model$rank
    p1 <- 1:p
    nacoefs <- which(is.na(coefficients(model)))
    XtXinv <- chol2inv(model$qr$qr[p1, p1, drop = FALSE])
    X <- model.matrix(terms(model), model.frame(model))
# fixed after looking at TOWN dummy in Boston data
    if (length(nacoefs) > 0L) X <- X[,-nacoefs]
    if (!is.null(wts <- weights(model))) {
	X <- sqrt(diag(wts)) %*% X
    }
    M <- diag(N) - X %*% XtXinv %*% t(X)
    U <- listw2mat(listw.U)
    if (is.null(Omega)) {
        MVM <- M %*% U %*% M
        MVM <- 0.5 * (t(MVM) + MVM)
        evalue <- eigen(MVM, only.values=TRUE)$values
        idxpos <- which(abs(evalue) < zero.tol)
        if (length(idxpos) != p)
            warning("number of zero eigenvalues greater than number of variables")
        idxpos <- idxpos[1] - 1
	if (idxpos < 1) {
            warning("first eigenvalue index zero")
	    tau <- c()
	} else tau <- evalue[1:idxpos]
        tau <- c(tau, evalue[(idxpos+1+p):N])
        mres <- moranSad(tau, I, tol, maxiter, tol.bounds,
            alternative=alternative, type="Global")
    } else {
        if (dim(Omega)[1] != N) stop("Omega of different size than data")
        mres <- moranSadAlt(I, M, U, Omega, N, tol=tol, maxiter=maxiter,
            tol.bounds=tol.bounds, alternative=alternative)
    }
    statistic <- mres$sad.p
    attr(statistic, "names") <- "Saddlepoint approximation"
    p.value <- mres$p.sad
    estimate <- c(I_save)
    attr(estimate, "names") <- "Observed Moran I"
    internal1 <- c(mres$omega, mres$sad.r, mres$sad.u)
    attr(internal1, "names") <- c("omega", "sad.r", "sad.u")
    internal2 <- unlist(mres$root)[2:4]
    attr(internal2, "names") <- c("f.root", "iter", "estim.prec")
    method <- paste("Saddlepoint approximation for global Moran's I",
        "(Barndorff-Nielsen formula)")
    data.name <- paste("\nmodel:", paste(strwrap(gsub("[[:space:]]+", " ", 
	    paste(deparse(model$call), sep="", collapse=""))), collapse="\n"),
    	    "\nweights: ", deparse(substitute(listw)), "\n", sep="")
    res <- list(statistic = statistic, p.value = p.value,
        estimate = estimate, method = method,
	alternative = alternative, data.name = data.name,
	internal1 = internal1, internal2 = internal2,
	df = (N-p), tau = tau)
    class(res) <- "moransad"
    if (!is.null(save.M)) res$M <- M
    if (!is.null(save.U)) res$U <- U
    return(res)
}

moranSadAlt <- function(I, M, U, Omega, n, tol=.Machine$double.eps^0.5,
    maxiter=1000, tol.bounds=0.0001, alternative="greater") {
    Omega <- chol(Omega)
    A <- Omega %*% M %*% (U- diag(I, n)) %*% M %*% t(Omega)
    tau <- sort(eigen(A)$values, decreasing=TRUE)
    obj <- moranSad(tau, I, tol=tol, maxiter=maxiter,
            tol.bounds=tol.bounds, alternative=alternative, type="Alternative")
    obj
}

moranSad <- function(tau, I, tol=.Machine$double.eps^0.5, maxiter=1000,
    tol.bounds=0.0001, alternative="greater", type="Global") {
    if (type == "Global") taumi <- tau - c(I)
    else if (type == "Alternative") taumi <- tau
    low <- (1 / (2*taumi[length(taumi)])) + tol.bounds
    high <- (1 / (2*taumi[1])) - tol.bounds
    if (!(low < high)) {
        omega <- root <- p.sad <- sad.p <- sad.r <- sad.u <- NA
        warning("low not less than high for uniroot in moranSad")
    } else {
        f <- function(omega, taumi) {sum(taumi/(1 - (2*omega*taumi)))}
        root <- uniroot(f, lower=low, upper=high, tol=tol, maxiter=maxiter,
            taumi=taumi)
        omega <- root$root
        if (omega < 0 ) sad.r <- -sqrt(sum(log(1 - 2*omega*taumi)))
        else sad.r <- sqrt(sum(log(1 - 2*omega*taumi)))
        sad.u <- omega * sqrt(2*sum(taumi^2 / (1 - (2*omega*taumi))^2))
        sad.p <- sad.r - ((1/sad.r)*log(sad.r/sad.u))
        if (alternative == "two.sided") p.sad <- 2 * pnorm(abs(sad.p), 
	    lower.tail=FALSE)
        else if (alternative == "greater")
            p.sad <- pnorm(sad.p, lower.tail=FALSE)
        else p.sad <- pnorm(sad.p)
        if (!is.finite(p.sad) || p.sad < 0 || p.sad > 1) 
	    warning("Out-of-range p-value: reconsider test arguments")
    }
    return(list(p.sad=p.sad, sad.p=sad.p, sad.r=sad.r, sad.u=sad.u,
        omega=omega, root=root))
}

print.moransad <- function(x, ...) {
    class(x) <- c("htest", "moransad")
    print(x, ...)
    invisible(x)
}

summary.moransad <- function(object, ...) {
    res <- object
    tau <- object$tau
    df <- object$df
    if (length(tau) == 2) tau <- c(tau[1], rep(0, df-2), tau[2])
    max.I <- tau[1]
    min.I <- tau[length(tau)]
    E.I <- sum(tau)/df
    tau <- tau - E.I
    V.I <- (2*sum(tau^2)) / (df*(df+2))
    Z.I <- (object$estimate - E.I) / sqrt(V.I)
    Sk.I <- ((8*sum(tau^3))/(df*(df+2)*(df+4))) / (V.I^(3/2))
    Kur.I <- ((48*sum(tau^4) + 12*(sum(tau^2))^2) /
        (df*(df+2)*(df+4)*(df+6))) / (V.I^2)
    res$xtra <- c(E.I, V.I, Z.I, Sk.I, Kur.I, min.I, max.I)
    attr(res$xtra, "names") <- c("Expectation", "Variance",
        "Std. deviate", "Skewness", "Kurtosis", "Minimum", "Maximum")
    class(res) <- c("summary.moransad", "moransad")
    res
}

print.summary.moransad <- function(x, ...) {
    class(x) <- c("htest", "summary.moransad", "moransad")
    print(x, ...)
    print(c(x$xtra, x$internal1), ...)
    if (!is.null(x$internal2)) print(x$internal2, ...)
    invisible(x)
}

