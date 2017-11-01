# Copyright 2010-2017 by Roger Bivand and Eric Blankmeyer

lagmess <- function(formula, data = list(), listw, zero.policy=NULL,
    na.action=na.fail, q=10, start=-2.5, control=list(), method="BFGS",
    verbose=NULL, use_expm=FALSE) {
    stopifnot(inherits(listw, "listw"))
    if (is.null(verbose)) verbose <- get("verbose", envir = .spdepOptions)
    stopifnot(is.logical(verbose))
        if (is.null(zero.policy))
            zero.policy <- get("zeroPolicy", envir = .spdepOptions)
        stopifnot(is.logical(zero.policy))
    if (listw$style != "W") warning("weights should be row-stochastic")
    if (class(formula) != "formula") formula <- as.formula(formula)
    mt <- terms(formula, data = data)
    mf <- lm(formula, data, na.action=na.action, method="model.frame")
    na.act <- attr(mf, "na.action")
    if (!is.null(na.act)) {
        subset <- !(1:length(listw$neighbours) %in% na.act)
        listw <- subset(listw, subset, zero.policy=zero.policy)
    }
    y <- model.extract(mf, "response")
    stopifnot(all(is.finite(y)))
    X <- model.matrix(mt, mf)
    stopifnot(all(is.finite(X)))
    env <- new.env()
    assign("y", y, envir=env)
    assign("X", X, envir=env)
    assign("n", length(y), envir=env)
    nullLL <- logLik(lm(formula, data, na.action=na.action))

    W <- as(listw, "CsparseMatrix")
    assign("W", W, envir=env)

    if (!use_expm) {
        Y <- powerWeightsMESS(W, y, q=q)
        assign("Y", Y, envir=env)
        v <- 0:(q-1)
        assign("v", v, envir=env)
        G1 <- diag(1/factorial(v))
        assign("G1", G1, envir=env)
    }

    if (use_expm) {
        bestmess <- optim(start, mymess1, gr=NULL, method=method, hessian=TRUE,
            control=control, env)
    } else {
        bestmess <- optim(start, mymess, gr=NULL, method=method, hessian=TRUE,
            control=control, env)
    }
    alpha <- bestmess$par[1]
    alphase <- 1.0/(bestmess$hessian[1,1])^0.5
    rho <- 1.0 - exp(alpha[1])

    if (use_expm) {
        Sy <- expAtv(alpha*W, y)$eAtv
    } else {
        va <- alpha^v
        Sy <- Y %*% G1 %*% va
    }
    data$Sy <- Sy
#    formula[[2]] <- formula(~ Sy)[[2]]
    lmobj <- lm(formula=update(formula, Sy ~ .), data=data)
    coefs <- c(alpha=alpha, coef(lmobj))
    if (use_expm) {
        mat <- optimHess(par=coefs, fn=mymess1_hess, env=env)
    } else {
        mat <- optimHess(par=coefs, fn=mymess_hess, env=env)
    }
    mess_hess <- solve(-(mat))

    call <- match.call()
    lmobj$call <- call

    res <- list(lmobj=lmobj, alpha=alpha, alphase=alphase, rho=rho,
        bestmess=bestmess, q=q, start=start, na.action=na.act,
        nullLL=nullLL, use_expm=use_expm, mess_hess=mess_hess)
    class(res) <- "lagmess"
    res
}

powerWeightsMESS <- function(W, y, q=10) {
        n <- dim(W)[1]
        res <- matrix(NA, nrow=n, ncol=q)
        res[,1] <- y
        last <- W %*% y
        res[,2] <- c(last[,1])
        for (i in 3:q) {
            last <- W %*% last
            res[,i] <- c(last[,1])
        }
        res
}

mymess <- function(alpha, env, verbose=FALSE) {
    va <- alpha^get("v", envir=env)
    Sy <- get("Y", envir=env) %*% get("G1", envir=env) %*% va
    lmobj <- lm(Sy ~ get("X", envir=env) - 1)
    res <- -c(logLik(lmobj))
    if (verbose) cat("res:", res, "\n")
    res
}

mymess_hess <- function(coefs, env) {
    alpha <- coefs[1]
    beta <- coefs[-1]
    va <- alpha^get("v", envir=env)
    Sy <- get("Y", envir=env) %*% get("G1", envir=env) %*% va
    res <- Sy - get("X", envir=env) %*% beta
    SSE <- c(crossprod(res))
    n <- get("n", envir=env)
    s2 <- SSE/n
    ret <- (0 - ((n/2) * log(2 * pi)) - (n/2) * log(s2) - 
        (1/(2 * s2)) * SSE)
    ret
}

mymess1_hess <- function(coefs, env) {
    alpha <- coefs[1]
    beta <- coefs[-1]
    Sy <- expAtv(alpha*get("W", envir=env), get("y", envir=env))$eAtv
    res <- Sy - get("X", envir=env) %*% beta
    SSE <- c(crossprod(res))
    n <- get("n", envir=env)
    s2 <- SSE/n
    ret <- (0 - ((n/2) * log(2 * pi)) - (n/2) * log(s2) - 
        (1/(2 * s2)) * SSE)
    ret
}

mymess1 <- function(alpha, env, verbose=FALSE) {
    Sy <- expAtv(alpha*get("W", envir=env), get("y", envir=env))$eAtv
    lmobj <- lm(Sy ~ get("X", envir=env) - 1)
    res <- -c(logLik(lmobj))
    if (verbose) cat("res:", res, "\n")
    res
}

print.lagmess <- function(x, ...) {
    print(x$lmobj, ...)
    cat("Alpha: ", x$alpha, "\n", sep="")
    invisible(x)
}

print.summary.lagmess <- function(x, digits = max(5, .Options$digits - 3),
    signif.stars = FALSE, ...) {
    cat("Matrix exponential spatial lag model:\n")
    if (x$use_expm) cat("(calculated with expm)\n")
    print(x$lmsum, signif.stars=signif.stars, digits=digits)
    cat("Alpha: ", format(signif(x$alpha, digits)), ", standard error: ",
        format(signif(x$alphase, digits)), "\n    z-value: ", 
        format(signif((x$alpha/x$alphase), digits)), ", p-value: ",
        format.pval(2 * (1 - pnorm(abs(x$alpha/x$alphase))), digits),
        "\n", sep="")
    res <- x$LR
    cat("LR test value: ", format(signif(res$statistic, digits)), 
        ", p-value: ", format.pval(res$p.value, digits), "\n", sep="")
    cat("Implied rho:", x$rho, "\n")
    cat("\n")
    invisible(x)
}

summary.lagmess <- function(object, ...) {
    object$lmsum <- summary(object$lmobj, ...)
    object$LR <- LR1.lagmess(object)
    class(object) <- "summary.lagmess"
    object
}

LR1.lagmess <- function(object) {
    LLx <- logLik(object)
    LLy <- object$nullLL
    statistic <- 2*(LLx - LLy)
    attr(statistic, "names") <- "Likelihood ratio"
    parameter <- abs(attr(LLx, "df") - attr(LLy, "df"))
    if (parameter < 1) 
	stop("non-positive degrees of freedom: no test possible")
    attr(parameter, "names") <- "df"
    p.value <- 1 - pchisq(abs(statistic), parameter)
    estimate <- c(LLx, LLy)
    attr(estimate, "names") <- c("Log likelihood of MESS fit",
        "Log likelihood of OLS fit")
    method <- "Likelihood Ratio diagnostics for spatial dependence"
    res <- list(statistic=statistic, parameter=parameter,
	p.value=p.value, estimate=estimate, method=method)
    class(res) <- "htest"
    res
}

residuals.lagmess <- function(object, ...) {
    object$lmobj$residuals
}

deviance.lagmess <- function(object, ...) {
    deviance(object$lmobj)
}

coef.lagmess <- function(object, ...) {
    ret <- NULL
    ap <- object$alpha
    names(ap) <- "alpha"
    ret <- c(ret, ap)
    ret <- c(ret, coef(object$lmobj))
    ret
}

fitted.lagmess <- function(object, ...) {
    object$lmobj$fitted.values
}

logLik.lagmess <- function (object, ...) 
{
    LL <- c(logLik(object$lmobj))
    class(LL) <- "logLik"
    N <- length(residuals(object))
    attr(LL, "nall") <- N
    attr(LL, "nobs") <- N
    attr(LL, "df") <- object$lmobj$rank + 2
    LL
}

#    res <- list(lmobj=lmobj, alpha=alpha, alphase=alphase, rho=rho, bestmess=bestmess, q=q, start=start, na.action=na.act, nullLL=nullLL, use_expm=use_expm, mess_hess=mess_hess)


impacts.lagmess <- function(obj, ..., R=NULL, listw=NULL, 
  tol=1e-6, empirical=FALSE) {
    if (!is.null(R)) stopifnot(!is.null(obj$mess_hess))
    stopifnot(!is.null(listw))
    timings <- list()
    .ptime_start <- proc.time()
    type <- class(obj)
    W <- as(listw, "CsparseMatrix")
    alpha <- obj$alpha
    S_W <- expm(-alpha*W)
    beta <- obj$lmobj$coefficients
    n <- length(obj$lmobj$residuals)
    icept <- grep("(Intercept)", names(beta))
    iicept <- length(icept) > 0L
    if (iicept) {
        P <- matrix(beta[-icept], ncol=1)
        bnames <- names(beta[-icept])
    } else {
        P <- matrix(beta, ncol=1)
        bnames <- names(beta)
    }
    res <- lagImpactsExact(S_W, P, n)
    timings[["expm_impacts"]] <- proc.time() - .ptime_start
    if (!is.null(R)) {
        .ptime_start <- proc.time()
        ialpha <- 1
        drop2beta <- 1
        samples <- mvrnorm(n=R, mu=c(alpha, beta), Sigma=obj$mess_hess,
            tol=tol, empirical=empirical)
        timings[["impacts_samples"]] <- proc.time() - .ptime_start
        .ptime_start <- proc.time()
        sres <- apply(samples, 1, processMessSample,
            drop2beta=drop2beta, type=type, iicept=iicept,
            icept=icept, n=n, W=W, ialpha=ialpha)
        timings[["process_samples"]] <- proc.time() - .ptime_start
        .ptime_start <- proc.time()
        if (length(bnames) == 1L) {
            direct <- as.mcmc(t(matrix(sapply(sres, function(x) x$direct),
                nrow=1)))
            indirect <- as.mcmc(t(matrix(sapply(sres,
                function(x) x$indirect), nrow=1)))
            total <- as.mcmc(t(matrix(sapply(sres, function(x) x$total),
                nrow=1)))
        } else {
            direct <- as.mcmc(t(sapply(sres, function(x) x$direct)))
            indirect <- as.mcmc(t(sapply(sres, function(x) x$indirect)))
            total <- as.mcmc(t(sapply(sres, function(x) x$total)))
        }
        colnames(direct) <- bnames
        colnames(indirect) <- bnames
        colnames(total) <- bnames
        timings[["postprocess_samples"]] <- proc.time() - .ptime_start
        res <- list(res=res, sres=list(direct=direct,
            indirect=indirect, total=total))
    }
    if (!is.null(R)) attr(res, "samples") <- list(samples=samples, irho=ialpha,
        drop2beta=drop2beta)
    attr(res, "timings") <- do.call("rbind", timings)[, c(1,3)]
    attr(res, "method") <- "exact"
    attr(res, "type") <- type
    attr(res, "bnames") <- bnames
    attr(res, "haveQ") <- FALSE
    class(res) <- "lagImpact"
    attr(res, "iClass") <- class(obj)
    res
}

processMessSample <- function(x, drop2beta, type, iicept, icept, n, W,
    ialpha) {
    alpha <- x[ialpha]
    S_W <- expm(-alpha*W)
    beta <- x[-drop2beta]
    if (iicept) {
      P <- matrix(beta[-icept], ncol=1)
    } else {
      P <- matrix(beta, ncol=1)
    }
    lagImpactsExact(S_W, P, n)
}
