# Copyright 2005-2012 by Roger Bivand
spautolm <- function(formula, data = list(), listw, weights,
    na.action, family="SAR", method="eigen", verbose=NULL, trs=NULL,
    interval=NULL, zero.policy=NULL, tol.solve=.Machine$double.eps, llprof=NULL,
    control=list()) {
    timings <- list()
    .ptime_start <- proc.time()
    con <- list(tol.opt=.Machine$double.eps^(2/3), 
        fdHess=NULL, optimHess=FALSE, optimHessMethod="optimHess",
        Imult=2, cheb_q=5, MC_p=16, MC_m=30, super=NULL, spamPivot="MMD",
        in_coef=0.1, type="MC",
        correct=TRUE, trunc=TRUE, SE_method="LU", nrho=200,
        interpn=2000, small_asy=TRUE, small=1500, SElndet=NULL,
        LU_order=FALSE, pre_eig=NULL)
    nmsC <- names(con)
    con[(namc <- names(control))] <- control
    if (length(noNms <- namc[!namc %in% nmsC])) 
        warning("unknown names in control: ", paste(noNms, collapse = ", "))
    if (!inherits(listw, "listw")) 
        stop("No neighbourhood list")
    if (is.null(verbose)) verbose <- get("verbose", envir = .spdepOptions)
    stopifnot(is.logical(verbose))
        if (is.null(zero.policy))
            zero.policy <- get("zeroPolicy", envir = .spdepOptions)
        stopifnot(is.logical(zero.policy))

    if (family == "SMA" && method != "eigen") stop("SMA only for eigen method")
    if (method == "spam" || method == "spam_update") stop("spam not supported as method")
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "weights", "na.action"), names(mf), 0)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    mt <- attr(mf, "terms")

#    mt <- terms(formula, data = data)
#    mf <- lm(formula, data, , weights, na.action=na.action,
#        method="model.frame")
    na.act <- attr(mf, "na.action")
    if (!is.null(na.act)) {
        subset <- !(1:length(listw$neighbours) %in% na.act)
        listw <- subset(listw, subset, zero.policy=zero.policy)
    }

    Y <- model.extract(mf, "response")
    if (any(is.na(Y))) stop("NAs in dependent variable")
    X <- model.matrix(mt, mf)
    if (any(is.na(X))) stop("NAs in independent variable")
    n <- nrow(X)
    if (n != length(listw$neighbours))
	 stop("Input data and neighbourhood list have different dimensions")
    weights <- as.vector(model.extract(mf, "weights"))
# set up default weights
    if (!is.null(weights) && !is.numeric(weights)) 
        stop("'weights' must be a numeric vector")
    if (is.null(weights)) weights <- rep(as.numeric(1), n)
    if (any(is.na(weights))) stop("NAs in weights")
    if (any(weights < 0)) stop("negative weights")
    lm.base <- lm(Y ~ X - 1, weights=weights)
    aliased <- is.na(coefficients(lm.base))
    cn <- names(aliased)
    names(aliased) <- substr(cn, 2, nchar(cn))
    if (any(aliased)) {
        nacoef <- which(aliased)
# bug x for X Bjarke Christensen 090924
	X <- X[,-nacoef]
    }
    can.sim <- FALSE
    if (listw$style %in% c("W", "S")) 
	can.sim <- can.be.simmed(listw)

    sum_lw <- sum(log(weights))
#    env <- new.env(parent=globalenv())
    env <- new.env()
    assign("Y", Y, envir=env)
    assign("X", X, envir=env)
    assign("n", n, envir=env)
    assign("weights", weights, envir=env)
    assign("can.sim", can.sim, envir=env)
    assign("family", family, envir=env)
    assign("method", method, envir=env)
    assign("verbose", verbose, envir=env)
    assign("listw", listw, envir=env)
    assign("sum_lw", sum_lw, envir=env)
    W <- as(listw, "CsparseMatrix")
    if (family == "CAR") if (!isTRUE(all.equal(W, t(W))))
        warning("Non-symmetric spatial weights in CAR model")
    assign("W", W, envir=env)
    I <- as_dsCMatrix_I(n)
    assign("I", I, envir=env)
    Sweights <- as(as(Diagonal(x=weights), "symmetricMatrix"), 
        "CsparseMatrix")
    assign("Sweights", Sweights, envir=env)
    timings[["set_up"]] <- proc.time() - .ptime_start
    .ptime_start <- proc.time()

    if (verbose) cat(paste("\nJacobian calculated using "))

    interval <- jacobianSetup(method, env, con, pre_eig=con$pre_eig, trs=trs,
        interval=interval)
    assign("interval", interval, envir=env)

# fix SMA bounds
    if (family == "SMA") interval <- -rev(interval)

    nm <- paste(method, "set_up", sep="_")
    timings[[nm]] <- proc.time() - .ptime_start
    .ptime_start <- proc.time()

    if (!is.null(llprof)) {
        if (length(llprof) == 1L)
            llprof <- seq(interval[1], interval[2], length.out=llprof)
        ll_prof <- numeric(length(llprof))
        for (i in seq(along=llprof)) 
            ll_prof[i] <- .opt.fit(llprof[i], env=env, tol.solve=tol.solve)
        nm <- paste(method, "profile", sep="_")
        timings[[nm]] <- proc.time() - .ptime_start
        .ptime_start <- proc.time()
    }

    opt <- optimize(.opt.fit, interval=interval, maximum=TRUE,
        tol = con$tol.opt, env=env, tol.solve=tol.solve)
    lambda <- opt$maximum
    if (isTRUE(all.equal(lambda, interval[1])) ||
        isTRUE(all.equal(lambda, interval[2]))) 
        warning("lambda on interval bound - results should not be used")
    names(lambda) <- "lambda"
    LL <- opt$objective
    nm <- paste(method, "opt", sep="_")
    timings[[nm]] <- proc.time() - .ptime_start
    .ptime_start <- proc.time()

# get GLS coefficients
    fit <- .SPAR.fit(lambda=lambda, env, out=TRUE, tol.solve=tol.solve)
# create residuals and fitted values (Cressie 1993, p. 564)
    fit$signal_trend <- drop(X %*% fit$coefficients)
    fit$signal_stochastic <- drop(lambda * W %*% (Y - fit$signal_trend))
    fit$fitted.values <- fit$signal_trend + fit$signal_stochastic
    fit$residuals <- drop(Y - fit$fitted.values)

# get null LL
    LL0 <- .opt.fit(lambda=0, env, tol.solve=tol.solve)
# NK null
    LLNullLlm <- logLik(lm(Y ~ 1, weights=weights))
    nm <- paste(method, "output", sep="_")
    timings[[nm]] <- proc.time() - .ptime_start
    .ptime_start <- proc.time()
#    if (method != "eigen") {
#        if (con$small >= n && con$small_asy) do_asy <- TRUE
#        else do_asy <- FALSE
#    } else do_asy <- TRUE
    do_asy <- FALSE
    if (is.null(con$fdHess)) {
        con$fdHess <-  !do_asy #&& method != "eigen"
        fdHess <- NULL
    }
    stopifnot(is.logical(con$fdHess))
    lambda.se <- NULL

    if (con$fdHess) {
        coefs <- c(lambda, fit$coefficients)
        fdHess <- getVcovmat(coefs, env, tol.solve=tol.solve,
            optim=con$optimHess, optimM=con$optimHessMethod)
        lambda.se <- sqrt(fdHess[1, 1])
    }

    timings[["fdHess"]] <- proc.time() - .ptime_start
    rm(env)
    GC <- gc()
    res <- list(fit=fit, lambda=lambda, LL=LL, LL0=LL0, call=match.call(),
        parameters=(ncol(X)+2), aliased=aliased, method=method, family=family,
        zero.policy=zero.policy, weights=weights, interval=interval, trs=trs,
        timings=do.call("rbind", timings)[, c(1, 3)], LLNullLlm=LLNullLlm,
        fdHess=fdHess, lambda.se=lambda.se, X=X, Y=Y)
    if (!is.null(na.act))
	res$na.action <- na.act
    if (is.null(llprof)) res$llprof <- llprof
    else {
        res$llprof <- list(lambda=llprof, ll=ll_prof)
    }
    if (zero.policy) {
        zero.regs <- attr(listw$neighbours, 
	    "region.id")[which(card(listw$neighbours) == 0)]
	if (length(zero.regs) > 0L)
	    attr(res, "zero.regs") <- zero.regs
	}

    class(res) <- "spautolm"
    res
}

.opt.fit <- function(lambda, env, tol.solve=.Machine$double.eps) {
# fitting function called from optimize()
    SSE <- .SPAR.fit(lambda=lambda, env=env, out=FALSE, tol.solve=tol.solve)
    n <- get("n", envir=env)
    s2 <- SSE/n
    ldet <- do_ldet(lambda, env)
    det <- ifelse(get("family", envir=env) == "CAR", 0.5*ldet, ldet)
    ret <- (det + (1/2)*get("sum_lw", envir=env) - ((n/2)*log(2*pi)) - 
        (n/2)*log(s2) - (1/(2*(s2)))*SSE)
    if (get("verbose", envir=env))  cat("lambda:", lambda, "function:", ret, "Jacobian", ldet, "SSE", SSE, "\n")
    ret
}


.SPAR.fit <- function(lambda, env, out=FALSE, tol.solve=.Machine$double.eps) {
    dmmf <- eval(parse(text=get("family", envir=env)))
    if (get("family", envir=env) == "SMA") IlW <- dmmf((get("I", envir=env) + 
        lambda * get("W", envir=env)), get("Sweights", envir=env))
    else IlW <- dmmf((get("I", envir=env) - lambda * get("W", envir=env)), 
        get("Sweights", envir=env))
    X <- get("X", envir=env)
    Y <- get("Y", envir=env)
    imat <- base::solve(crossprod(X, as.matrix(IlW %*% X)), tol=tol.solve)
    coef <- crossprod(imat, crossprod(X, as.matrix(IlW %*% Y)))
    fitted <- X %*% coef
    residuals <- Y - fitted
    SSE <- c(crossprod(residuals, as.matrix(IlW %*% residuals)))
    if (!out) return(SSE)

    n <- get("n", envir=env)
    s2 <- SSE/n
#    var <- s2 * diag(imat)
    coef <- c(coef)
    names(coef) <- colnames(X)
    res <- list(coefficients=coef, SSE=c(SSE), s2=c(s2), imat=imat,
        N=length(residuals))
    res
}

# Simultaneous autoregressive
SAR <- function(IlW, weights) {
    t(IlW) %*% weights %*% IlW
}

# Conditional  autoregressive
CAR <- function(IlW, weights) {
    IlW %*% weights
}

# Spatial moving average
SMA <- function(IlW, weights) {
    IlW <- solve(IlW)
    t(IlW) %*% weights %*% IlW
}


print.spautolm <- function(x, ...) {
        if (isTRUE(all.equal(x$lambda, x$interval[1])) ||
            isTRUE(all.equal(x$lambda, x$interval[2]))) 
            warning("lambda on interval bound - results should not be used")
	cat("\nCall:\n")
	print(x$call)
	cat("\nCoefficients:\n")
	print(coef(x))
	cat("\nLog likelihood:", logLik(x), "\n")
	invisible(x)
    
}

residuals.spautolm <- function(object, ...) {
	if (is.null(object$na.action))
		object$fit$residuals
	else napredict(object$na.action, object$fit$residuals)
}

fitted.spautolm <- function(object, ...) {
	if (is.null(object$na.action))
		object$fit$fitted.values
	else napredict(object$na.action, object$fit$fitted.values)
}

deviance.spautolm <- function(object, ...) {
	object$SSE
}

coef.spautolm <- function(object, ...) {
	c(object$fit$coefficients, object$lambda)
}


logLik.spautolm <- function(object, ...) {
	LL <- c(object$LL)
	class(LL) <- "logLik"
	N <- object$fit$N
	attr(LL, "nall") <- N
	attr(LL, "nobs") <- N
	attr(LL, "df") <- object$parameters
	LL
}

LR1.spautolm <- function(object)
{
	if (!inherits(object, "spautolm")) stop("Not a spautolm object")
	LLx <- logLik(object)
	LLy <- object$LL0
	statistic <- 2*(LLx - LLy)
	attr(statistic, "names") <- "Likelihood ratio"
	parameter <- 1
	attr(parameter, "names") <- "df"
	p.value <- 1 - pchisq(abs(statistic), parameter)
	estimate <- c(LLx, LLy)
	attr(estimate, "names") <- c(paste("Log likelihood of spatial regression fit"), paste("Log likelihood of OLS fit",
		deparse(substitute(y))))
	method <- "Likelihood Ratio diagnostics for spatial dependence"
	res <- list(statistic=statistic, parameter=parameter,
		p.value=p.value, estimate=estimate, method=method)
	class(res) <- "htest"
	res
}

summary.spautolm <- function(object, correlation = FALSE, adj.se=FALSE,
 Nagelkerke=FALSE, ...) {
	N <- object$fit$N
	adj <- ifelse (adj.se, N/(N-length(object$fit$coefficients)), 1) 
	object$fit$s2 <- object$fit$s2*adj
	object$resvar <- object$fit$s2*object$fit$imat
	rownames(object$resvar) <- colnames(object$resvar) <- 
		names(object$fit$coefficients)
	object$adj.se <- adj.se

	object$rest.se <- sqrt(diag(object$resvar))
	object$Coef <- cbind(object$fit$coefficients, object$rest.se, 
		object$fit$coefficients/object$rest.se,
		2*(1-pnorm(abs(object$fit$coefficients/object$rest.se))))
	colnames(object$Coef) <- c("Estimate", "Std. Error", 
		ifelse(adj.se, "t value", "z value"), "Pr(>|z|)")
        if (Nagelkerke) {
            nk <- NK.sarlm(object)
            if (!is.null(nk)) object$NK <- nk
        }
	if (correlation) {
		object$correlation <- diag((diag(object$resvar))
			^(-1/2)) %*% object$resvar %*% 
			diag((diag(object$resvar))^(-1/2))
		dimnames(object$correlation) <- dimnames(object$resvar)
	}
	object$LR1 <- LR1.spautolm(object)
	rownames(object$Coef) <- names(object$fit$coefficients)
	structure(object, class=c("summary.spautolm", class(object)))
}

print.summary.spautolm <- function(x, digits = max(5, .Options$digits - 3),
	signif.stars = FALSE, ...)
{
	cat("\nCall: ", deparse(x$call),	sep = "", fill=TRUE)
        if (isTRUE(all.equal(x$lambda, x$interval[1])) ||
            isTRUE(all.equal(x$lambda, x$interval[2]))) 
            warning("lambda on interval bound - results should not be used")
	cat("\nResiduals:\n")
	resid <- residuals(x)
	nam <- c("Min", "1Q", "Median", "3Q", "Max")
	rq <- if (length(dim(resid)) == 2L) 
		structure(apply(t(resid), 1, quantile), dimnames = list(nam, 
			dimnames(resid)[[2]]))
	else structure(quantile(resid), names = nam)
	print(rq, digits = digits, ...)
	if (x$zero.policy) {
		zero.regs <- attr(x, "zero.regs")
		if (!is.null(zero.regs))
			cat("\nRegions with no neighbours included:\n",
			zero.regs, "\n")
	}
	cat("\nCoefficients:", x$coeftitle, "\n")
	coefs <- x$Coef
	if (!is.null(aliased <- x$aliased) && any(x$aliased)){
		cat("    (", table(aliased)["TRUE"], 
			" not defined because of singularities)\n", sep = "")
		cn <- names(aliased)
		coefs <- matrix(NA, length(aliased), 4, dimnames = list(cn, 
                	colnames(x$Coef)))
            	coefs[!aliased, ] <- x$Coef
	}
	printCoefmat(coefs, signif.stars=signif.stars, digits=digits,
		na.print="NA")
	res <- x$LR1
	cat("\nLambda:", format(signif(x$lambda, digits)),
		"LR test value:", format(signif(res$statistic, digits)),
		"p-value:", format.pval(res$p.value, digits), 
		"\n")
        if (!is.null(x$lambda.se))
            cat("Numerical Hessian standard error of lambda:",
                format(signif(x$lambda.se, digits)), "\n")
	cat("\nLog likelihood:", logLik(x), "\n")
	if (x$adj.se) cat("Residual variance (sigma squared): ") 
	else cat("ML residual variance (sigma squared): ") 
	cat(format(signif(x$fit$s2, digits)), ", (sigma: ", 
		format(signif(sqrt(x$fit$s2), digits)), ")\n", sep="")
	cat("Number of observations:", x$fit$N, "\n")
	cat("Number of parameters estimated:", x$parameters, "\n")
	cat("AIC: ", format(signif(AIC(x), digits)), "\n", sep="")
        if (!is.null(x$NK)) cat("Nagelkerke pseudo-R-squared:",
            format(signif(x$NK, digits)), "\n")
    	correl <- x$correlation
    	if (!is.null(correl)) {
        	p <- NCOL(correl)
        	if (p > 1) {
            		cat("\nCorrelation of Coefficients:\n")
                	correl <- format(round(correl, 2), nsmall = 2, 
                  	digits = digits)
                	correl[!lower.tri(correl)] <- ""
                	print(correl[-1, -p, drop = FALSE], quote = FALSE)
            	}
    	}
    	cat("\n")
        invisible(x)
}

getVcovmat <- function(coefs, env, tol.solve=.Machine$double.eps, optim=FALSE,
    optimM="optimHess") {
    if (optim) {
      if (optimM == "nlm") {
           options(warn=-1)
           opt <- nlm(f=f_spautolm_hess_nlm, p=coefs, env=env, hessian=TRUE)
           options(warn=0)
           mat <- opt$hessian
#        opt <- optimHess(par=coefs, fn=f_laglm_hess, env=env)
#        mat <- opt
       } else if (optimM == "optimHess") {
           mat <- optimHess(par=coefs, fn=f_spautolm_hess, env=env)
       } else {
           opt <- optim(par=coefs, fn=f_spautolm_hess, env=env, method=optimM,
           hessian=TRUE)
           mat <- opt$hessian
      }
#        opt <- optimHess(par=coefs, fn=f_spautolm_hess, env=env)
#        mat <- opt
    } else {
        fd <- fdHess(coefs, f_spautolm_hess, env)
        mat <- fd$Hessian
    }
    res <- solve(-(mat), tol.solve=tol.solve)
    res
}

f_spautolm_hess_nlm <- function(coefs, env) {
    ret <- f_spautolm_hess(coefs, env)
    -ret
}

f_spautolm_hess <- function(coefs, env) {
    lambda <- coefs[1]
    int <- get("interval", envir=env)
    if (lambda <= int[1] || lambda >= int[2]) return(-Inf)
    beta <- coefs[-1]
    X <- get("X", envir=env)
    Y <- get("Y", envir=env)
    fitted <- X %*% beta
    residuals <- Y - fitted
    dmmf <- eval(parse(text=get("family", envir=env)))
    if (get("family", envir=env) == "SMA") IlW <- dmmf((get("I", envir=env) + 
        lambda * get("W", envir=env)), get("Sweights", envir=env))
    else IlW <- dmmf((get("I", envir=env) - lambda * get("W", envir=env)), 
        get("Sweights", envir=env))
    SSE <- c(crossprod(residuals, as.matrix(IlW %*% residuals)))
    n <- get("n", envir=env)
    s2 <- SSE/n
    ldet <- do_ldet(lambda, env)
    det <- ifelse(get("family", envir=env) == "CAR", 0.5*ldet, ldet)
    ret <- (det + (1/2)*get("sum_lw", envir=env) - ((n/2)*log(2*pi)) - 
        (n/2)*log(s2) - (1/(2*(s2)))*SSE)
    if (get("verbose", envir=env))
        cat("lambda:", lambda, "function:", ret, "Jacobian", ldet, "SSE",
            SSE, "\n")
    if (!is.finite(ret)) return(-Inf)
    ret
}





