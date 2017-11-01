# Copyright 1998-2016 by Roger Bivand (non-W styles Rein Halbersma)
#

errorsarlm <- function(formula, data = list(), listw, na.action, weights=NULL,
        etype="error", method="eigen", quiet=NULL, zero.policy=NULL,
        interval=NULL, tol.solve=1.0e-10, trs=NULL, control=list()) {
        timings <- list()
        .ptime_start <- proc.time()
        con <- list(tol.opt=.Machine$double.eps^0.5, returnHcov=TRUE,
            pWOrder=250, fdHess=NULL, optimHess=FALSE,
            optimHessMethod="optimHess", LAPACK=FALSE,
            compiled_sse=FALSE, Imult=2, cheb_q=5, MC_p=16, MC_m=30,
            super=NULL, spamPivot="MMD", in_coef=0.1, type="MC",
            correct=TRUE, trunc=TRUE, SE_method="LU", nrho=200,
            interpn=2000, small_asy=TRUE, small=1500, SElndet=NULL,
            LU_order=FALSE, pre_eig=NULL)
        nmsC <- names(con)
        con[(namc <- names(control))] <- control
        if (length(noNms <- namc[!namc %in% nmsC])) 
            warning("unknown names in control: ", paste(noNms, collapse = ", "))
        if (is.null(quiet)) quiet <- !get("verbose", envir = .spdepOptions)
        stopifnot(is.logical(quiet))
        if (is.null(zero.policy))
            zero.policy <- get("zeroPolicy", envir = .spdepOptions)
        stopifnot(is.logical(zero.policy))
	if (!inherits(listw, "listw")) stop("No neighbourhood list")
        stopifnot(is.logical(con$optimHess))
        stopifnot(is.logical(con$LAPACK))
#        stopifnot(is.logical(con$super))
        stopifnot(is.logical(con$compiled_sse))
        stopifnot(is.character(con$spamPivot))
        if (class(formula) != "formula") formula <- as.formula(formula)
#	mt <- terms(formula, data = data)
#	mf <- lm(formula, data, na.action=na.action, method="model.frame")
#
        mf <- match.call(expand.dots = FALSE)
        m <- match(c("formula", "data", "weights", "na.action"), names(mf), 0)
        mf <- mf[c(1, m)]
        mf$drop.unused.levels <- TRUE
        mf[[1]] <- as.name("model.frame")
        mf <- eval(mf, parent.frame())
        mt <- attr(mf, "terms")
#
	na.act <- attr(mf, "na.action")
	if (!is.null(na.act)) {
	    subset <- !(1:length(listw$neighbours) %in% na.act)
	    listw <- subset(listw, subset, zero.policy=zero.policy)
	}

	switch(etype, error = if (!quiet)
                cat("\nSpatial autoregressive error model\n"),
	    emixed = if (!quiet)
                cat("\nSpatial mixed autoregressive error model\n"),
	    stop("\nUnknown model type\n"))

	y <- model.response(mf, "numeric")
	if (any(is.na(y))) stop("NAs in dependent variable")
	x <- model.matrix(mt, mf)
	if (any(is.na(x))) stop("NAs in independent variable")
        n <- nrow(x)
	if (n != length(listw$neighbours))
	    stop("Input data and neighbourhood list have different dimensions")
#
        weights <- as.vector(model.extract(mf, "weights"))
        if (!is.null(weights) && !is.numeric(weights)) 
            stop("'weights' must be a numeric vector")
        if (is.null(weights)) weights <- rep(as.numeric(1), n)
        if (any(is.na(weights))) stop("NAs in weights")
        if (any(weights < 0)) stop("negative weights")
#
	m <- NCOL(x)
	xxcolnames <- colnames(x)

        stopifnot(is.logical(con$small_asy))
        if (method != "eigen") {
            if (con$small >= n && con$small_asy) do_asy <- TRUE
            else do_asy <- FALSE
        } else do_asy <- TRUE
        if (is.null(con$fdHess)) {
            con$fdHess <- method != "eigen" && !do_asy
            fdHess <- NULL
        }
        stopifnot(is.logical(con$fdHess))
	if (etype == "emixed") {
                WX <- create_WX(x, listw, zero.policy=zero.policy,
                    prefix="lag")
		x <- cbind(x, WX)
		m <- NCOL(x)
		rm(WX)
	}
# added aliased after trying boston with TOWN dummy
	lm.base <- lm(y ~ x - 1, weights=weights)
	aliased <- is.na(coefficients(lm.base))
	cn <- names(aliased)
	names(aliased) <- substr(cn, 2, nchar(cn))
	if (any(aliased)) {
		nacoef <- which(aliased)
		x <- x[,-nacoef]
	}
#
        sw <- sqrt(weights)
#
        LL_null_lm <- NULL
	if ("(Intercept)" %in% colnames(x)) LL_null_lm <- logLik(lm(y ~ 1))
	m <- NCOL(x)
	xcolnames <- colnames(x)
	K <- ifelse(xcolnames[1] == "(Intercept)", 2, 1)
	wy <- lag.listw(listw, y, zero.policy=zero.policy)
	if (any(is.na(wy)))
	    stop("NAs in lagged dependent variable")
# added no intercept Guillaume Blanchet 091103
	if (m > 1 || (m == 1 && K == 1)) {
	    WX <- matrix(nrow=n,ncol=(m-(K-1)))
	    for (k in K:m) {
		wx <- lag.listw(listw, x[,k], zero.policy=zero.policy)
		if (any(is.na(wx)))
		    stop("NAs in lagged independent variable")
		WX[,(k-(K-1))] <- wx
	    }
	}
	if (K == 2) {
# modified to meet other styles, email from Rein Halbersma
		wx1 <- as.double(rep(1, n))
		wx <- lag.listw(listw, wx1, zero.policy=zero.policy)
		if (m > 1) WX <- cbind(wx, WX)
		else WX <- matrix(wx, nrow=n, ncol=1)
	}
	colnames(WX) <- xcolnames
	rm(wx)

	can.sim <- FALSE
	if (listw$style %in% c("W", "S")) 
		can.sim <- can.be.simmed(listw)
        sum_lw <- sum(log(weights))

#        env <- new.env(parent=globalenv())
        env <- new.env()
        assign("y", y, envir=env)
        assign("x", x, envir=env)
        assign("wy", wy, envir=env)
        assign("WX", WX, envir=env)
        assign("n", n, envir=env)
        assign("p", m, envir=env)
        assign("verbose", !quiet, envir=env)
        assign("family", "SAR", envir=env)
        assign("compiled_sse", con$compiled_sse, envir=env)
        assign("first_time", TRUE, envir=env)
        assign("LAPACK", con$LAPACK, envir=env)
        assign("can.sim", can.sim, envir=env)
        assign("listw", listw, envir=env)
        assign("similar", FALSE, envir=env)
        assign("f_calls", 0L, envir=env)
        assign("hf_calls", 0L, envir=env)
        assign("sum_lw", sum_lw, envir=env)
        assign("sw", sw, envir=env)
        timings[["set_up"]] <- proc.time() - .ptime_start
        .ptime_start <- proc.time()

	if (!quiet) cat(paste("\nJacobian calculated using "))

        interval <- jacobianSetup(method, env, con, pre_eig=con$pre_eig,
            trs=trs, interval=interval)
        assign("interval", interval, envir=env)

        nm <- paste(method, "set_up", sep="_")
        timings[[nm]] <- proc.time() - .ptime_start
        .ptime_start <- proc.time()
        if (con$compiled_sse) {
             ptr <- .Call("opt_error_init", PACKAGE="spdep")
             assign("ptr", ptr, envir=env)
        }
	opt <- optimize(sar.error.f, interval=interval, 
		maximum=TRUE, tol=con$tol.opt, env=env)
	lambda <- opt$maximum
        if (isTRUE(all.equal(lambda, interval[1])) ||
            isTRUE(all.equal(lambda, interval[2]))) 
            warning("lambda on interval bound - results should not be used")
	names(lambda) <- "lambda"
	LL <- opt$objective
        if (con$compiled_sse) {
             .Call("opt_error_free", get("ptr", envir=env), PACKAGE="spdep")
        }
        nm <- paste(method, "opt", sep="_")
        timings[[nm]] <- proc.time() - .ptime_start
        .ptime_start <- proc.time()
	lm.target <- lm(I(y - lambda*wy) ~ I(x - lambda*WX) - 1,
            weights=weights)
	r <- as.vector(residuals(lm.target))
	fit <- as.vector(y - r)
	p <- lm.target$rank
	SSE <- deviance(lm.target)
	s2 <- SSE/n
	rest.se <- (summary(lm.target)$coefficients[,2])*sqrt((n-p)/n)
	coef.lambda <- coefficients(lm.target)
	names(coef.lambda) <- xcolnames
        sum_lm_target <- summary.lm(lm.target, correlation = FALSE)
        emixedImps <- NULL
	if (etype == "emixed") {
            odd <- (m%/%2) > 0
            if (odd) {
                m2 <- (m-1)/2
            } else {
                m2 <- m/2
            }
            if (K == 1 && odd) {
                warning("model configuration issue: no total impacts")
            } else {
                cm <- matrix(0, ncol=m, nrow=m2)
                if (K == 2) {
                    if (odd) {
                        rownames(cm) <- xxcolnames[2:(m2+1)]
                    } else {
                        rownames(cm) <- xxcolnames[1:m2]
                    }
                    for (i in 1:m2) cm[i, c(i+1, i+(m2+1))] <- 1
# drop bug fix 2016-09-21 Philipp Hunziker
                    dirImps <- sum_lm_target$coefficients[2:(m2+1), 1:2, drop=FALSE]
                    rownames(dirImps) <- rownames(cm)
                    indirImps <- sum_lm_target$coefficients[(m2+2):m, 1:2, drop=FALSE]
                    rownames(indirImps) <- rownames(cm)
                } else {
                    rownames(cm) <- xxcolnames[1:m2]
                    for (i in 1:m2) cm[i, c(i, i+m2)] <- 1
                    dirImps <- sum_lm_target$coefficients[1:m2, 1:2, drop=FALSE]
                    rownames(dirImps) <- rownames(cm)
                    indirImps <- sum_lm_target$coefficients[(m2+1):m, 1:2, drop=FALSE]
                    rownames(indirImps) <- rownames(cm)
                }
            }
            totImps <- as.matrix(estimable(lm.target, cm)[, 1:2, drop=FALSE])
            emixedImps <- list(dirImps=dirImps, indirImps=indirImps,
                totImps=totImps)
        }
        Vs <- sum_lm_target$cov.unscaled
        tarX <- model.matrix(lm.target)
        tary <- model.response(model.frame(lm.target))
	lm.model <- lm(y ~ x - 1, weights=weights)
        logLik_lm.model <- logLik(lm.model)
        AIC_lm.model <- AIC(lm.model)
	ase <- FALSE
	lambda.se <- NULL
	LMtest <- NULL
	asyvar1 <- FALSE
        force_assign_eigen <- FALSE
        Hcov <- NULL
        pWinternal <- NULL
        timings[["coefs"]] <- proc.time() - .ptime_start
        .ptime_start <- proc.time()
        assign("first_time", TRUE, envir=env)
	if (method == "eigen" || do_asy) {
		tr <- function(A) sum(diag(A))
		W <- listw2mat(listw)
		A <- solve(diag(n) - lambda*W)
                if (con$returnHcov) {
                    pp <- lm.model$rank
                    p1 <- 1L:pp
                    R <- chol2inv(lm.model$qr$qr[p1, p1, drop = FALSE])
                    B <- tcrossprod(R, (sw*x)) %*% A
#                    A <- solve(diag(n) - lambda*t(W))
                    C <- A %*% (sw*x) %*% R
                    Hcov <- B %*% C
                    attr(Hcov, "method") <- method
                    timings[["eigen_hcov"]] <- proc.time() - .ptime_start
                    .ptime_start <- proc.time()
                }

		WA <- W %*% A
		asyvar <- matrix(0, nrow=2+p, ncol=2+p)
		asyvar[1,1] <- n / (2*(s2^2))
#		asyvar[1,1] <- n / (2*(s2))
		asyvar[2,1] <- asyvar[1,2] <- tr(WA) / s2
#		asyvar[2,1] <- asyvar[1,2] <- tr(WA)
		asyvar[2,2] <- tr(WA %*% WA) + tr(crossprod(WA))
# bug found 100224 German Muchnik Izon
#		asyvar[3:(p+2),3:(p+2)] <- s2*(t(x - lambda*WX) %*% 
                xl <- sw * (x - lambda*WX)
#		asyvar[3:(p+2),3:(p+2)] <- crossprod(xl)
		asyvar[3:(p+2),3:(p+2)] <- crossprod(xl)/s2
		asyvar1 <- try(solve(asyvar, tol=tol.solve), silent=TRUE)
                if (class(asyvar1) == "try-error") {
                    timings[["eigen_se"]] <- proc.time() - .ptime_start
                    .ptime_start <- proc.time()
                    con$fdHess <- TRUE
                    force_assign_eigen <- TRUE
                    warning(paste("inversion of asymptotic covariance",
                        "matrix failed for tol.solve =", tol.solve,
                        "\n", strsplit(attr(asyvar1, "condition")$message,
                            ":")[[1]][2], "- using numerical Hessian."))
                } else {
		    rownames(asyvar1) <- colnames(asyvar1) <- 
			c("sigma", "lambda", xcolnames)
		
		    lambda.se <- sqrt(asyvar1[2,2])
#		lambda.se <- sqrt(s2*asyvar1[2,2])
                    timings[["eigen_se"]] <- proc.time() - .ptime_start
                    .ptime_start <- proc.time()
		    ase <- TRUE
                }
	} else {
                if (con$returnHcov) {
                    pp <- lm.model$rank
                    p1 <- 1L:pp
                    R <- chol2inv(lm.model$qr$qr[p1, p1, drop = FALSE])
                    B <- tcrossprod(R, (sw*x))
                    W <- as(get("listw", envir=env), "CsparseMatrix")
                    B0 <- powerWeights(W=W, rho=lambda, order=con$pWOrder,
                        X=B, tol=tol.solve)
                    if (!is.null(attr(B0, "internal")) &&
                        !attr(B0, "internal")$conv)
                        pWinternal <- c(pWinternal, attr(B0, "internal"))
                    B1 <- as(B0, "matrix")
                    C <- (sw*x) %*% R
                    C0 <- powerWeights(W=t(W), rho=lambda, order=con$pWOrder,
                        X=C, tol=tol.solve)
                    if (!is.null(attr(C0, "internal")) &&
                        !attr(C0, "internal")$conv)
                        pWinternal <- c(pWinternal, attr(C0, "internal"))
                    C1 <- as(C0, "matrix")
                    Hcov <- B1 %*% C1
                    attr(Hcov, "method") <- method
                    timings[["sparse_hcov"]] <- proc.time() - .ptime_start
                    .ptime_start <- proc.time()
                }
        }
        if (con$fdHess) {
            coefs <- c(lambda, coef.lambda)
            if (con$compiled_sse) {
                ptr <- .Call("hess_error_init", PACKAGE="spdep")
                assign("ptr", ptr, envir=env)
            }
            fdHess <- getVmate(coefs, env, s2, trs, tol.solve=tol.solve,
                optim=con$optimHess, optimM=con$optimHessMethod)
            if (con$compiled_sse) {
                .Call("hess_error_free", get("ptr", envir=env),
                    PACKAGE="spdep")
            }
            if (is.null(trs)) {
                rownames(fdHess) <- colnames(fdHess) <- 
                    c("lambda", colnames(x))
                if (method != "eigen" || force_assign_eigen) {
                    lambda.se <- sqrt(fdHess[1, 1])
                }
            } else {
                rownames(fdHess) <- colnames(fdHess) <- 
                    c("sigma2", "lambda", colnames(x))
                if (method != "eigen" || force_assign_eigen) {
                    lambda.se <- sqrt(fdHess[2, 2])
                }
            }
            nm <- paste(method, "fdHess", sep="_")
            timings[[nm]] <- proc.time() - .ptime_start
        }
	call <- match.call()
        if (method=="SE_classic") {
            iC <- get("intern_classic", envir=env)
        } else iC <- NULL
	names(r) <- names(y)
	names(fit) <- names(y)
	ret <- structure(list(type="error", etype=etype, lambda=lambda,
		coefficients=coef.lambda, rest.se=rest.se, 
		LL=LL, s2=s2, SSE=SSE, parameters=(m+2), #lm.model=lm.model, 
                logLik_lm.model=logLik_lm.model, AIC_lm.model=AIC_lm.model,
                coef_lm.model=coef(lm.model),
                tarX=tarX, tary=tary, y=y, X=x,
		method=method, call=call, residuals=r, #lm.target=lm.target,
		opt=opt, fitted.values=fit, ase=ase, #formula=formula,
		se.fit=NULL, resvar=asyvar1, similar=get("similar", envir=env),
		lambda.se=lambda.se, LMtest=LMtest, zero.policy=zero.policy, 
		aliased=aliased, LLNullLlm=LL_null_lm, Hcov=Hcov, Vs=Vs,
                interval=interval, fdHess=fdHess,
                optimHess=con$optimHess, insert=!is.null(trs), trs=trs,
                timings=do.call("rbind", timings)[, c(1, 3)], 
                f_calls=get("f_calls", envir=env),
                hf_calls=get("hf_calls", envir=env), intern_classic=iC,
                pWinternal=pWinternal, weights=weights, emixedImps=emixedImps),
                class=c("sarlm"))
        rm(env)
        GC <- gc()
	if (zero.policy) {
		zero.regs <- attr(listw$neighbours, 
			"region.id")[which(card(listw$neighbours) == 0)]
		if (length(zero.regs) > 0L)
			attr(ret, "zero.regs") <- zero.regs
	}
	if (!is.null(na.act))
		ret$na.action <- na.act
	ret
}

sar_error_sse <- function(lambda, env) {
    if (get("compiled_sse", envir=env)) {
        ft <- get("first_time", envir=env)
        SSE <- .Call("R_ml_sse_env", env, lambda, PACKAGE="spdep")
        if (ft) assign("first_time", FALSE, envir=env)
    } else {
        yl <- get("y", envir=env) - lambda * get("wy", envir=env)
        yl <- get("sw", envir=env) * yl
        xl <- get("x", envir=env) - lambda * get("WX", envir=env)
        xl <- get("sw", envir=env) * xl
	xl.q <- qr.Q(qr(xl, LAPACK=get("LAPACK", envir=env)))
	xl.q.yl <- crossprod(xl.q, yl)
#        xl.q.yl <- qr.qty(qr(xl, LAPACK=get("LAPACK", envir=env)), yl)
	SSE <- crossprod(yl) - crossprod(xl.q.yl)
    }
    SSE
}


sar.error.f <- function(lambda, env) {
    SSE <- sar_error_sse(lambda, env)
    n <- get("n", envir=env)
    s2 <- SSE/n
    ldet <- do_ldet(lambda, env)
    ret <- (ldet + (1/2)*get("sum_lw", envir=env) - ((n/2)*log(2*pi)) - 
        (n/2)*log(s2) - (1/(2*(s2)))*SSE)
    if (get("verbose", envir=env)) cat("lambda:", lambda, " function:", ret, " Jacobian:", ldet, " SSE:", SSE, "\n")
    assign("f_calls", get("f_calls", envir=env)+1L, envir=env)
    ret
}

lmSLX <- function(formula, data = list(), listw, na.action, weights=NULL, zero.policy=NULL) {
        if (is.null(zero.policy))
            zero.policy <- get("zeroPolicy", envir = .spdepOptions)
        stopifnot(is.logical(zero.policy))
        if (class(formula) != "formula") formula <- as.formula(formula)
#	mt <- terms(formula, data = data)
#	mf <- lm(formula, data, na.action=na.action, weights=weights,
#            method="model.frame")
        mf <- match.call(expand.dots = FALSE)
        m <- match(c("formula", "data", "weights", "na.action"), names(mf), 0)
        mf <- mf[c(1, m)]
        mf$drop.unused.levels <- TRUE
        mf[[1]] <- as.name("model.frame")
        mf <- eval(mf, parent.frame())
        mt <- attr(mf, "terms")

	na.act <- attr(mf, "na.action")
	if (!inherits(listw, "listw")) stop("No neighbourhood list")
	if (!is.null(na.act)) {
	    subset <- !(1:length(listw$neighbours) %in% na.act)
	    listw <- subset(listw, subset, zero.policy=zero.policy)
	}
        
	y <- model.response(mf, "numeric")
	if (any(is.na(y))) stop("NAs in dependent variable")
	x <- model.matrix(mt, mf)
	if (any(is.na(x))) stop("NAs in independent variable")
        n <- nrow(x)
        if (n != length(listw$neighbours))
            stop("listw and data of different lengths")
        nclt <- colnames(x)

        weights <- as.vector(model.extract(mf, "weights"))
        if (!is.null(weights) && !is.numeric(weights)) 
            stop("'weights' must be a numeric vector")
        if (is.null(weights)) weights <- rep(as.numeric(1), n)
        if (any(is.na(weights))) stop("NAs in weights")
        if (any(weights < 0)) stop("negative weights")

        WX <- create_WX(x, listw, zero.policy=zero.policy, prefix="lag")
        x <- cbind(x, WX)
        lm.model <- lm(y ~ x - 1, weights=weights)

        sum_lm_model <- summary.lm(lm.model, correlation = FALSE)
        mixedImps <- NULL
	K <- ifelse(isTRUE(grep("\\(Intercept\\)",
            names(coefficients(lm.model))[1]) == 1L), 2, 1)
        m <- length(coefficients(lm.model))
        odd <- (m%/%2) > 0
        if (odd) {
            m2 <- (m-1)/2
        } else {
            m2 <- m/2
        }
        if (K == 1 && odd) {
            warning("model configuration issue: no total impacts")
        } else {
            cm <- matrix(0, ncol=m, nrow=m2)
            if (K == 2) {
                if (odd) {
                    rownames(cm) <- nclt[2:(m2+1)]
                } else {
                    rownames(cm) <- nclt[1:m2]
                 }
                for (i in 1:m2) cm[i, c(i+1, i+(m2+1))] <- 1
# drop bug fix 2016-09-21 Philipp Hunziker
                dirImps <- sum_lm_model$coefficients[2:(m2+1), 1:2, drop=FALSE]
                rownames(dirImps) <- rownames(cm)
                indirImps <- sum_lm_model$coefficients[(m2+2):m, 1:2, drop=FALSE]
                rownames(indirImps) <- rownames(cm)
            } else {
                rownames(cm) <- nclt[1:m2] # FIXME
                for (i in 1:m2) cm[i, c(i, i+m2)] <- 1
                dirImps <- sum_lm_model$coefficients[1:m2, 1:2, drop=FALSE]
                rownames(dirImps) <- rownames(cm)
                indirImps <- sum_lm_model$coefficients[(m2+1):m, 1:2, drop=FALSE]
                rownames(indirImps) <- rownames(cm)
            }
            totImps <- as.matrix(estimable(lm.model, cm)[, 1:2, drop=FALSE])
            mixedImps <- list(dirImps=dirImps, indirImps=indirImps,
                totImps=totImps)
        }
        attr(lm.model, "mixedImps") <- mixedImps
        class(lm.model) <- c("SLX", class(lm.model))
        lm.model
}


predict.SLX <- function(object, newdata, listw, zero.policy=NULL, ...) {
    if (is.null(zero.policy))
        zero.policy <- get("zeroPolicy", envir = .spdepOptions)
    stopifnot(is.logical(zero.policy))
    if (missing(newdata)) {
        return(fitted(object))
    }
    if (!inherits(listw, "listw")) stop("No neighbourhood list")
    if (is(newdata, "Spatial")) newdata <- as(newdata, "data.frame")
    if (!inherits(newdata, "data.frame"))
        stop("newdata must be a Spatial*DataFrame or a data.frame")
    vars <- rownames(attr(object, "mixedImps")$dirImps)
    f <- formula(paste("~", paste(vars, collapse=" + ")))
    mf <- lm(f, newdata, method="model.frame")
    mt <- attr(mf, "terms")
    x <- model.matrix(mt, mf)
    if (any(is.na(x))) stop("NAs in independent variable")
    n <- nrow(x)
    if (n != length(listw$neighbours))
        stop("listw and data of different lengths")
    WX <- create_WX(x, listw, zero.policy=zero.policy, prefix="lag")
    x <- cbind(x, WX)
    res <- as.vector(x %*% coef(object))
    names(res) <- row.names(newdata)
    res
}

create_WX <- function(x, listw, zero.policy=NULL, prefix="") {
        if (is.null(zero.policy))
            zero.policy <- get("zeroPolicy", envir = .spdepOptions)
        stopifnot(is.logical(zero.policy))
	if (!inherits(listw, "listw")) stop("No neighbourhood list")
	if (NROW(x) != length(listw$neighbours))
	    stop("Input data and neighbourhood list have different dimensions")
	n <- NROW(x)
	m <- NCOL(x)
	# check if there are enough regressors
	xcolnames <- colnames(x)
        stopifnot(!is.null(xcolnames))
	K <- ifelse(xcolnames[1] == "(Intercept)", 2, 1)
        Wvars <- NULL
        wxI <- NULL
        WX <- NULL
	if (K == 2) {
        # unnormalized weight matrices
               	if (!(listw$style == "W")) {
 			intercept <- as.double(rep(1, n))
       	       		wxI <- lag.listw(listw, intercept, 
				zero.policy = zero.policy)
                        Wvars <- paste(prefix, ".(Intercept)", sep="")
               	} 
        }   
	if (m > 1 || (m == 1 && K == 1)) {
                WX <- matrix(as.numeric(NA), nrow=n,
                    ncol=ifelse(m==1, 1, (m-(K-1))))
		for (k in K:m) {
                        j <- ifelse(k==1, 1, k-(K-1))
			WX[,j] <- lag.listw(listw, x[,xcolnames[k]], 
			    zero.policy=zero.policy)
			if (any(is.na(WX[,j]))) 
			    stop("NAs in lagged independent variable")
                        Wvars <- c(Wvars, paste(prefix, ".",
                            xcolnames[k], sep=""))
		}
	}
        if (!is.null(wxI)) WX <- cbind(wxI, WX)
        colnames(WX) <- Wvars
        WX
}

