# Copyright 2010-12 by Roger Bivand
sacsarlm <- function(formula, data = list(), listw, listw2=NULL, na.action, 
	type="sac", method="eigen", quiet=NULL, zero.policy=NULL, 
	tol.solve=1.0e-10, llprof=NULL, interval1=NULL, interval2=NULL,
        trs1=NULL, trs2=NULL, control=list()) {
        timings <- list()
        .ptime_start <- proc.time()
        con <- list(fdHess=NULL, #LAPACK=FALSE,
           Imult=2L, cheb_q=5L, MC_p=16L, MC_m=30L, spamPivot="MMD",
           in_coef=0.1, super=NULL, opt_method="nlminb", opt_control=list(),
           pars=NULL, npars=4L, pre_eig1=NULL, pre_eig2=NULL)
        nmsC <- names(con)
        con[(namc <- names(control))] <- control
        if (length(noNms <- namc[!namc %in% nmsC])) 
            warning("unknown names in control: ", paste(noNms, collapse = ", "))
        if (is.null(quiet)) quiet <- !get("verbose", envir = .spdepOptions)
	switch(type, sac = if (!quiet) cat("\nSpatial ARAR model\n"),
	    sacmixed = if (!quiet) cat("\nSpatial ARAR mixed model (Manski)\n"),
	    stop("\nUnknown model type\n"))
        stopifnot(is.logical(quiet))
        if (is.null(zero.policy))
            zero.policy <- get("zeroPolicy", envir = .spdepOptions)
        stopifnot(is.logical(zero.policy))
        if (class(formula) != "formula") formula <- as.formula(formula)
	mt <- terms(formula, data = data)
	mf <- lm(formula, data, na.action=na.action, method="model.frame")
	na.act <- attr(mf, "na.action")
	if (!inherits(listw, "listw")) stop("No neighbourhood list")
        if (is.null(listw2)) listw2 <- listw
        if (!is.null(con$pre_eig1) && is.null(con$pre_eig2))
            con$pre_eig2 <- con$pre_eig1
        else if (!inherits(listw2, "listw")) stop("No 2nd neighbourhood list")
        if (is.null(con$fdHess)) con$fdHess <- method != "eigen"
        if (!is.null(con$pars)) {
            stopifnot(is.numeric(con$pars))
        }
#        stopifnot(is.integer(con$npars))
#        stopifnot(is.logical(con$fdHess))
#        stopifnot(is.logical(con$LAPACK))
#        stopifnot(is.logical(con$super))
	can.sim <- FALSE
	if (listw$style %in% c("W", "S")) 
		can.sim <- can.be.simmed(listw)
	if (!is.null(na.act)) {
	    subset <- !(1:length(listw$neighbours) %in% na.act)
	    listw <- subset(listw, subset, zero.policy=zero.policy)
	}
	can.sim2 <- FALSE
	if (listw2$style %in% c("W", "S")) 
		can.sim2 <- can.be.simmed(listw2)
	if (!is.null(na.act)) {
	    subset <- !(1:length(listw2$neighbours) %in% na.act)
	    listw2 <- subset(listw2, subset, zero.policy=zero.policy)
	}
	y <- model.response(mf, "numeric")
	if (any(is.na(y))) stop("NAs in dependent variable")
	x <- model.matrix(mt, mf)
	if (any(is.na(x))) stop("NAs in independent variable")
	if (NROW(x) != length(listw$neighbours))
	    stop("Input data and neighbourhood list have different dimensions")
	wy <- lag.listw(listw, y, zero.policy=zero.policy)
	n <- NROW(x)
	m <- NCOL(x)
	if (type != "sac") {
                WX <- create_WX(x, listw, zero.policy=zero.policy,
                    prefix="lag")
		x <- cbind(x, WX)
		m <- NCOL(x)
		rm(WX)
	}
	if (NROW(x) != length(listw2$neighbours))
	    stop("Input data and neighbourhood list2 have different dimensions")
	w2y <- lag.listw(listw2, y, zero.policy=zero.policy)
	w2wy <- lag.listw(listw2, wy, zero.policy=zero.policy)
	lm.base <- lm(y ~ x - 1)
	aliased <- is.na(coefficients(lm.base))
	cn <- names(aliased)
	names(aliased) <- substr(cn, 2, nchar(cn))
	if (any(aliased)) {
		nacoef <- which(aliased)
		x <- x[,-nacoef]
	}
        LL_null_lm <- NULL
	if ("(Intercept)" %in% colnames(x)) LL_null_lm <- logLik(lm(y ~ 1))
	m <- NCOL(x)
	xcolnames <- colnames(x)
	K <- ifelse(xcolnames[1] == "(Intercept)", 2, 1)
	if (any(is.na(wy)))
	    stop("NAs in lagged dependent variable")
	if (m > 1 || (m == 1 && K == 1)) {
	    WX <- matrix(nrow=n,ncol=(m-(K-1)))
	    for (k in K:m) {
		wx <- lag.listw(listw2, x[,k], zero.policy=zero.policy)
		if (any(is.na(wx)))
		    stop("NAs in lagged independent variable")
		WX[,(k-(K-1))] <- wx
	    }
	}
	if (K == 2) {
		wx1 <- as.double(rep(1, n))
		wx <- lag.listw(listw2, wx1, zero.policy=zero.policy)
		if (m > 1) WX <- cbind(wx, WX)
		else WX <- matrix(wx, nrow=n, ncol=1)
	}
	colnames(WX) <- xcolnames
	rm(wx)

#        env <- new.env(parent=globalenv())
        env <- new.env()
        assign("y", y, envir=env)
        assign("x", x, envir=env)
        assign("wy", wy, envir=env)
        assign("w2y", w2y, envir=env)
        assign("w2wy", w2wy, envir=env)
        assign("WX", WX, envir=env)
        assign("n", n, envir=env)
        assign("p", m, envir=env)
        assign("verbose", !quiet, envir=env)
        assign("family", "SAR", envir=env)
        assign("first_time", TRUE, envir=env)
#        assign("LAPACK", con$LAPACK, envir=env)
        assign("can.sim", can.sim, envir=env)
        assign("can.sim2", can.sim2, envir=env)
        assign("listw", listw, envir=env)
        assign("listw2", listw2, envir=env)
        assign("similar", FALSE, envir=env)
        assign("similar2", FALSE, envir=env)
        timings[["set_up"]] <- proc.time() - .ptime_start
        .ptime_start <- proc.time()

	if (!quiet) cat(paste("\nSpatial autoregressive error model\n", 
		"Jacobian calculated using "))

        interval1 <- jacobianSetup(method, env, con, pre_eig=con$pre_eig1,
            trs=trs1, interval=interval1, which=1)
        assign("interval1", interval1, envir=env)
        interval2 <- jacobianSetup(method, env, con, pre_eig=con$pre_eig2,
            trs=trs2, interval=interval2, which=2)
        assign("interval2", interval2, envir=env)

        nm <- paste(method, "set_up", sep="_")
        timings[[nm]] <- proc.time() - .ptime_start
        .ptime_start <- proc.time()

        pars <- con$pars
        lower <- c(interval1[1], interval2[1])
        upper <- c(interval1[2], interval2[2])

        if (!is.null(llprof)) {
            llrho <- NULL
            lllambda <- NULL
            if (length(llprof) == 1L) {
                llrho <- seq(lower[1], upper[1], length.out=llprof)
                lllambda <- seq(lower[2], upper[2], length.out=llprof)
                llprof <- as.matrix(expand.grid(llrho, lllambda))
            }
            ll_prof <- numeric(nrow(llprof))
            for (i in 1:nrow(llprof)) 
                ll_prof[i] <- sacsar.f(llprof[i,], env=env)
            nm <- paste(method, "profile", sep="_")
            timings[[nm]] <- proc.time() - .ptime_start
            .ptime_start <- proc.time()
        }

        if (is.null(pars)) {
          if (con$npars == 4L) {
             xseq <- c(lower[1], 0, upper[1], upper[1])*0.8
             yseq <- c(upper[2], 0, upper[2], lower[2])*0.8
             mpars <- cbind(xseq, yseq)
          } else {
             xseq <- seq(lower[1], upper[1], (upper[1]-lower[1])/2)*0.8
             yseq <- seq(lower[2], upper[2], (upper[2]-lower[2])/2)*0.8
             mpars <- as.matrix(expand.grid(xseq, yseq))
          }
        } else {
            mxs <- NULL
        }
        if (con$opt_method == "nlminb") {
            if (is.null(pars)) {
                mxs <- apply(mpars, 1, function(pp) -nlminb(pp, sacsar.f, 
                    env=env, control=con$opt_control, lower=lower, 
                    upper=upper)$objective)
                pars <- mpars[which.max(mxs),]
                optres <- nlminb(pars, sacsar.f, env=env,
                    control=con$opt_control, lower=lower, upper=upper)
            } else {
                optres <- nlminb(pars, sacsar.f, env=env,
                    control=con$opt_control, lower=lower, upper=upper)
            }
        } else if (con$opt_method == "L-BFGS-B"){
            if (is.null(pars)) {
                mxs <- apply(mpars, 1, function(pp) -optim(pars, sacsar.f, 
                    env=env, method="L-BFGS-B", control=con$opt_control, 
                    lower=lower, upper=upper)$objective)
                pars <- mpars[which.max(mxs),]
	        optres <- optim(pars, sacsar.f, env=env,
                    method="L-BFGS-B", control=con$opt_control,
                    lower=lower, upper=upper)
            } else {
	        optres <- optim(pars, sacsar.f, env=env,
                    method="L-BFGS-B", control=con$opt_control,
                    lower=lower, upper=upper)
            }
        } else {
            if (is.null(pars)) {
                mxs <- apply(mpars, 1, function(pp) -optim(pars, sacsar.f, 
                    env=env, method=con$opt_method, 
                    control=con$opt_control)$objective)
                pars <- mpars[which.max(mxs),]
	        optres <- optim(pars, sacsar.f, env=env,
                    method=con$opt_method, control=con$opt_control)
            } else {
	        optres <- optim(pars, sacsar.f, env=env,
                    method=con$opt_method, control=con$opt_control)
            }
        }
	LL <- -optres$objective
        if (optres$convergence != 0)
            warning(paste("convergence failure:", optres$message))
	rho <- optres$par[1]
	names(rho) <- "rho"
	lambda <- optres$par[2]
	names(lambda) <- "lambda"
        nm <- paste(method, "opt", sep="_")
        timings[[nm]] <- proc.time() - .ptime_start
        .ptime_start <- proc.time()

	lm.target <- lm(I(y - rho*wy - lambda*w2y + rho*lambda*w2wy) ~ 
            I(x - lambda*WX) - 1)
	r <- as.vector(residuals(lm.target))
	fit <- as.vector(y - r)
	p <- lm.target$rank
	SSE <- deviance(lm.target)
	s2 <- SSE/n
	coef.sac <- coefficients(lm.target)
        tarX <- model.matrix(lm.target)
        tary <- model.response(model.frame(lm.target))
	names(coef.sac) <- xcolnames
	lm.model <- lm(formula, data)
        logLik_lm.model <- logLik(lm.model)
        AIC_lm.model <- AIC(lm.model)
        ase <- FALSE
	asyvar1 <- FALSE
        force_assign_eigen <- FALSE
        if (method == "eigen") {
# taken from spatial/sac_models/sac.m
	    tr <- function(A) sum(diag(A))
            W1 <- listw2mat(listw)
            W2 <- listw2mat(listw2)
            A <- diag(n) - rho * W1
            AI <- solve(A)
            WA <- W1 %*% AI
            B <- diag(n) - lambda * W2
            BI <- solve(B)
            WB <- W2 %*% BI
            omeg <- s2*diag(n)
            omegi <- (1/s2)*diag(n)
            bhat <- coef.sac
            p3 <- p+3
            asyvar <- matrix(0.0, ncol=p3, nrow=p3)
            Bx <- B %*% x
            asyvar[4:p3, 4:p3] <- (1/s2)*crossprod(Bx)
            asyvar[1, 1] <- n/(2*s2*s2)
            term1 <- tr(WA %*% WA)
            BWABI <- B %*% WA %*% BI
            term2 = tr(omeg %*% t(BWABI) %*% omegi %*% BWABI)
            BWAxbhat <- B %*% WA %*% x %*% bhat
            term3 = t(BWAxbhat) %*% omegi %*% BWAxbhat
            asyvar[2, 2] <- term1 + term2 + term3
            term1 <- tr(WB %*% WB)
            term2 <- tr(omeg %*% t(WB) %*% omegi %*% WB)
            asyvar[3, 3] <- term1 + term2
            asyvar[2, 4:p3] <- (1/s2)*(t(Bx) %*% BWAxbhat)
            asyvar[4:p3, 2] <- asyvar[2, 4:p3]
            asyvar[2, 1] <- (1/s2)*tr(W1 %*% AI)
            asyvar[1, 2] <- asyvar[2, 1]
            asyvar[3, 1] <- (1/s2)*tr(W2 %*% BI)
            asyvar[1, 3] <- asyvar[3, 1]
            term1 = tr(t(WB) %*% omegi %*% BWABI %*% omeg)
            term2 = tr(W2 %*% WA %*% BI)
            asyvar[3, 2] <- term1 + term2
            asyvar[2, 3] <- asyvar[3, 2]
            asyvar1 <- try(solve(asyvar, tol.solve=tol.solve), silent=TRUE)
            if (class(asyvar1) == "try-error") {
                timings[["eigen_se"]] <- proc.time() - .ptime_start
                .ptime_start <- proc.time()
                con$fdHess <- TRUE
                force_assign_eigen <- TRUE
                warning(paste("inversion of asymptotic covariance",
                    "matrix failed for tol.solve =", tol.solve,
                    "\n", strsplit(attr(asyvar1, "condition")$message,
                        ":")[[1]][2], "- using numerical Hessian."))
                ase=FALSE
            } else {
                rownames(asyvar1) <- colnames(asyvar1) <- 
			c("sigma", "rho", "lambda", xcolnames)
                ase=TRUE
                rho.se <- sqrt(asyvar1[2, 2])
                lambda.se <- sqrt(asyvar1[3, 3])
                rest.se <- sqrt(diag(asyvar1)[-c(1:3)])
                nm <- "asymptotic vcov"
                timings[[nm]] <- proc.time() - .ptime_start
                .ptime_start <- proc.time()
            }
        }
        if (con$fdHess) {
            coefs <- c(rho, lambda, coef.sac)
            fdHess <- getVmatsac(coefs, env, tol.solve=tol.solve)
            rownames(fdHess) <- colnames(fdHess) <- c("rho", "lambda",
                xcolnames)
            if (!ase) {
                rho.se <- sqrt(fdHess[1, 1])
                lambda.se <- sqrt(fdHess[2, 2])
                rest.se <- sqrt(diag(fdHess)[-c(1,2)])
            }
            nm <- paste(method, "fdHess", sep="_")
            timings[[nm]] <- proc.time() - .ptime_start
        }
	call <- match.call()
	names(r) <- names(y)
	names(fit) <- names(y)
	ret <- structure(list(type=type, rho=rho, lambda=lambda,
	    coefficients=coef.sac, rest.se=rest.se, ase=ase,
	    LL=LL, s2=s2, SSE=SSE, parameters=(p+3), 
            logLik_lm.model=logLik_lm.model, AIC_lm.model=AIC_lm.model,
            #lm.model=lm.model, 
	    method=method, call=call, residuals=r, #lm.target=lm.target,
            tarX=tarX, tary=tary, y=y, X=x, W2X=WX, trs1=trs1, trs2=trs2,
	    opt=optres, pars=pars, mxs=mxs, fitted.values=fit, #formula=formula,
	    similar=get("similar", envir=env), rho.se=rho.se,
	    lambda.se=lambda.se, zero.policy=zero.policy, 
	    aliased=aliased, LLNullLlm=LL_null_lm,
            fdHess=fdHess, resvar=asyvar1, listw_style=listw$style,
            optimHess=FALSE, insert=FALSE, interval1=interval1,
            interval2=interval2, timings=do.call("rbind", timings)[, c(1, 3)]),
            class=c("sarlm"))
        rm(env)
        GC <- gc()
        if (is.null(llprof)) ret$llprof <- llprof
        else {
            ret$llprof <- list(grd=llprof, ll=ll_prof, xseq=llrho,
                yseq=lllambda)
        }
	if (!is.null(na.act))
		ret$na.action <- na.act
	ret

}


sacsar_sse <- function(coefs, env) {
    rho <- coefs[1]
    lambda <- coefs[2]
    yl <- get("y", envir=env) - rho * get("wy", envir=env) - 
        lambda * get("w2y", envir=env) + rho * lambda * get("w2wy", envir=env)
    xl <- get("x", envir=env) - lambda * get("WX", envir=env)
    xl.q <- qr.Q(qr(xl, #LAPACK=get("LAPACK", envir=env)
))
    xl.q.yl <- crossprod(xl.q, yl)
    SSE <- crossprod(yl) - crossprod(xl.q.yl)
    SSE
}


sacsar.f <- function(coefs, env) {
    SSE <- sacsar_sse(coefs, env)
    n <- get("n", envir=env)
    s2 <- SSE/n
    ldet1 <- do_ldet(coefs[1], env, which=1)
    ldet2 <- do_ldet(coefs[2], env, which=2)
    ret <- (ldet1 + ldet2 - ((n/2)*log(2*pi)) - (n/2)*log(s2) - 
        (1/(2*(s2)))*SSE)
    if (get("verbose", envir=env)) cat("rho:", coefs[1], " lambda:", coefs[2],
        " function:", ret, " Jacobian1:", ldet1, " Jacobian2:", ldet2,
        " SSE:", SSE, "\n")
    -ret
}


getVmatsac <- function(coefs, env, tol.solve=1.0e-10) {
    fd <- fdHess(coefs, f_sac_hess, env)
    mat <- fd$Hessian
    res <- solve(-(mat), tol.solve=tol.solve)
    res
}

sar_sac_hess_sse <- function(rho, lambda, beta, env) {
    yl <- get("y", envir=env) - rho * get("wy", envir=env) - 
        lambda * get("w2y", envir=env) + rho * lambda * get("w2wy", envir=env)
    xl <- get("x", envir=env) - lambda * get("WX", envir=env)
    res <- yl - (xl %*% beta)
    SSE <- c(crossprod(res))
    SSE
}

f_sac_hess <- function(coefs, env) {
    rho <- coefs[1]
    int <- get("interval1", envir=env)
    if (rho <= int[1] || rho >= int[2]) return(-Inf)
    lambda <- coefs[2]
    int <- get("interval2", envir=env)
    if (lambda <= int[1] || lambda >= int[2]) return(-Inf)
    beta <- coefs[-(1:2)]
    SSE <- sar_sac_hess_sse(rho, lambda, beta, env)
    n <- get("n", envir=env)
    s2 <- SSE/n
    ldet1 <- do_ldet(rho, env, which=1)
    ldet2 <- do_ldet(lambda, env, which=2)
    ret <- (ldet1 + ldet2 - ((n/2) * log(2 * pi)) - (n/2) * log(s2) - 
        (1/(2 * s2)) * SSE)
    if (get("verbose", envir=env)) cat("rho:", rho, "lambda:", lambda,
        " function:", ret, " Jacobian1:", ldet1, " Jacobian2:",
        ldet2, " SSE:", SSE, "\n")
    if (!is.finite(ret)) return(-Inf)
   ret
}


