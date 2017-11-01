# Copyright 1998-2012 by Roger Bivand and Andrew Bernat
#

lagsarlm <- function(formula, data = list(), listw, 
	na.action, type="lag", method="eigen", quiet=NULL, 
	zero.policy=NULL, interval=NULL, tol.solve=1.0e-10, 
	trs=NULL, control=list()) {
        timings <- list()
        .ptime_start <- proc.time()
        con <- list(tol.opt=.Machine$double.eps^0.5,
            fdHess=NULL, optimHess=FALSE, optimHessMethod="optimHess",
            compiled_sse=FALSE, Imult=2,
            cheb_q=5, MC_p=16L, MC_m=30L, super=NULL, spamPivot="MMD",
            in_coef=0.1, type="MC", correct=TRUE, trunc=TRUE,
            SE_method="LU", nrho=200, interpn=2000, small_asy=TRUE,
            small=1500, SElndet=NULL, LU_order=FALSE, pre_eig=NULL)
        nmsC <- names(con)
        con[(namc <- names(control))] <- control
        if (length(noNms <- namc[!namc %in% nmsC])) 
            warning("unknown names in control: ", paste(noNms, collapse = ", "))
        if (is.null(quiet)) quiet <- !get("verbose", envir = .spdepOptions)
        stopifnot(is.logical(quiet))
        if (is.null(zero.policy))
            zero.policy <- get.ZeroPolicyOption()
        stopifnot(is.logical(zero.policy))
        if (class(formula) != "formula") formula <- as.formula(formula)
	mt <- terms(formula, data = data)
	mf <- lm(formula, data, na.action=na.action, 
		method="model.frame")
	na.act <- attr(mf, "na.action")
	if (!inherits(listw, "listw")) stop("No neighbourhood list")
	can.sim <- FALSE
	if (listw$style %in% c("W", "S")) 
		can.sim <- can.be.simmed(listw)
	if (!is.null(na.act)) {
	    subset <- !(1:length(listw$neighbours) %in% na.act)
	    listw <- subset(listw, subset, zero.policy=zero.policy)
	}
        if (type == "Durbin") type <- "mixed"
	switch(type, lag = if (!quiet) cat("\nSpatial lag model\n"),
	    mixed = if (!quiet) cat("\nSpatial mixed autoregressive model\n"),
	    stop("\nUnknown model type\n"))
	y <- model.extract(mf, "response")
	x <- model.matrix(mt, mf)
	if (NROW(x) != length(listw$neighbours))
		stop("Input data and weights have different dimensions")
	n <- NROW(x)
	m <- NCOL(x)
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
	xcolnames <- colnames(x)
	K <- ifelse(xcolnames[1] == "(Intercept)", 2, 1)
	wy <- lag.listw(listw, y, zero.policy=zero.policy)
	if (any(is.na(wy))) stop("NAs in lagged dependent variable")
	if (type != "lag") {
                WX <- create_WX(x, listw, zero.policy=zero.policy,
                    prefix="lag")
		x <- cbind(x, WX)
		m <- NCOL(x)
		rm(WX)
	}
# added aliased after trying boston with TOWN dummy
	lm.base <- lm(y ~ x - 1)
	aliased <- is.na(coefficients(lm.base))
	cn <- names(aliased)
	names(aliased) <- substr(cn, 2, nchar(cn))
	if (any(aliased)) {
		nacoef <- which(aliased)
		x <- x[,-nacoef]
	}
	LL_null_lm <- logLik(lm(y ~ 1))
	m <- NCOL(x)
	similar <- FALSE
	lm.null <- lm(y ~ x - 1)
        logLik_lm.model <- logLik(lm.null)
        AIC_lm.model <- AIC(lm.null)
	lm.w <- lm.fit(x, wy)
	e.null <- lm.null$residuals
	e.w <- lm.w$residuals
	e.a <- t(e.null) %*% e.null
	e.b <- t(e.w) %*% e.null
	e.c <- t(e.w) %*% e.w
#        env <- new.env(parent=globalenv())
        env <- new.env()
        assign("y", y, envir=env)
        assign("wy", wy, envir=env)
        assign("x", x, envir=env)
        assign("n", n, envir=env)
        assign("m", m, envir=env)
        assign("K", K, envir=env)
        assign("e.a", e.a, envir=env)
        assign("e.b", e.b, envir=env)
        assign("e.c", e.c, envir=env)
        assign("family", "SAR", envir=env)
        assign("verbose", !quiet, envir=env)
        assign("compiled_sse", con$compiled_sse, envir=env)
        assign("can.sim", can.sim, envir=env)
        assign("listw", listw, envir=env)
        assign("similar", FALSE, envir=env)
        assign("f_calls", 0L, envir=env)
        assign("hf_calls", 0L, envir=env)
        timings[["set_up"]] <- proc.time() - .ptime_start
        .ptime_start <- proc.time()
	if (!quiet) cat("Jacobian calculated using ")

        interval <- jacobianSetup(method, env, con, pre_eig=con$pre_eig,
            trs=trs, interval=interval)
        assign("interval", interval, envir=env)

        nm <- paste(method, "set_up", sep="_")
        timings[[nm]] <- proc.time() - .ptime_start
        .ptime_start <- proc.time()
	opt <- optimize(sar.lag.mixed.f, interval=interval, 
		maximum=TRUE, tol=con$tol.opt, env=env)
	rho <- opt$maximum
        if (isTRUE(all.equal(rho, interval[1])) ||
            isTRUE(all.equal(rho, interval[2]))) 
            warning("rho on interval bound - results should not be used")
	names(rho) <- "rho"
	LL <- opt$objective
	optres <- opt
        nm <- paste(method, "opt", sep="_")
        timings[[nm]] <- proc.time() - .ptime_start
        .ptime_start <- proc.time()
	lm.lag <- lm((y - rho*wy) ~ x - 1)
	r <- residuals(lm.lag)
	fit <- y - r
	names(r) <- names(fit)
	coef.rho <- coefficients(lm.lag)
        tarX <- model.matrix(lm.lag)
        tary <- model.response(model.frame(lm.lag))
	names(coef.rho) <- colnames(x)
	SSE <- deviance(lm.lag)
	s2 <- SSE/n
        timings[["coefs"]] <- proc.time() - .ptime_start
        .ptime_start <- proc.time()
        assign("first_time", TRUE, envir=env)
        LMtest <- NULL
	varb <- FALSE
	ase <- FALSE
        force_assign_eigen <- FALSE
	if (method == "eigen" || do_asy) {
		rest.se <- NULL
		rho.se <- NULL
		tr <- function(A) sum(diag(A))
# beware of complex eigenvalues!
                if (do_asy && method != "eigen") eigen_setup(env)
                eig <- get("eig", envir=env)
		O <- (eig/(1-rho*eig))^2
		omega <- sum(O)
		if (is.complex(omega)) omega <- Re(omega)
		W <- listw2mat(listw)
		A <- solve(diag(n) - rho*W)
		AW <- A %*% W
		zero <- rbind(rep(0,length(coef.rho)))
		xtawxb <- s2*(t(x) %*% AW %*% x %*% coef.rho)
#		V <- s2*(s2*tr(t(AW) %*% AW) +
#			t(AW %*% x %*% coef.rho) %*%
#			(AW %*% x %*% coef.rho)) + omega*s2^2
		V <- s2*(s2*tr(crossprod(AW)) +
			crossprod(AW %*% x %*% coef.rho)) + omega*s2^2
		inf1 <- rbind(n/2, s2*tr(AW), t(zero))
		inf2 <- rbind(s2*tr(AW), V, xtawxb)
#		xtx <- s2*t(x) %*% x
		xtx <- s2*crossprod(x)
		inf3 <- rbind(zero, t(xtawxb), xtx)
		inf <- cbind(inf1, inf2, inf3)
		varb <- try(solve(inf, tol=tol.solve), silent=TRUE)
                if (class(varb) == "try-error") {
                    timings[["eigen_se"]] <- proc.time() - .ptime_start
                    .ptime_start <- proc.time()
                    con$fdHess <- TRUE
                    force_assign_eigen <- TRUE
                    warning(paste("inversion of asymptotic covariance",
                        "matrix failed for tol.solve =", tol.solve,
                        "\n", strsplit(attr(varb, "condition")$message,
                            ":")[[1]][2], "- using numerical Hessian."))
                } else {
                    varb <- (s2^2) * varb
		    rownames(varb) <- colnames(varb) <- 
			c("sigma", "rho", colnames(x))
		    rest.se <- sqrt(diag(varb))[-c(1:2)]
		    rho.se <- sqrt(varb[2,2])
		    TW <- (W %*% W) + crossprod(W)
		    T22 <- sum(diag(TW))
		    T21A <- sum(diag(TW %*% A))
		    LMtest <- ((t(r) %*% W %*% r)/s2)^2
		    LMtest <- LMtest/(T22 - ((T21A^2)*(rho.se^2)))
		    ase <- TRUE
                    timings[["eigen_se"]] <- proc.time() - .ptime_start
                    .ptime_start <- proc.time()
                }
	}

        if (con$fdHess) {
            coefs <- c(rho, coef.rho)
            if (con$compiled_sse) {
               ptr <- .Call("hess_lag_init", PACKAGE="spdep")
               assign("ptr", ptr, envir=env)
            }
            fdHess <- getVmatl(coefs, env,
               s2, trs, tol.solve=tol.solve, optim=con$optimHess,
               optimM=con$optimHessMethod)
            if (con$compiled_sse) {
                .Call("hess_lag_free", get("ptr", envir=env),
                     PACKAGE="spdep")
            }
            if (is.null(trs)) {
                rownames(fdHess) <- colnames(fdHess) <- 
                    c("rho", colnames(x))
            } else {
                rownames(fdHess) <- colnames(fdHess) <- 
                    c("sigma2", "rho", colnames(x))
            }
            if (is.null(trs)) {
 	        if (method != "eigen" || force_assign_eigen) {
                    rest.se <- sqrt(diag(fdHess)[-1])
		    rho.se <- sqrt(fdHess[1,1])
                }
            } else {
 	        if (method != "eigen" || force_assign_eigen) {
 	            rest.se <- sqrt(diag(fdHess)[-c(1,2)])
	            rho.se <- sqrt(fdHess[2,2])
                }
            }

            nm <- paste(method, "fdHess", sep="_")
            timings[[nm]] <- proc.time() - .ptime_start
            .ptime_start <- proc.time()
        } else fdHess <- FALSE
	call <- match.call()
        if (method=="SE_classic") {
            iC <- get("intern_classic", envir=env)
        } else iC <- NULL
	ret <- structure(list(type=type, rho=rho, 
		coefficients=coef.rho, rest.se=rest.se, 
		LL=LL, s2=s2, SSE=SSE, parameters=(m+2), #lm.model=lm.null,
                logLik_lm.model=logLik_lm.model, AIC_lm.model=AIC_lm.model,
		method=method, call=call, residuals=r, opt=optres,
                tarX=tarX, tary=tary, y=y, X=x,
		#lm.target=lm.lag, 
                fitted.values=fit,
		se.fit=NULL, #formula=formula,
                similar=similar,
		ase=ase, rho.se=rho.se, LMtest=LMtest, 
		resvar=varb, zero.policy=zero.policy, aliased=aliased,
                listw_style=listw$style, interval=interval, fdHess=fdHess,
                optimHess=con$optimHess, insert=!is.null(trs), trs=trs,
                LLNullLlm=LL_null_lm,
                timings=do.call("rbind", timings)[, c(1, 3)], 
                f_calls=get("f_calls", envir=env),
                hf_calls=get("hf_calls", envir=env), intern_classic=iC),
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

sar.lag.mixed.f <- function(rho, env) {
        e.a <- get("e.a", envir=env)
        e.b <- get("e.b", envir=env)
        e.c <- get("e.c", envir=env)
	SSE <- e.a - 2*rho*e.b + rho*rho*e.c
        n <- get("n", envir=env)
	s2 <- SSE/n
	ldet <- do_ldet(rho, env)
	ret <- (ldet - ((n/2)*log(2*pi)) - (n/2)*log(s2)
		- (1/(2*s2))*SSE)
	if (get("verbose", envir=env)) cat("rho:\t", rho, "\tfunction value:\t", ret, "\n")
        assign("f_calls", get("f_calls", envir=env)+1L, envir=env)

	ret
}




