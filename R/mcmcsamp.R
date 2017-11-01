MCMCsamp <- function(object, mcmc = 1L, verbose = NULL, ...) UseMethod("MCMCsamp")
# from lme4/R/AllGeneric.R
MCMCsamp.spautolm <- function(object, mcmc = 1L, verbose = NULL, ...,
    burnin=0L, scale=1, listw, control=list()) {
    con <- list(Imult=2, cheb_q=5, MC_p=16, MC_m=30, super=NULL,
        spamPivot="MMD", in_coef=0.1, type="MC",
        correct=TRUE, trunc=TRUE, SE_method="LU", nrho=200,
        interpn=2000, small_asy=TRUE, small=1500, SElndet=NULL,
        LU_order=FALSE)
    timings <- list()
    .ptime_start <- proc.time()
    nmsC <- names(con)
    con[(namc <- names(control))] <- control
    if (length(noNms <- namc[!namc %in% nmsC])) 
        warning("unknown names in control: ", paste(noNms, collapse = ", "))
    if (is.null(verbose)) verbose <- get("verbose", envir = .spdepOptions)
    stopifnot(is.logical(verbose))
    if (!inherits(listw, "listw")) 
        stop("No neighbourhood list")
    method <- object$method
    family <- object$family
    if (family == "SMA" && method != "eigen") stop("SMA only for eigen method")
    X <- object$X
    N <- nrow(X)
    if (N != length(listw$neighbours))
	 stop("Input data and neighbourhood list have different dimensions")
    stopifnot(ncol(X) == length(object$fit$coefficients))
    weights <- object$weights
    stopifnot(length(weights) == N)
    can.sim <- FALSE
    if (listw$style %in% c("W", "S")) 
	can.sim <- can.be.simmed(listw)
    sum_lw <- sum(log(weights))
    env <- new.env()
    assign("Y", object$Y, envir=env)
    assign("X", X, envir=env)
    assign("n", N, envir=env)
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
    I <- as_dsCMatrix_I(N)
    assign("I", I, envir=env)
    Sweights <- as(as(Diagonal(x=weights), "symmetricMatrix"), 
        "CsparseMatrix")
    assign("Sweights", Sweights, envir=env)
    timings[["data_setup"]] <- proc.time() - .ptime_start
    .ptime_start <- proc.time()

    if (verbose) cat(paste("\nJacobian calculated using "))

    interval <- jacobianSetup(method, env, con, trs=object$trs,
        interval=object$interval)
    assign("interval", interval, envir=env)
    timings[["logdet_setup"]] <- proc.time() - .ptime_start
    .ptime_start <- proc.time()

    start <- c(object$lambda, object$fit$coefficients)
    V <- object$fdHess
    stopifnot(nrow(V) == ncol(V))
    stopifnot(nrow(V) == length(start))
    res0 <- rwmetrop(logpost=f_spautolm_hess, start=start,
        proposal=list(var=V, scale=scale), m=(mcmc+burnin), env=env)
#        MCMCmetrop1R(fun=f_spautolm_hess, theta.init=start, burnin=burnin,
#        mcmc=mcmc, thin=thin, verbose=ifelse(verbose, verbose_step, 0),
#        seed=seed, logfun=TRUE, V=V, env=env)
    res <- as.mcmc(res0$par[(burnin+1):(mcmc+burnin),])
    attr(res, "accept") <- res0$accept
    colnames(res) <- c(names(object$lambda), names(object$fit$coefficients))
    attr(res, "type") <- family
    timings[["rwmetrop"]] <- proc.time() - .ptime_start
    attr(res, "timings") <- do.call("rbind", timings)[, c(1, 3)]
    res
}

MCMCsamp.sarlm <- function(object, mcmc = 1L, verbose = NULL, ...,
    burnin=0L, scale=1, listw, listw2=NULL, control=list()) {
    timings <- list()
    .ptime_start <- proc.time()
    con <- list(compiled_sse=FALSE, Imult=2, cheb_q=5, MC_p=16, MC_m=30,
        super=NULL, spamPivot="MMD", in_coef=0.1, type="MC",
        correct=TRUE, trunc=TRUE, SE_method="LU", nrho=200,
        interpn=2000, small_asy=TRUE, small=1500, SElndet=NULL,
        LU_order=FALSE)
    nmsC <- names(con)
    con[(namc <- names(control))] <- control
    if (length(noNms <- namc[!namc %in% nmsC])) 
        warning("unknown names in control: ", paste(noNms, collapse = ", "))
    if (is.null(verbose)) verbose <- get("verbose", envir = .spdepOptions)
    stopifnot(is.logical(verbose))
    if (!inherits(listw, "listw")) 
        stop("No neighbourhood list")
    type <- object$type
    if (length(grep("sac", type) > 0)) Type <- "SAC"
    else if (length(grep("error", type) > 0)) Type <- "ERROR"
    else Type <- "LAG"
    method <- object$method
    X <- object$X
    N <- nrow(X)
    if (N != length(listw$neighbours))
	 stop("Input data and neighbourhood list have different dimensions")
    stopifnot(ncol(X) == length(object$coefficients))
    can.sim <- FALSE
    if (listw$style %in% c("W", "S")) 
	can.sim <- can.be.simmed(listw)
    weights <- object$weights
    if (is.null(weights)) weights <- rep(1, N)
    stopifnot(length(weights) == N)
    sum_lw <- sum(log(weights))
    sw <- sqrt(weights)

    env <- new.env()
    assign("y", object$y, envir=env)
    assign("x", X, envir=env)
    assign("n", N, envir=env)
    assign("listw", listw, envir=env)
    assign("sum_lw", sum_lw, envir=env)
    assign("sw", sw, envir=env)
    assign("can.sim", can.sim, envir=env)
    assign("method", method, envir=env)
    assign("verbose", verbose, envir=env)
    assign("family", "SAR", envir=env)
    W <- as(listw, "CsparseMatrix")
    wy <- c(as.matrix(W %*% c(object$y)))
    assign("wy", wy, envir=env)

    if (is.numeric(object$resvar)) V <- object$resvar[-1,-1]
    else V <- object$fdHess
    stopifnot(nrow(V) == ncol(V))
    timings[["data_setup"]] <- proc.time() - .ptime_start
    .ptime_start <- proc.time()

    if (verbose) cat(paste("\nJacobian calculated using "))

    if (Type == "SAC") {

        interval1 <- jacobianSetup(method, env, con, trs=object$trs1,
            interval=object$interval1, which=1)
        assign("interval1", interval1, envir=env)

        if (is.null(listw2)) listw2 <- listw
        else if (!inherits(listw2, "listw")) stop("No 2nd neighbourhood list")
        can.sim2 <- FALSE
	if (listw2$style %in% c("W", "S")) 
		can.sim2 <- can.be.simmed(listw2)
        assign("listw2", listw2, envir=env)
        assign("can.sim2", can.sim2, envir=env)
        interval2 <- jacobianSetup(method, env, con, trs=object$trs2,
            interval=object$interval2, which=2)
        assign("interval2", interval2, envir=env)

    } else {

        interval <- jacobianSetup(method, env, con, trs=object$trs,
            interval=object$interval)
        assign("interval", interval, envir=env)
        assign("hf_calls", 0L, envir=env)
        assign("compiled_sse", con$compiled_sse, envir=env)
        assign("first_time", TRUE, envir=env)

    }
    timings[["logdet_setup"]] <- proc.time() - .ptime_start
    .ptime_start <- proc.time()

    if (Type == "LAG") {
        if (con$compiled_sse) {
           ptr <- .Call("hess_lag_init", PACKAGE="spdep")
           assign("ptr", ptr, envir=env)
        }
        start <- c(object$rho, object$coefficients)
        stopifnot(nrow(V) == length(start))
        res0 <- rwmetrop(logpost=f_laglm_hess, start=start,
            proposal=list(var=V, scale=scale), m=(mcmc+burnin), env=env)
#        res <- MCMCmetrop1R(fun=f_laglm_hess, theta.init=start,
#            burnin=burnin, mcmc=mcmc, thin=thin, verbose=ifelse(verbose,
#            verbose_step, 0), seed=seed, logfun=TRUE, V=V, env=env)

        if (con$compiled_sse) {
            .Call("hess_lag_free", get("ptr", envir=env),
                 PACKAGE="spdep")
        }
        res <- as.mcmc(res0$par[(burnin+1):(mcmc+burnin),])
        attr(res, "accept") <- res0$accept
        colnames(res) <- c(names(object$rho), names(object$coefficients))
        attr(res, "type") <- object$type

    } else if (Type == "ERROR") {

        WX <- as.matrix(W %*% X)
        assign("WX", WX, envir=env)

        if (con$compiled_sse) {
            ptr <- .Call("hess_error_init", PACKAGE="spdep")
            assign("ptr", ptr, envir=env)
        }
        start <- c(object$lambda, object$coefficients)
        stopifnot(nrow(V) == length(start))
        res0 <- rwmetrop(logpost=f_errlm_hess, start=start,
            proposal=list(var=V, scale=scale), m=(mcmc+burnin), env=env)
#        res <- MCMCmetrop1R(fun=f_errlm_hess, theta.init=start,
#            burnin=burnin, mcmc=mcmc, thin=thin, verbose=ifelse(verbose,
#            verbose_step, 0), seed=seed, logfun=TRUE, V=V, env=env)

        if (con$compiled_sse) {
            .Call("hess_error_free", get("ptr", envir=env),
                PACKAGE="spdep")
        }
        res <- as.mcmc(res0$par[(burnin+1):(mcmc+burnin),])
        attr(res, "accept") <- res0$accept
        colnames(res) <- c(names(object$lambda), names(object$coefficients))
        attr(res, "type") <- object$type

    } else if (Type == "SAC") {
        assign("WX", object$W2X, envir=env)
        W2 <- as(listw2, "CsparseMatrix")
        w2y <- c(as.matrix(W2 %*% object$y))
        assign("w2y", w2y, envir=env)
        w2wy <- c(as.matrix(W2 %*% wy))
        assign("w2wy", w2wy, envir=env)
        start <- c(object$rho, object$lambda, object$coefficients)
        stopifnot(nrow(V) == length(start))   
        res0 <- rwmetrop(logpost=f_sac_hess, start=start,
            proposal=list(var=V, scale=scale), m=(mcmc+burnin), env=env)
#        res <- MCMCmetrop1R(fun=f_sac_hess, theta.init=start,
#            burnin=burnin, mcmc=mcmc, thin=thin, verbose=ifelse(verbose,
#            verbose_step, 0), seed=seed, logfun=TRUE, V=V, env=env)
        res <- as.mcmc(res0$par[(burnin+1):(mcmc+burnin),])
        attr(res, "accept") <- res0$accept
        colnames(res) <- c(names(object$rho), names(object$lambda),
            names(object$coefficients))
        attr(res, "type") <- object$type
    }
    timings[["rwmetrop"]] <- proc.time() - .ptime_start
    attr(res, "timings") <- do.call("rbind", timings)[, c(1, 3)]
    if (attr(res, "accept") < 0.05)
        warning("MCMCsamp: very low acceptance rate")
    res
}

