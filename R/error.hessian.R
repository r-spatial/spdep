# Copyright 2009-2013 by Roger Bivand

getVmate <- function(coefs, env, s2, trs, tol.solve=1.0e-10, optim=FALSE,
    optimM="optimHess") {
    if (optim) {
      if (optimM == "nlm") {
           options(warn=-1)
           opt <- nlm(f=f_laglm_hess_nlm, p=coefs, env=env, hessian=TRUE)
           options(warn=0)
           mat <- opt$hessian
#        opt <- optimHess(par=coefs, fn=f_laglm_hess, env=env)
#        mat <- opt
       } else if (optimM == "optimHess") {
           mat <- optimHess(par=coefs, fn=f_laglm_hess, env=env)
       } else {
           opt <- optim(par=coefs, fn=f_laglm_hess, env=env, method=optimM,
           hessian=TRUE)
           mat <- opt$hessian
      }
#        opt <- optimHess(par=coefs, fn=f_errlm_hess, env=env)
#        mat <- opt
    } else {
        fd <- fdHess(coefs, f_errlm_hess, env)
        mat <- fd$Hessian
    }
    if (!is.null(trs)) {
         mat <- insert_asye(coefs, env, s2, mat, trs)
    }
    res <- solve(-(mat), tol.solve=tol.solve)
    res
}

sar_error_hess_sse <- function(lambda, beta, env) {
    if (get("compiled_sse", envir=env)) {
        ft <- get("first_time", envir=env)
        SSE <- .Call("R_ml1_sse_env", env, lambda, beta, PACKAGE="spdep")
        if (ft) assign("first_time", FALSE, envir=env)
    } else {
        yl <- get("y", envir=env) - lambda * get("wy", envir=env)
        xl <- get("x", envir=env) - lambda * get("WX", envir=env)
        res <- get("sw", envir=env) * (yl - (xl %*% beta))
        SSE <- c(crossprod(res))
    }
    SSE
}

f_errlm_hess <- function(coefs, env) {
    lambda <- coefs[1]
    int <- get("interval", envir=env)
    if (lambda <= int[1] || lambda >= int[2]) return(-Inf)
    beta <- coefs[-1]
    SSE <- sar_error_hess_sse(lambda, beta, env)
    n <- get("n", envir=env)
    s2 <- SSE/n
    det <- do_ldet(lambda, env)
    ret <- (det + (1/2)*get("sum_lw", envir=env) - ((n/2) * log(2 * pi)) - 
        (n/2) * log(s2) - (1/(2 * s2)) * SSE)
    if (get("verbose", envir=env)) cat("lambda:", lambda, " function:", ret,
        " Jacobian:", det, " SSE:", SSE, "\n")
    assign("hf_calls", get("hf_calls", envir=env)+1L, envir=env)
    if (!is.finite(ret)) return(-Inf)
    ret
}

insert_asye <- function(coefs, env, s2, mat, trs) {
    lambda <- coefs[1]
    p <- length(coefs)-1L
    p2 <- p+2
    omat <- matrix(0, nrow=p2, ncol=p2)
    LX <- get("sw", envir=env) * (get("x", envir=env) - lambda *
        get("WX", envir=env))
#    omat[3:p2, 3:p2] <- -crossprod(LX)*s2
#    omat[3:p2, 3:p2] <- -crossprod(LX)
    omat[3:p2, 3:p2] <- -crossprod(LX)/s2
    omat[2, 2] <- mat[1, 1]
    n <- get("n", envir=env)
    omat[1, 1] <- -n/(2*(s2^2))
#    omat[1, 1] <- -n/(2*(s2))
    omat[1, 2] <- omat[2, 1] <- -trB(lambda, trs)/s2
#    omat[1, 2] <- omat[2, 1] <- -trB(lambda, trs)
    omat
}

