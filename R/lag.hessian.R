# Copyright 2009-2013 by Roger Bivand
#

getVmatl <- function(coefs, env, s2, trs, tol.solve=1.0e-10, optim=FALSE,
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
    } else {
        fd <- fdHess(coefs, f_laglm_hess, env)
        mat <- fd$Hessian
    }
    if (!is.null(trs)) {
         mat <- insert_asy(coefs, env, s2, mat, trs)
    }
    res <- solve(-(mat), tol.solve=tol.solve)
    res
}

sar_lag_hess_sse <- function(rho, beta, env) {
    if (get("compiled_sse", envir=env)) {
        ft <- get("first_time", envir=env)
        SSE <- .Call("R_ml2_sse_env", env, rho, beta, PACKAGE="spdep")
        if (ft) assign("first_time", FALSE, envir=env)
    } else {
        res <- (get("y", envir=env) - rho * get("wy", envir=env)) - 
            get("x", envir=env) %*% beta
        SSE <- c(crossprod(res))
    }
    SSE
}

f_laglm_hess <- function(coefs, env) {
    rho <- coefs[1]
    int <- get("interval", envir=env)
    if (rho <= int[1] || rho >= int[2]) return(-Inf)
    beta <- coefs[-1]
    SSE <- sar_lag_hess_sse(rho, beta, env)
    n <- get("n", envir=env)
    s2 <- SSE/n
    det <- do_ldet(rho, env)
    ret <- (det - ((n/2) * log(2 * pi)) - (n/2) * log(s2) - 
        (1/(2 * s2)) * SSE)
    if (get("verbose", envir=env)) cat("Hessian: rho:\t", rho, "\tfunction value:\t", ret, "\n")
    assign("hf_calls", get("hf_calls", envir=env)+1L, envir=env)
    if (!is.finite(ret)) return(-Inf)
    ret
}

f_laglm_hess_nlm <- function(coefs, env) {
    ret <- f_laglm_hess(coefs, env)
    -ret
}

#f_errlm_hess <- function(coefs, env) {
#    lambda <- coefs[1]
#    beta <- coefs[-1]
#    SSE <- sar_error_hess_sse(lambda, beta, env)
#    n <- get("n", envir=env)
#    s2 <- SSE/n
#    det <- do_ldet(lambda, env)
#    ret <- (det - ((n/2) * log(2 * pi)) - (n/2) * log(s2) - 
#        (1/(2 * s2)) * SSE)
#    if (get("verbose", envir=env)) cat("lambda:", lambda, " function:", ret,
#        " Jacobian:", det, " SSE:", SSE, "\n")
#   ret
#}

trB <- function(rho, tr)  sum(sapply(0:(length(tr)-1L),
    function(i) rho^i * tr[i+1]))

insert_asy <- function(coefs, env, s2, mat, trs) {
    p <- length(coefs)-1L
    p2 <- p+2
    n <- get("n", envir=env)
    omat <- matrix(0, nrow=p2, ncol=p2)
    omat[3:p2, 3:p2] <- -crossprod(get("x", envir=env))/s2
    omat[2, 2] <- mat[1, 1]
    omat[2, 3:p2] <- omat[3:p2, 2] <- -c(crossprod(get("wy", envir=env),
        get("x", envir=env))/s2)
    omat[1, 1] <- -n/(2*(s2^2))
    omat[1, 2] <- omat[2, 1] <- -trB(coefs[1], trs)/s2
    omat
}

