# Copyright 2009-2017 by Roger Bivand

trW <- function(W=NULL, m=30, p=16, type="mult", listw=NULL, momentsSymmetry=TRUE) {
# returns traces
    timings <- list()
    .ptime_start <- proc.time()
    if (type == "mult") {
        stopifnot(!is.null(W))
        stopifnot(inherits(W, "sparseMatrix"))
        n <- dim(W)[1]
        iW <- W
        tr <- numeric(m)
        for (i in 1:m) {
            tr[i] <- sum(diag(iW))
            iW <- W %*% iW
        }
    } else if (type == "MC") {
        stopifnot(!is.null(W))
        stopifnot(inherits(W, "sparseMatrix"))
        n <- dim(W)[1]
        tr <- numeric(m)
# return sd of traces 111126
        sdtr <- numeric(m)
        x <- matrix(rnorm(n*p), nrow=n, ncol=p)
        xx <- x
        for (i in 1:m) {
            xx <- W %*% xx
# return sd of traces 111126
            v <- apply(x * as.matrix(xx), 2, sum)
            tr[i] <- mean(v)
            sdtr[i] <- sd(v)/sqrt(p)
#            tr[i] <- sum(apply(x * as.matrix(xx), 2,  function(y) sum(y)/p))
# mean replaced by sum(y)/p 091012, 0.4-47
        }
        tr[1] <- 0.0
        tr[2] <- sum(t(W) * W)
# return sd of traces 111126
        sdtr[1:2] <- NA
        attr(tr, "sd") <- sdtr
    } else if (type == "moments") {
        if (!is.null(W) && is.null(listw)) {
            if (momentsSymmetry && !is(W, "symmetricMatrix"))
                stop("moments require symmetric W")
            listw <- mat2listw(W)
        }
        tr <- mom_calc(listw, m)
        n <- length(listw$neighbours)
    } else stop("unknown type")
    timings[["make_traces"]] <- proc.time() - .ptime_start
    attr(tr, "timings") <- do.call("rbind", timings)[, c(1, 3)]
    attr(tr, "type") <- type
    attr(tr, "n") <- n
    tr
}

mom_calc_int <- function(is, m, W, eta0) {
    Omega <- rep(0.0, m)
    for (i in is) {
        eta <- eta0
        eta[i] <- 1
        for (j in seq(2, m, 2)) {
            zeta <- W %*% eta
            Omega[j-1] <- Omega[j-1] + crossprod(zeta, eta)[1,1]
            Omega[j] <- Omega[j] + crossprod(zeta, zeta)[1,1]
            eta <- zeta
        }
    }
    Omega
}

mom_calc_int2 <- function(is, m, nb, weights, Card) {
    Omega <- .Call("mom_calc_int2", is, as.integer(m), nb, weights, Card, PACKAGE="spdep")
    Omega
}

mom_calc <- function(lw, m) {
    stopifnot((m %% 2) == 0)
    nb <- lw$neighbours
    n <- length(nb)
    weights <- lw$weights
    Card <- card(nb)

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

    if (parallel == "snow") {
      if (requireNamespace("parallel", quietly = TRUE)) {
#        require(parallel)
        lis <- parallel::splitIndices(n, length(cl))
        lOmega <- parallel::parLapply(cl, lis, mom_calc_int2, m, nb, weights, Card)
        Omega <- apply(do.call("cbind", lOmega), 1, sum)
      } else {
        stop("parallel not available")
      }
    } else if (parallel == "multicore") {
      if (requireNamespace("parallel", quietly = TRUE)) {
#        require(parallel)
        lis <- parallel::splitIndices(n, ncpus)
        lOmega <- parallel::mclapply(lis, mom_calc_int2, m, nb, weights, Card,
            mc.set.seed=FALSE, mc.cores=ncpus)
        Omega <- apply(do.call("cbind", lOmega), 1, sum)
      } else {
        stop("parallel not available")
      }
    } else {
        Omega <- mom_calc_int2(is=1:n, m=m, nb=nb, weights=weights, Card=Card)
    }
    Omega
}

impacts <- function(obj, ...)
    UseMethod("impacts", obj)

impacts.SLX <- function(obj, ...) {
    stopifnot(!is.null(attr(obj, "mixedImps")))
    n <- nrow(obj$model)
    k <- obj$qr$rank
    impactsWX(attr(obj, "mixedImps"), n, k, type="SLX")
}

impactSDEM <- function(obj) {
    n <- nrow(obj$tarX)
    k <- ncol(obj$tarX)
    impactsWX(obj$emixedImps, n, k, type="SDEM")
}

impactsWX <- function(obj, n, k, type="SLX") {
    imps <- lapply(obj, function(x) x[, 1])
    names(imps) <- c("direct", "indirect", "total")
    attr(imps, "bnames") <- rownames(obj[[1]])
    ses <- lapply(obj, function(x) x[, 2])
    names(ses) <- c("direct", "indirect", "total")
    attr(ses, "bnames") <- rownames(obj[[1]])
    res <- list(impacts=imps, se=ses)
    attr(res, "n") <- n
    attr(res, "k") <- k
    attr(res, "type") <- type
    attr(res, "method") <- "estimable"
    attr(res, "bnames") <- rownames(obj[[1]])
    class(res) <- "WXImpact"
    res
}


print.WXImpact <- function(x, ...) {
    mat <- lagImpactMat(x$impacts)
    cat("Impact measures (", attr(x, "type"), ", ",
        attr(x, "method"), "):\n", sep="")
    print(mat, ...)
    cat("\n")
    invisible(x)
}

print.summary.WXImpact <- function(x, ...) {
    mat <- x$mat
    cat("Impact measures (", attr(x, "type"), ", ",
        attr(x, "method"), "):\n", sep="")
    print(mat, ...)
    cat("========================================================\n")
    mat <- x$semat
    cat("Standard errors:\n", sep="")
    print(mat, ...)
    cat("========================================================\n")
    cat("Z-values:\n")
    mat <- x$zmat
    rownames(mat) <- attr(x, "bnames")
    print(mat, ...)
    cat("\np-values:\n")
    xx <- apply(x$pzmat, 2, format.pval)
# 100928 Eelke Folmer
    if (length(attr(x, "bnames")) == 1L) {
        xx <- matrix(xx, ncol=3)
        colnames(xx) <- c("Direct", "Indirect", "Total")
    }
    rownames(xx) <- attr(x, "bnames")
    print(xx, ..., quote=FALSE)
    cat("\n")
    invisible(x)
}

summary.WXImpact <- function(object, ..., adjust_k=FALSE) {
    stopifnot(is.logical(adjust_k))
    stopifnot(length(adjust_k) == 1L)
    object$mat <- lagImpactMat(object$impacts)
    object$semat <- lagImpactMat(object$se)
    if (adjust_k) {
        object$semat <- object$semat * (attr(object, "n")/
            (attr(object, "n") - attr(object, "k")))
        attr(object, "method") <- paste(attr(object, "method"),
            ", n-k", sep="")
    }
    object$zmat <- object$mat/object$semat
    object$pzmat <- 2*(1-pnorm(abs(object$zmat)))
    class(object) <- c("summary.WXImpact", class(object))
    object
}



impacts.stsls <- function(obj, ..., tr=NULL, R=NULL, listw=NULL,
  tol=1e-6, empirical=FALSE, Q=NULL) {
    if (is.null(listw) && !is.null(obj$listw_style) && 
            obj$listw_style != "W")
            stop("Only row-standardised weights supported")
    rho <- obj$coefficients[1]
    beta <- obj$coefficients[-1]
    icept <- grep("(Intercept)", names(beta))
    iicept <- length(icept) > 0
    if (iicept) {
        P <- matrix(beta[-icept], ncol=1)
        bnames <- names(beta[-icept])
    } else {
        P <- matrix(beta, ncol=1)
        bnames <- names(beta)
    }
    p <- length(beta)
    n <- length(obj$residuals)
    mu <- c(rho, beta)
    Sigma <- obj$var
    irho <- 1
    drop2beta <- 1
    res <- intImpacts(rho=rho, beta=beta, P=P, n=n, mu=mu, Sigma=Sigma,
        irho=irho, drop2beta=drop2beta, bnames=bnames, interval=NULL,
        type="lag", tr=tr, R=R, listw=listw, tol=tol, empirical=empirical,
        Q=Q, icept=icept, iicept=iicept, p=p)
    attr(res, "iClass") <- class(obj)
    if (!is.null(obj$robust)) {
        attr(res, "robust") <- obj$robust
        attr(res, "HC") <- obj$HC
    }
    res
}

impacts.gmsar <- function(obj, ..., n=NULL, tr=NULL, R=NULL, listw=NULL,
  tol=1e-6, empirical=FALSE, Q=NULL) {
    stopifnot(obj$type == "SARAR") 
    if (is.null(listw) && !is.null(obj$listw_style) && 
            obj$listw_style != "W")
            stop("Only row-standardised weights supported")
    rho <- obj$coefficients[1]
    beta <- obj$coefficients[-1]
    icept <- grep("(Intercept)", names(beta))
    iicept <- length(icept) > 0
    if (iicept) {
        P <- matrix(beta[-icept], ncol=1)
        bnames <- names(beta[-icept])
    } else {
        P <- matrix(beta, ncol=1)
        bnames <- names(beta)
    }
    p <- length(beta)
# allow n passthrough 140305 Angela Parenti
    if (is.null(n)) n <- length(obj$residuals)
    stopifnot(is.integer(n))
    stopifnot(length(n) == 1)
    stopifnot(is.finite(n))
    mu <- c(rho, beta)
    Sigma <- obj$secstep_var
    irho <- 1
    drop2beta <- 1
    res <- intImpacts(rho=rho, beta=beta, P=P, n=n, mu=mu, Sigma=Sigma,
        irho=irho, drop2beta=drop2beta, bnames=bnames, interval=NULL,
        type="lag", tr=tr, R=R, listw=listw, tol=tol, empirical=empirical,
        Q=Q, icept=icept, iicept=iicept, p=p)
    attr(res, "iClass") <- class(obj)
    res
}


lagImpacts <- function(T, g, P) {
    PT <- P %*% T
    direct <- apply(apply(PT, 1, function(x) x*g), 2, sum)
    total <- c(apply(P, 1, sum) * sum(g))
    indirect <- total - direct
    names(direct) <- names(total)
    list(direct=direct, indirect=indirect, total=total)
}

lagDistrImpacts <- function(T, g, P, q=10) {
    PT <- P %*% T
    direct <- apply(PT, 1, function(x) x * g)[1:q, ]
    if (nrow(P) == 1) {
        total <- sapply(g, function(x) apply(P, 1, sum)*x)[1:q]
    } else {
        total <- t(sapply(g, function(x) apply(P, 1, sum)*x))[1:q, ]
    }
    indirect <- total - direct
    list(direct=direct, indirect=indirect, total=total)
}

processSample <- function(x, irho, drop2beta, type, iicept, icept, T, Q, q) {
    g <- x[irho]^(0:q)
    beta <- x[-drop2beta]
    if (type == "lag" || type == "sac") {
      if (iicept) {
        P <- matrix(beta[-icept], ncol=1)
      } else {
        P <- matrix(beta, ncol=1)
      }
    } else if (type == "mixed" || type == "sacmixed") {
        if (iicept) {
          b1 <- beta[-icept]
        } else {
          b1 <- beta
        }
        p <- length(b1)
        if (p %% 2 != 0) stop("non-matched coefficient pairs")
        P <- cbind(b1[1:(p/2)], b1[((p/2)+1):p])
    }
    res <- lagImpacts(T, g, P)
    if (!is.null(Q)) {
        Qres <- lagDistrImpacts(T, g, P, q=as.integer(Q))
        attr(res, "Qres") <- Qres
    }
    res
}

lagImpactsExact <- function(SW, P, n) {
    direct <- sapply(P, function(x) sum(diag(x*SW))/n)
    total <- sapply(P, function(x) sum(x*SW)/n)
    indirect <- total - direct
    list(direct=direct, indirect=indirect, total=total)
}

mixedImpactsExact <- function(SW, P, n, listw) {
    p <- dim(P)[1]
    direct <- numeric(p)
    total <- numeric(p)
    W <- listw2mat(listw)
    for (i in 1:p) {
        SWr <- SW %*% (P[i,1]*diag(n) + P[i,2]*W)
        direct[i] <- sum(diag(SWr))/n
        total[i] <- sum(SWr)/n
    }
    indirect <- total - direct
    list(direct=direct, indirect=indirect, total=total)
}

processXSample <- function(x, drop2beta, type, iicept, icept, n, listw,
    irho) {
    rho <- x[irho]
    SW <- invIrW(listw, rho)
    beta <- x[-drop2beta]
    if (type == "lag" || type == "sac") {
        if (iicept) {
          P <- matrix(beta[-icept], ncol=1)
        } else {
          P <- matrix(beta, ncol=1)
        }
        return(lagImpactsExact(SW, P, n))
    } else if (type == "mixed" || type == "sacmixed") {
        if (iicept) {
            b1 <- beta[-icept]
        } else {
            b1 <- beta
        }
        p <- length(b1)
        if (p %% 2 != 0) stop("non-matched coefficient pairs")
        P <- cbind(b1[1:(p/2)], b1[((p/2)+1):p])
        return(mixedImpactsExact(SW, P, n, listw))
    }
}

intImpacts <- function(rho, beta, P, n, mu, Sigma, irho, drop2beta, bnames,
    interval, type, tr, R, listw, tol, empirical, Q, icept, iicept, p,
    mess=FALSE, samples=NULL) {
    if (is.null(listw) && is.null(tr))
        stop("either tr or listw must be given")
    timings <- list()
    .ptime_start <- proc.time()
    if (is.null(listw)) {

        q <- length(tr)-1L
        g <- rho^(0:q)
        T <- matrix(c(1, tr[-(q+1)]/n), nrow=1)
        if (type == "mixed" || type == "sacmixed") {
            T <- rbind(T, tr/n)
        }
        res <- lagImpacts(T, g, P)
        if (!is.null(Q)) {
            if (!is.numeric(Q) || length(Q) > 1L) stop("Invalid Q argument")
            if (Q > length(tr)) stop("Q larger than length of tr")
            Qres <- lagDistrImpacts(T, g, P, q=as.integer(Q))
            attr(res, "Qres") <- Qres
        }
        timings[["trace_impacts"]] <- proc.time() - .ptime_start
        .ptime_start <- proc.time()
        if (!is.null(R)) {
            if (is.null(samples)) {
                samples <- mvrnorm(n=R, mu=mu, Sigma=Sigma, tol=tol,
                    empirical=empirical)
                if (mess) samples[,irho] <- 1 - exp(samples[,irho])
            }
            if (!is.null(interval)) {
                check <- ((samples[,irho] > interval[1]) & 
                    (samples[,irho] < interval[2]))
                if (any(!check)) samples <- samples[check,]
            }
            timings[["impacts_samples"]] <- proc.time() - .ptime_start
            .ptime_start <- proc.time()
# type, iicept, icept, T, Q
            sres <- apply(samples, 1, processSample, irho=irho,
                drop2beta=drop2beta, type=type, iicept=iicept,
                icept=icept, T=T, Q=Q, q=q)
            timings[["process_samples"]] <- proc.time() - .ptime_start
            .ptime_start <- proc.time()
# 100928 Eelke Folmer
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
            ssres <- list(direct=direct, indirect=indirect, total=total)
            if (!is.null(Q)) {
                Qdirect <- as.mcmc(t(sapply(sres, function(x)
                    attr(x, "Qres")$direct)))
                Qindirect <- as.mcmc(t(sapply(sres, function(x) 
                    attr(x, "Qres")$indirect)))
                Qtotal <- as.mcmc(t(sapply(sres, function(x) 
                    attr(x, "Qres")$total)))
                Qnames <- c(sapply(bnames, function(x) 
                    paste(x, 1:Q, sep="__Q")))
                if (length(Qnames) == 1L) {
                    Qdirect <- t(Qdirect)
                    Qindirect <- t(Qindirect)
                    Qtotal <- t(Qtotal)
                }
                colnames(Qdirect) <- Qnames
                colnames(Qindirect) <- Qnames
                colnames(Qtotal) <- Qnames
                Qmcmc <- list(direct=Qdirect, indirect=Qindirect, total=Qtotal)
                attr(ssres, "Qmcmc") <- Qmcmc
            }
            timings[["postprocess_samples"]] <- proc.time() - .ptime_start
            res <- list(res=res, sres=ssres)
        }
        attr(res, "method") <- "trace"
    } else {
# added checks 140304
        stopifnot(length(listw$neighbours) == n)
        V <- listw2mat(listw)
        e <- eigen(V, only.values = TRUE)$values
        if (is.complex(e)) interval <- 1/(range(Re(e)))
	else interval <- 1/(range(e))
        SW <- invIrW(listw, rho)
        if (type == "lag" || type == "sac") res <- lagImpactsExact(SW, P, n)
        else if (type == "mixed" || type == "sacmixed")
            res <- mixedImpactsExact(SW, P, n, listw)
        timings[["weights_impacts"]] <- proc.time() - .ptime_start
        .ptime_start <- proc.time()
        if (!is.null(R)) {
            samples <- mvrnorm(n=R, mu=mu, Sigma=Sigma, tol=tol,
                empirical=empirical)
            check <- ((samples[,irho] > interval[1]) & 
                (samples[,irho] < interval[2]))
            if (any(!check)) samples <- samples[check,]
            timings[["impacts_samples"]] <- proc.time() - .ptime_start
            .ptime_start <- proc.time()
# type, iicept, icept, SW, n, listw
            sres <- apply(samples, 1, processXSample,
                drop2beta=drop2beta, type=type, iicept=iicept,
                icept=icept, n=n, listw=listw, irho=irho)
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
        attr(res, "method") <- "exact"
    }
    if (!is.null(R)) attr(res, "samples") <- list(samples=samples, irho=irho,
        drop2beta=drop2beta)
    attr(res, "type") <- type
    attr(res, "bnames") <- bnames
    attr(res, "haveQ") <- !is.null(Q)
    attr(res, "timings") <- do.call("rbind", timings)[, c(1,3)]
    class(res) <- "lagImpact"
    res
}

impacts.sarlm <- function(obj, ..., tr=NULL, R=NULL, listw=NULL, useHESS=NULL,
  tol=1e-6, empirical=FALSE, Q=NULL) {
    if (obj$type == "error") {
        if (obj$etype == "emixed") {
            return(impactSDEM(obj))
        } else {
            stop("impact measures not for error models")
        }
    }
    if (is.null(listw) && !is.null(obj$listw_style) && 
            obj$listw_style != "W")
            stop("Only row-standardised weights supported")
    rho <- obj$rho
    beta <- obj$coefficients
    s2 <- obj$s2
    if (obj$type == "sac" || obj$type == "sacmixed") lambda <- obj$lambda
    usingHESS <- NULL
    iNsert <- obj$insert
    if (!is.null(R)) {
        resvar <- obj$resvar
        usingHESS <- FALSE
        irho <- 2
        drop2beta <- 1:2
        if (obj$type == "sac" || obj$type == "sacmixed")
            drop2beta <- c(drop2beta, 3)
        if (is.logical(resvar)) {
            fdHess <- obj$fdHess
            if (is.logical(fdHess)) 
                stop("coefficient covariance matrix not available")
            usingHESS <- TRUE
            if (!iNsert) {
                irho <- 1
                drop2beta <- 1
                if (obj$type == "sac" || obj$type == "sacmixed")
                    drop2beta <- c(drop2beta, 2)
            }
        }
        if (!is.null(useHESS) && useHESS) {
            fdHess <- obj$fdHess
            if (is.logical(fdHess)) 
                stop("Hessian matrix not available")
            usingHESS <- TRUE
            if (!iNsert) {
                irho <- 1
                drop2beta <- 1
                if (obj$type == "sac" || obj$type == "sacmixed")
                    drop2beta <- c(drop2beta, 2)
            }
        }
        interval <- obj$interval
        if (is.null(interval)) interval <- c(-1,0.999)
    }
    icept <- grep("(Intercept)", names(beta))
    iicept <- length(icept) > 0L
    if (obj$type == "lag" || obj$type == "sac") {
      if (iicept) {
        P <- matrix(beta[-icept], ncol=1)
        bnames <- names(beta[-icept])
      } else {
        P <- matrix(beta, ncol=1)
        bnames <- names(beta)
      }
      p <- length(beta)
    } else if (obj$type == "mixed" || obj$type == "sacmixed") {
      if (iicept) {
        b1 <- beta[-icept]
      } else {
        b1 <- beta
      }
      p <- length(b1)
      if (p %% 2 != 0) stop("non-matched coefficient pairs")
      P <- cbind(b1[1:(p/2)], b1[((p/2)+1):p])
      bnames <- names(b1[1:(p/2)])
    }
    n <- length(obj$residuals)
    mu <- NULL
    Sigma <- NULL
    if (!is.null(R)) {
        if (usingHESS && !iNsert) {
            mu <- c(rho, beta)
            if (obj$type == "sac" || obj$type == "sacmixed")
                mu <- c(rho, lambda, beta)
            Sigma <- fdHess
        } else {
            mu <- c(s2, rho, beta)
            if (obj$type == "sac" || obj$type == "sacmixed")
                mu <- c(s2, rho, lambda, beta)
            if (usingHESS) {
                Sigma <- fdHess
            } else {
                Sigma <- resvar
            }
        }
    }
    res <- intImpacts(rho=rho, beta=beta, P=P, n=n, mu=mu, Sigma=Sigma,
        irho=irho, drop2beta=drop2beta, bnames=bnames, interval=interval,
        type=obj$type, tr=tr, R=R, listw=listw, tol=tol, empirical=empirical,
        Q=Q, icept=icept, iicept=iicept, p=p)
    attr(res, "useHESS") <- usingHESS
    attr(res, "insert") <- iNsert
    attr(res, "iClass") <- class(obj)
    res
}

lagImpactMat <- function(x, reportQ=NULL) {
    if (is.null(x$res)) {
        direct <- x$direct
        indirect <- x$indirect
        total <- x$total
    } else {
        direct <- x$res$direct
        indirect <- x$res$indirect
        total <- x$res$total
    }
    mat <- cbind(direct, indirect, total)
    colnames(mat) <- c("Direct", "Indirect", "Total")
    rownames(mat) <- attr(x, "bnames")
    if (!is.null(reportQ) && reportQ) {
        if (is.null(x$res)) {
            Qobj <- attr(x, "Qres")
        } else {
            Qobj <- attr(x$res, "Qres")
        }
        if (is.null(Qobj)) warning("No impact components to report")
        else {
# 100928 Eelke Folmer
            if (length(attr(x, "bnames")) == 1L) {
                Qobj$direct <- matrix(Qobj$direct, ncol=1)
                Qobj$indirect <- matrix(Qobj$indirect, ncol=1)
                Qobj$total <- matrix(Qobj$total, ncol=1)
            }
            colnames(Qobj$direct) <- attr(x, "bnames")
            colnames(Qobj$indirect) <- attr(x, "bnames")
            colnames(Qobj$total) <- attr(x, "bnames")
            rownames(Qobj$direct) <- paste("Q", 1:nrow(Qobj$direct), sep="")
            rownames(Qobj$indirect) <- paste("Q", 1:nrow(Qobj$indirect), sep="")
            rownames(Qobj$total) <- paste("Q", 1:nrow(Qobj$total), sep="")
            attr(mat, "Qobj") <- Qobj
        }
    }
    mat
}


print.lagImpact <- function(x, ..., reportQ=NULL) {
    mat <- lagImpactMat(x, reportQ=reportQ)
    Qobj <- attr(mat, "Qobj")
    cat("Impact measures (", attr(x, "type"), ", ", attr(x, "method"), "):\n", sep="")
    attr(mat, "Qobj") <- NULL
    print(mat)
    if (!is.null(reportQ) && reportQ) {
        if (is.null(Qobj)) warning("No impact components to report")
        else {
            cat("=================================\nImpact components\n")
            print(Qobj)
        }
    }
    invisible(x)
}

summary.lagImpact <- function(object, ..., zstats=FALSE, short=FALSE, reportQ=NULL) {
    if (is.null(object$sres)) stop("summary method unavailable")
# pass coda arguments 101006
    direct_sum <- summary(object$sres$direct, ...)
    indirect_sum <- summary(object$sres$indirect, ...)
    total_sum <- summary(object$sres$total, ...)
# 101109 Eelke Folmer
    if (length(attr(object, "bnames")) == 1L) {
        scnames <- names(direct_sum$statistics)
        qcnames <- names(direct_sum$quantiles)
        direct_sum$statistics <- matrix(direct_sum$statistics, nrow=1)
        rownames(direct_sum$statistics) <- attr(object, "bnames")[1]
        colnames(direct_sum$statistics) <- scnames
        direct_sum$quantiles <- matrix(direct_sum$quantiles, nrow=1)
        rownames(direct_sum$quantiles) <- attr(object, "bnames")[1]
        colnames(direct_sum$quantiles) <- qcnames
        indirect_sum$statistics <- matrix(indirect_sum$statistics, nrow=1)
        rownames(indirect_sum$statistics) <- attr(object, "bnames")[1]
        colnames(indirect_sum$statistics) <- scnames
        indirect_sum$quantiles <- matrix(indirect_sum$quantiles, nrow=1)
        rownames(indirect_sum$quantiles) <- attr(object, "bnames")[1]
        colnames(indirect_sum$quantiles) <- qcnames
        total_sum$statistics <- matrix(total_sum$statistics, nrow=1)
        rownames(total_sum$statistics) <- attr(object, "bnames")[1]
        colnames(total_sum$statistics) <- scnames
        total_sum$quantiles <- matrix(total_sum$quantiles, nrow=1)
        rownames(total_sum$quantiles) <- attr(object, "bnames")[1]
        colnames(total_sum$quantiles) <- qcnames
    }
    Qmcmc <- NULL
    if (!is.null(attr(object$sres, "Qmcmc")) && !is.null(reportQ) && reportQ) {
        Qdirect_sum <- summary(attr(object$sres, "Qmcmc")$direct, ...)
        Qindirect_sum <- summary(attr(object$sres, "Qmcmc")$indirect, ...)
        Qtotal_sum <- summary(attr(object$sres, "Qmcmc")$total, ...)
        Qmcmc <- list(Qdirect_sum=Qdirect_sum, Qindirect_sum=Qindirect_sum,
           Qtotal_sum=Qtotal_sum)
    }
    lres <- list(direct_sum=direct_sum, indirect_sum=indirect_sum,
        total_sum=total_sum)
    res <- c(object, lres, Qmcmc)
    if (zstats) {
# 100928 Eelke Folmer
        if (length(attr(object, "bnames")) == 1L) {
            zmat <- sapply(lres, function(x) x$statistics[1]/x$statistics[2])
            zmat <- matrix(zmat, ncol=3)
            colnames(zmat) <- c("Direct", "Indirect", "Total")
        } else {
            zmat <- sapply(lres, function(x) x$statistics[,1]/x$statistics[,2])
            colnames(zmat) <- c("Direct", "Indirect", "Total")
        }
        pzmat <- 2*(1-pnorm(abs(zmat)))
        res <- c(res, list(zmat=zmat, pzmat=pzmat))
        if (!is.null(Qmcmc) && !is.null(reportQ) && reportQ) {
            Qzmats <- lapply(Qmcmc, function(x) {
                Qm <- matrix(x$statistics[,1]/x$statistics[,2],
                    ncol=length(attr(object, "bnames")))
                colnames(Qm) <- attr(object, "bnames")
                rownames(Qm) <- paste("Q", 1:nrow(Qm), sep="")
                Qm
            })
            names(Qzmats) <- c("Direct", "Indirect", "Total")
            Qpzmats <- lapply(Qzmats, function(x) {
                xo <- 2*(1-pnorm(abs(x)))
                rownames(xo) <- paste("Q", 1:nrow(xo), sep="")
                xo
            })
            res <- c(res, list(Qzmats=Qzmats, Qpzmats=Qpzmats))
        }
    }
    attr(res, "useHESS") <- attr(object, "useHESS")
    attr(res, "bnames") <- attr(object, "bnames")
    attr(res, "method") <- attr(object, "method")
    attr(res, "insert") <- attr(object, "insert")
    attr(res, "type") <- attr(object, "type")
    attr(res, "short") <- short
    attr(res, "reportQ") <- reportQ
    tp <- NULL
    if ("sarlm" %in% attr(object, "iClass")) {
       tp <- ifelse(attr(object,
           "useHESS"), ifelse(attr(object, "insert"),
           "mixed Hessian approximation", "numerical Hessian approximation"),
           "asymptotic")
    } else if ("lagmess" %in% attr(object, "iClass")) {
        tp <- "numerical Hessian approximation"
    } else if ("stsls" %in% attr(object, "iClass")) {
        tp <- "asymptotic IV"
        if (!is.null(attr(object, "robust")) && attr(object, "robust")) {
            HC <- attr(object, "HC")
            if (is.null(HC)) HC <- "HC0"
            tp <- paste(HC, "IV")
        }
    }
    if ("sphet" %in% attr(object, "iClass")) {
            tp <- "IV HAC"
            if ("gstsls" %in% attr(object, "iClass")) 
                tp <- "GSTSLS"
    }
    if ("MCMC_sar_g" %in% attr(object, "iClass")) tp <- "MCMC samples"
    attr(res, "tp") <- tp
    class(res) <- "summary.lagImpact"
    res
}

print.summary.lagImpact <- function(x, ...) {
    reportQ <- attr(x, "reportQ")
    mat <- lagImpactMat(x, reportQ)
    Qobj <- attr(mat, "Qobj")
    attr(mat, "Qobj") <- NULL
    cat("Impact measures (", attr(x, "type"), ", ", attr(x, "method"),
        "):\n", sep="")
    print(mat)
    if (!is.null(reportQ) && reportQ) {
        if (is.null(Qobj)) warning("No impact components to report")
        else {
            cat("=================================\nImpact components\n")
            print(Qobj)
        }
    }
    cat("========================================================\n")

    if (!is.null(attr(x, "tp")) && attr(x, "tp") == "MCMC samples") {
        cat("MCMC sample impact results:\n")
    } else {
        cat("Simulation results (", attr(x, "tp"), " variance matrix):\n",
            sep="")
    }
    if (!attr(x, "short")) {
        cat("Direct:\n")
        print(x$direct_sum)
        cat("========================================================\n")
        cat("Indirect:\n")
        print(x$indirect_sum)
        cat("========================================================\n")
        cat("Total:\n")
        print(x$total_sum)
        if (!is.null(reportQ) && reportQ && !is.null(x$Qdirect_sum)) {
            cat("========================================================\n")
            cat("Direct impact components:\n")
            print(x$Qdirect_sum)
            cat("========================================================\n")
            cat("Indirect impact components:\n")
            print(x$Qindirect_sum)
            cat("========================================================\n")
            cat("Total impact components:\n")
            print(x$Qtotal_sum)
        }
    }
    if (!is.null(x$zmat)) {
        cat("========================================================\n")
        cat("Simulated z-values:\n")
        mat <- x$zmat
        rownames(mat) <- attr(x, "bnames")
        print(mat)
        cat("\nSimulated p-values:\n")
        xx <- apply(x$pzmat, 2, format.pval)
# 100928 Eelke Folmer
        if (length(attr(x, "bnames")) == 1L) {
            xx <- matrix(xx, ncol=3)
            colnames(xx) <- c("Direct", "Indirect", "Total")
        }
        rownames(xx) <- attr(x, "bnames")
        print(xx, quote=FALSE)
        if (!is.null(x$Qzmats)) {
            cat("========================================================\n")
            cat("Simulated impact components z-values:\n")
            print(x$Qzmats)
            cat("\nSimulated impact components p-values:\n")
            xx <- lapply(x$Qpzmats, function(y) {
                xo <- apply(y, 2, format.pval)
                rownames(xo) <- paste("Q", 1:nrow(xo), sep="")
                xo
            })
            print(xx, quote=FALSE)
        }
    }
    invisible(x)
}

plot.lagImpact <- function(x, ..., choice="direct", trace=FALSE,
    density=TRUE) {
    if (is.null(x$sres)) stop("plot method unavailable")
    plot(x$sres[[choice]], trace=trace, density=density, sub=choice)
    invisible(x)
}

HPDinterval.lagImpact <- function(obj, prob = 0.95, ..., choice="direct") {
    if (is.null(obj$sres)) stop("HPDinterval method unavailable")
    res <- HPDinterval(obj$sres[[choice]], prob=prob)
    res
}

