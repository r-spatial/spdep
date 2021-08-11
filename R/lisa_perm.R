localmoran_perm <- function(x, listw, nsim=499L, zero.policy=NULL,
    na.action=na.fail, alternative = "greater", p.adjust.method="none",
    mlvar=TRUE, spChk=NULL, adjust.x=FALSE, sample_Ei=TRUE, iseed=NULL) {
    stopifnot(is.vector(x))
    if (!inherits(listw, "listw"))
	stop(paste(deparse(substitute(listw)), "is not a listw object"))
    if (is.null(zero.policy))
        zero.policy <- get("zeroPolicy", envir = .spdepOptions)
    stopifnot(is.logical(zero.policy))
    if (!is.null(attr(listw$neighbours, "self.included")) &&
	attr(listw$neighbours, "self.included"))
	stop("Self included among neighbours")
    if (is.null(spChk)) spChk <- get.spChkOption()
    if (spChk && !chkIDs(x, listw))
	stop("Check of data and weights ID integrity failed")
    if (!is.numeric(x))
	stop(paste(deparse(substitute(x)), "is not a numeric vector"))
    NAOK <- deparse(substitute(na.action)) == "na.pass"
    x <- na.action(x)
    na.act <- attr(x, "na.action")
    rn <- attr(listw, "region.id")
    if (!is.null(na.act)) {
        subset <- !(1:length(listw$neighbours) %in% na.act)
        listw <- subset(listw, subset, zero.policy=zero.policy)
        excl <- class(na.act) == "exclude"
    }
    n <- length(listw$neighbours)
    if (n != length(x))stop("Different numbers of observations")
    res <- matrix(nrow=n, ncol=6)
    if (alternative == "two.sided") Prname <- "Pr(z != 0)"
    else if (alternative == "greater") Prname <- "Pr(z > 0)"
    else Prname <- "Pr(z < 0)"
    colnames(res) <- c("Ii", "E.Ii", "Var.Ii", "Z.Ii", Prname, "p_sim")
    if (adjust.x) {
        nc <- card(listw$neighbours) > 0L
	xx <- mean(x[nc], na.rm=NAOK)
    } else {
        xx <- mean(x, na.rm=NAOK)
    }
    z <- x - xx 
    lz <- lag.listw(listw, z, zero.policy=zero.policy, NAOK=NAOK)

    if (mlvar) {
        if (adjust.x) {
            s2 <- sum(z[nc]^2, na.rm=NAOK)/sum(nc)
        } else {
            s2 <- sum(z^2, na.rm=NAOK)/n
        }
    } else {
        if (adjust.x) {
            s2 <- sum(z[nc]^2, na.rm=NAOK)/(sum(nc)-1) 
        } else {
            s2 <- sum(z^2, na.rm=NAOK)/(n-1) 
        }
    }
    res[,1] <- (z/s2) * lz

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
    if (!is.null(iseed)) {
        stopifnot(is.numeric(iseed))
        stopifnot(length(iseed) == 1L)
    }

    crd <- card(listw$neighbours)
    permI_int <- function(i, zi, z_i, crdi, wtsi, nsim, Ii) {
        if (crdi > 0) {
            sz_i <- matrix(sample(z_i, size=crdi*nsim, replace=TRUE),
                ncol=crdi, nrow=nsim)
            lz_i <- sz_i %*% wtsi
            res_p <- (zi/s2)*lz_i
            res <- c(mean(res_p), var(res_p), as.integer(sum(res_p >= Ii)))
        } else {
            res <- c(as.numeric(NA), as.numeric(NA), as.integer(NA))
        }
        res
    }

    lww <- listw$weights
    Iis <- res[,1]

    if (parallel == "snow") {
      if (requireNamespace("parallel", quietly = TRUE)) {
        sI <- parallel::splitIndices(n, length(cl))
        env <- new.env()
        assign("z", z, envir=env)
        assign("crd", crd, envir=env)
        assign("lww", lww, envir=env)
        assign("nsim", nsim, envir=env)
        assign("Iis", Iis, envir=env)
        parallel::clusterExport(cl, varlist=c("z", "crd", "lww", "nsim", "Iis"),
            envir=env)
        if (!is.null(iseed)) parallel::clusterSetRNGStream(cl, iseed = iseed)
        oo <- parallel::clusterApply(cl, x = sI, fun=lapply, function(i) {
 	    permI_int(i, z[i], z[-i], crd[i], lww[[i]], nsim, Iis[i])})
        out <- do.call("rbind", do.call("c", oo))
        rm(env)
      } else {
        stop("parallel not available")
      }
    } else if (parallel == "multicore") {
      if (requireNamespace("parallel", quietly = TRUE)) {
        sI <- parallel::splitIndices(n, ncpus)
        oldRNG <- RNGkind()
        RNGkind("L'Ecuyer-CMRG")
        if (!is.null(iseed)) set.seed(iseed)
        oo <- parallel::mclapply(sI, FUN=lapply, function(i) {permI_int(i,
            z[i], z[-i], crd[i], lww[[i]], nsim, Iis[i])}, mc.cores=ncpus)
        RNGkind(oldRNG[1])
        out <- do.call("rbind", do.call("c", oo))
      } else {
        stop("parallel not available")
      }
    } else {
        if (!is.null(iseed)) set.seed(iseed)
        oo <- lapply(1:n, function(i) permI_int(i, z[i], z[-i], 
            crd[i], lww[[i]], nsim, Iis[i]))
        out <- do.call("rbind", oo)
    }

    if (!sample_Ei) {
        res[,2] <- -sapply(listw$weights, sum) / (n-1)
    } else {
        res[,2] <- out[,1]
    }
    res[,3] <- out[,2]

    res[,4] <- (res[,1] - res[,2]) / sqrt(res[,3])
    if (alternative == "two.sided") pv <- 2 * pnorm(abs(res[,4]), 
        lower.tail=FALSE)
    else if (alternative == "greater")
        pv <- pnorm(res[,4], lower.tail=FALSE)
    else pv <- pnorm(res[,4])
    res[,5] <- p.adjustSP(pv, listw$neighbours, method=p.adjust.method)
    low_extreme <- (nsim - out[,3]) < out[,3]
    out[low_extreme, 3] <- nsim - out[low_extreme, 3]
    res[,6] <- (out[,3] + 1.0) / (nsim + 1.0)
    if (!is.null(na.act) && excl) {
	res <- naresid(na.act, res)
    }
    if (!is.null(rn)) rownames(res) <- rn
    attr(res, "call") <- match.call()
    if (!is.null(na.act)) attr(res, "na.action") <- na.act
    class(res) <- c("localmoran", class(res))
    res
}



localG_perm <- function(x, listw, nsim=499, zero.policy=NULL, spChk=NULL, return_internals=FALSE, iseed=NULL) {
    if (!inherits(listw, "listw"))
	stop(paste(deparse(substitute(listw)), "is not a listw object"))
    if (!is.numeric(x))
	stop(paste(deparse(substitute(x)), "is not a numeric vector"))
    if (is.null(zero.policy))
        zero.policy <- get("zeroPolicy", envir = .spdepOptions)
    stopifnot(is.logical(zero.policy))
    stopifnot(is.vector(x))
    if (any(is.na(x))) stop(paste("NA in ", deparse(substitute(x))))
    n <- length(listw$neighbours)
    if (n != length(x))stop("Different numbers of observations")
    if (is.null(spChk)) spChk <- get.spChkOption()
    if (spChk && !chkIDs(x, listw))
	stop("Check of data and weights ID integrity failed")
    gstari <- FALSE
    if (!is.null(attr(listw$neighbours, "self.included")) &&
	attr(listw$neighbours, "self.included")) gstari <- TRUE
    lx <- lag.listw(listw, x, zero.policy=zero.policy)
    x_star <- sum(x)
    if (gstari) {
        G <- lx/x_star
    } else {
        G <- lx/(x_star-c(x))
    }

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
    if (!is.null(iseed)) {
        stopifnot(is.numeric(iseed))
        stopifnot(length(iseed) == 1L)
    }

    permG_int <- function(i, xi, x_i, crdi, wtsi, nsim) {
        if (crdi > 0) {
            sx_i <- matrix(sample(x_i, size=crdi*nsim, replace=TRUE),
                ncol=crdi, nrow=nsim)
            lx_i <- sx_i %*% wtsi
            res_p <- lx_i/(x_star-xi)
            res <- c(mean(res_p), var(res_p))
        } else {
            res <- c(as.numeric(NA), as.numeric(NA))
        }
        res
    }

    crd <- card(listw$neighbours)
    lww <- listw$weights

    if (parallel == "snow") {
      if (requireNamespace("parallel", quietly = TRUE)) {
        sI <- parallel::splitIndices(n, length(cl))
        env <- new.env()
        assign("x", x, envir=env)
        assign("crd", crd, envir=env)
        assign("lww", lww, envir=env)
        assign("nsim", nsim, envir=env)
        parallel::clusterExport(cl, varlist=c("x", "crd", "lww", "nsim"),
            envir=env)
        if (!is.null(iseed)) parallel::clusterSetRNGStream(cl, iseed = iseed)
        oo <- parallel::clusterApply(cl, x = sI, fun=lapply, function(i) {
 	    permG_int(i, x[i], x[-i], crd[i], lww[[i]], nsim)})
        out <- do.call("rbind", do.call("c", oo))
        rm(env)
      } else {
        stop("parallel not available")
      }
    } else if (parallel == "multicore") {
      if (requireNamespace("parallel", quietly = TRUE)) {
        sI <- parallel::splitIndices(n, ncpus)
        oldRNG <- RNGkind()
        RNGkind("L'Ecuyer-CMRG")
        if (!is.null(iseed)) set.seed(iseed)
        oo <- parallel::mclapply(sI, FUN=lapply, function(i) {permG_int(i,
            x[i], x[-i], crd[i], lww[[i]], nsim)}, mc.cores=ncpus)
        RNGkind(oldRNG[1])
        out <- do.call("rbind", do.call("c", oo))
      } else {
        stop("parallel not available")
      }
    } else {
        if (!is.null(iseed)) set.seed(iseed)
        oo <- lapply(1:n, function(i) permG_int(i, x[i], x[-i], 
            crd[i], lww[[i]], nsim))
        out <- do.call("rbind", oo)
    }
    EG <- out[,1]
    VG <- out[,2]

    res <- (G - EG)
    res <- res / sqrt(VG)
    if (return_internals) {
        attr(res, "internals") <- cbind(G=G, EG=EG, VG=VG)
    }
    attr(res, "gstari") <- gstari
    attr(res, "call") <- match.call()
    class(res) <- "localG"
    res
}
