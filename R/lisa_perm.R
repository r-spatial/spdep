parallel_setup <- function(iseed) {
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
    list(parallel=parallel, ncpus=ncpus, cl=cl)
}

spdep_splitIndices <- function(idx, lcl) {
    n <- length(idx)
    stopifnot(n == length(unique(idx)))
    splits <- parallel::splitIndices(n, lcl)
    res <- lapply(splits, function(i) idx[i])
    res
}

run_perm <- function(fun, idx, env, iseed, varlist) {
    p_setup <- parallel_setup(iseed)
    parallel <- p_setup$parallel
    ncpus <- p_setup$ncpus
    cl <- p_setup$cl
    if (parallel == "snow") {
      if (requireNamespace("parallel", quietly = TRUE)) {
        sI <- spdep_splitIndices(idx, length(cl))
        parallel::clusterExport(cl, varlist=varlist, envir=env)
        if (!is.null(iseed)) parallel::clusterSetRNGStream(cl, iseed = iseed)
        oo <- parallel::clusterApply(cl, x = sI, fun=lapply, function(i) {
 	    fun(i, env)})
        out <- do.call("rbind", do.call("c", oo))
      } else {
        stop("parallel not available")
      }
    } else if (parallel == "multicore") {
      if (requireNamespace("parallel", quietly = TRUE)) {
        sI <- spdep_splitIndices(idx, ncpus)
        oldRNG <- RNGkind()
        RNGkind("L'Ecuyer-CMRG")
        if (!is.null(iseed)) set.seed(iseed)
        oo <- parallel::mclapply(sI, FUN=lapply, function(i) {fun(i,
            env)}, mc.cores=ncpus)
        RNGkind(oldRNG[1])
        out <- do.call("rbind", do.call("c", oo))
      } else {
        stop("parallel not available")
      }
    } else {
        if (!is.null(iseed)) set.seed(iseed)
        oo <- lapply(idx, function(i) fun(i, env))
        out <- do.call("rbind", oo)
    }
    out
}

probs_lut <- function(stat="I", nsim, alternative) {
    gr <- punif((1:(nsim+1))/(nsim+1), 0, 1)
    ls <- rev(gr)
    ts <- (ifelse(gr > ls, ls, gr))*2
    if (alternative == "two.sided") {
        probs <- ts
        Prname <- paste0("Pr(z != E(", stat, "i))")
    } else if (alternative == "greater") {
        Prname <- paste0("Pr(z > E(", stat, "i))")
        probs <- gr
    } else {
        Prname <- paste0("Pr(z < E(", stat, "i))")
        probs <- ls
    }
    attr(probs, "Prname") <- Prname
    probs
}

localmoran_perm <- function(x, listw, nsim=499L, zero.policy=attr(listw, "zero.policy"),
    na.action=na.fail, alternative = "two.sided",
    mlvar=TRUE, spChk=NULL, adjust.x=FALSE, sample_Ei=TRUE, iseed=NULL,
    no_repeat_in_row=FALSE) {
    alternative <- match.arg(alternative, c("greater", "less", "two.sided"))
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
    if (n != length(x)) stop("Different numbers of observations")
    res <- matrix(nrow=n, ncol=9)
    if (adjust.x) {
        nc <- card(listw$neighbours) > 0L
	xx <- mean(x[nc], na.rm=NAOK)
    } else {
        xx <- mean(x, na.rm=NAOK)
    }
    z <- x - xx 
    EIc <- -(z^2 * sapply(listw$weights, sum)) / ((n - 1) * (sum(z * z) / n))
    lz <- lag.listw(listw, z, zero.policy=zero.policy, NAOK=NAOK)
    lbs <- c("Low", "High")
# https://github.com/pysal/esda/blob/4a63e0b5df1e754b17b5f1205b8cadcbecc5e061/esda/moran.py#L1068-L1081
    quadr_ps <- interaction(cut(z, c(-Inf, 0, Inf), labels=lbs), 
        cut(lz, c(-Inf, 0, Inf), labels=lbs), sep="-")
    lx <- lag.listw(listw, x, zero.policy=zero.policy, NAOK=NAOK)
    lxx <- mean(lx, na.rm=NAOK)
    quadr <- interaction(cut(x, c(-Inf, xx, Inf), labels=lbs), 
        cut(lx, c(-Inf, lxx, Inf), labels=lbs), sep="-")
    xmed <- median(x, na.rm=NAOK)
    lxmed <- median(lx, na.rm=NAOK)
    quadr_med <- interaction(cut(x, c(-Inf, xmed, Inf), labels=lbs),
        cut(lx, c(-Inf, lxmed, Inf), labels=lbs), sep="-")

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

    crd <- card(listw$neighbours)
    lww <- listw$weights
    Iis <- res[,1]

    env <- new.env()
    assign("z", z, envir=env)
    assign("crd", crd, envir=env)
    assign("lww", lww, envir=env)
    assign("nsim", nsim, envir=env)
    assign("Iis", Iis, envir=env)
    assign("s2", s2, envir=env)
    assign("n", n, envir=env)
    assign("no_repeat_in_row", no_repeat_in_row, envir=env)
    varlist <- ls(envir = env)

    permI_int <- function(i, env) {
        res_i <- rep(as.numeric(NA), 6) # initialize output
        crdi <- get("crd", envir=env)[i]
        no_repeat_in_row <- get("no_repeat_in_row", envir=env)
        if (crdi > 0) { # if i has neighbours
            nsim <- get("nsim", envir=env)
            n_i <- get("n", envir=env) - 1L
            zi <- get("z", envir=env)[i]
            z_i <- get("z", envir=env)[-i]
            if (no_repeat_in_row) {
              samples <- .Call("perm_no_replace", as.integer(nsim),
                as.integer(n_i), as.integer(crdi), PACKAGE="spdep")
              sz_i <- matrix(z_i[samples], ncol=crdi, nrow=nsim)
            } else {
              sz_i <- matrix(sample(z_i, size=crdi*nsim, replace=TRUE),
                ncol=crdi, nrow=nsim) # permute nsim*#neighbours from z[-i]
            }
            wtsi <- get("lww", envir=env)[[i]]
            lz_i <- sz_i %*% wtsi # nsim by 1 = nsim by crdi %*% crdi by 1
            # create nsim samples of Ii at i
            s2 <- get("s2", envir=env)
            res_p <- (zi/s2)*lz_i # nsim by 1 = scalar/scalar * nsim by 1
            res_i[1] <- mean(res_p)
            res_i[2] <- var(res_p)
            Ii <- get("Iis", envir=env)[i]
            xrank <- rank(c(res_p, Ii))[(nsim + 1L)]
	    res_i[3] <- xrank
            rnk0 <- as.integer(sum(res_p >= Ii))
            drnk0 <- nsim - rnk0
            rnk <- ifelse(drnk0 < rnk0, drnk0, rnk0)
            res_i[4] <- rnk0
            res_i[5] <- e1071::skewness(res_p)
            res_i[6] <- e1071::kurtosis(res_p)
        }
        res_i
    }

    out <- run_perm(fun=permI_int, idx=1:n, env=env, iseed=iseed, varlist=varlist)

    if (sample_Ei) res[,2] <- out[,1]
    else  res[,2] <- EIc
    res[,3] <- out[,2]
    res[,4] <- (res[,1] - res[,2])/sqrt(res[,3])
    if (alternative == "two.sided") 
        res[,5] <- 2 * pnorm(abs(res[,4]), lower.tail=FALSE)
    else if (alternative == "greater") 
        res[,5] <- pnorm(res[,4], lower.tail=FALSE)
    else res[,5] <- pnorm(res[,4])
# look-up table
    probs <- probs_lut(stat="I", nsim=nsim, alternative=alternative)
    Prname <- attr(probs, "Prname")
    Prname_rank <- paste0(Prname, " Sim")
    Prname_sim <- "Pr(folded) Sim"
    res[,6] <- probs[as.integer(out[,3])]
# 210811 from https://github.com/pysal/esda/blob/4a63e0b5df1e754b17b5f1205b8cadcbecc5e061/esda/crand.py#L211-L213
    rnk0 <- as.integer(out[,4])
    drnk0 <- nsim - rnk0
    rnk <- ifelse(drnk0 < rnk0, drnk0, rnk0)
# folded
    res[,7] <- (rnk + 1.0) / (nsim + 1.0)
# skewness
    res[,8] <- out[,5]
# kurtosis
    res[,9] <- out[,6]
    colnames(res) <- c("Ii", "E.Ii", "Var.Ii", "Z.Ii", Prname, Prname_rank, 
        Prname_sim, "Skewness", "Kurtosis")
    if (!is.null(na.act) && excl) {
	res <- naresid(na.act, res)
    }
    if (!is.null(rn)) rownames(res) <- rn
    attr(res, "call") <- match.call()
    if (!is.null(na.act)) attr(res, "na.action") <- na.act
    attr(res, "quadr") <- data.frame(mean=quadr, median=quadr_med, 
        pysal=quadr_ps)
    class(res) <- c("localmoran", class(res))
    res
}

# "localmoran" quadr mean/median/pysal "Low-Low", "Low-High", "High-Low", "High-High"


localG_perm <- function(x, listw, nsim=499, zero.policy=attr(listw, "zero.policy"), spChk=NULL, alternative = "two.sided", iseed=NULL, fix_i_in_Gstar_permutations=TRUE, no_repeat_in_row=FALSE) {
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

# #124 repeat localG standard deviate direct output
    if (gstari) {
	xibar <- rep(mean(x), n)
	si2 <- rep(sum(scale(x, scale=FALSE)^2)/n, n)
    } else {
	xibar <- (rep(sum(x),n) - x) / (n - 1)
	si2 <- ((rep(sum(x^2),n) - x^2) / (n - 1)) - xibar^2
    }
    Wi <- sapply(listw$weights, sum)
    S1i <- sapply(listw$weights, function(x) sum(x^2))
    EG <- Wi*xibar
    if (gstari) {
        res <- G - (EG/x_star)
    } else {
        res <- G - (EG/(x_star-c(x)))
    }
    if (gstari) {
        VG <- si2*((n*S1i - Wi^2)/(n-1))
    } else {
        VG <- si2*(((n-1)*S1i - Wi^2)/(n-2))
    }
    if (gstari) {
        res <- res / sqrt(VG/x_star^2)
    } else {
        res <- res / sqrt(VG/(x_star-c(x))^2)
    }

    crd <- card(listw$neighbours)
    lww <- listw$weights

    env <- new.env()
    assign("x", x, envir=env)
    assign("crd", crd, envir=env)
    assign("lww", lww, envir=env)
    assign("nsim", nsim, envir=env)
    assign("G", G, envir=env)
    assign("x_star", x_star, envir=env)
    assign("gstari", gstari, envir=env)
    assign("n", n, envir=env)
    assign("no_repeat_in_row", no_repeat_in_row, envir=env)
    varlist <- ls(envir = env)

    permG_int <- function(i, env) {
        res_i <- rep(as.numeric(NA), 6)
        crdi <- get("crd", envir=env)[i]
        n_i <- get("n", envir=env) - 1L
        no_repeat_in_row <- get("no_repeat_in_row", envir=env)
        if (crdi > 0) { # if i has neighbours
            nsim <- get("nsim", envir=env)
            xi <- get("x", envir=env)[i]
            x_i <- get("x", envir=env)[-i]
            if (no_repeat_in_row) {
              samples <- .Call("perm_no_replace", as.integer(nsim),
                as.integer(n_i), as.integer(crdi), PACKAGE="spdep")
              sx_i <- matrix(x_i[samples], ncol=crdi, nrow=nsim)
            } else {
              sx_i <- matrix(sample(x_i, size=crdi*nsim, replace=TRUE),
                ncol=crdi, nrow=nsim) # permute nsim*#neighbours from x[-i]
            }
            wtsi <- get("lww", envir=env)[[i]]
            lx_i <- sx_i %*% wtsi # nsim by 1 = nsim by crdi %*% crdi by 1
            # create nsim samples of Gi at i
            x_star <- get("x_star", envir=env)
            gstari <- get("gstari", envir=env)
            # nsim by 1 = nsim by 1 / scalar
            if (gstari) res_p <- lx_i/x_star
            else res_p <- lx_i/(x_star-xi)
            res_i[1] <- mean(res_p)
            res_i[2] <- var(res_p)
            Gi <- get("G", envir=env)[i]
	    res_i[3] <- rank(c(res_p, Gi))[(nsim + 1L)]
            res_i[4] <- as.integer(sum(res_p >= Gi))
            res_i[5] <- e1071::skewness(res_p)
            res_i[6] <- e1071::kurtosis(res_p)
        }
        res_i
    }

    permGstar_int <- function(i, env) {
        res_i <- rep(as.numeric(NA), 6)
        crdi <- get("crd", envir=env)[i]
        n_i <- get("n", envir=env) - 1L
        no_repeat_in_row <- get("no_repeat_in_row", envir=env)
        if (crdi > 0) { # if i has neighbours
            nsim <- get("nsim", envir=env)
            xi <- get("x", envir=env)[i]
            x_i <- get("x", envir=env)[-i]
            if (no_repeat_in_row) {
              samples <- .Call("perm_no_replace", as.integer(nsim),
                as.integer(n_i), as.integer(crdi), PACKAGE="spdep")
              sx_i <- matrix(x_i[samples], ncol=crdi, nrow=nsim)
            } else {
              sx_i <- matrix(sample(x_i, size=crdi*nsim, replace=TRUE),
                ncol=crdi, nrow=nsim) # permute nsim*#neighbours from x[-i]
            }
            wtsi <- get("lww", envir=env)[[i]]
            nbsi <- get("nbs", envir=env)[[i]]
            ithnb <- which(nbsi == i)
            sx_i[, ithnb] <- xi
            lx_i <- sx_i %*% wtsi # nsim by 1 = nsim by crdi %*% crdi by 1
            # create nsim samples of Gi at i
            x_star <- get("x_star", envir=env)
            # nsim by 1 = nsim by 1 / scalar
            res_p <- lx_i/x_star
            res_i[1] <- mean(res_p)
            res_i[2] <- var(res_p)
            Gi <- get("G", envir=env)[i]
	    res_i[3] <- rank(c(res_p, Gi))[(nsim + 1L)]
            res_i[4] <- as.integer(sum(res_p >= Gi))
            res_i[5] <- e1071::skewness(res_p)
            res_i[6] <- e1071::kurtosis(res_p)
        }
        res_i
    }

    thisfun <- permG_int
    if (gstari && fix_i_in_Gstar_permutations) {
        thisfun <- permGstar_int
        assign("nbs", listw$neighbours, envir=env)
    }

    out <- run_perm(fun=thisfun, idx=1:n, env=env, iseed=iseed, varlist=varlist)

# add simulated standard deviate direct output
    res_sim <- (G - out[,1])
    res_sim <- res_sim / sqrt(out[,2])
    probs <- probs_lut(stat="G", nsim=nsim, alternative=alternative)
    Prname <- attr(probs, "Prname")
    Prname_rank <- paste0(Prname, " Sim")
    Prname_sim <- "Pr(folded) Sim"
    if (alternative == "two.sided") 
        pv <- 2 * pnorm(abs(res_sim), lower.tail=FALSE)
    else if (alternative == "greater") 
        pv <- pnorm(res_sim, lower.tail=FALSE)
    else pv <- pnorm(res_sim)
# look-up table
    pu <- probs[as.integer(out[,3])]
# 210811 from https://github.com/pysal/esda/blob/4a63e0b5df1e754b17b5f1205b8cadcbecc5e061/esda/crand.py#L211-L213
    rnk0 <- as.integer(out[,4])
    drnk0 <- nsim - rnk0
    rnk <- ifelse(drnk0 < rnk0, drnk0, rnk0)
# folded
    pr <- (rnk + 1.0) / (nsim + 1.0)
    resint <- cbind(G=G, EG_sim=out[,1], VG_sim=out[,2], res_sim=res_sim,
        pv=pv, pu=pu, pr=pr, sk=out[,5], ku=out[,6])
    colnames(resint) <- c("Gi", "E.Gi", "Var.Gi", "StdDev.Gi", Prname,
        Prname_rank, Prname_sim, "Skewness", "Kurtosis")
    attr(res, "internals") <- resint
    attr(res, "cluster") <- cut(x, c(-Inf, mean(x), Inf), labels = c("Low", "High"))
    attr(res, "gstari") <- gstari
    attr(res, "call") <- match.call()
    class(res) <- "localG"
    res
}

# "localG" cluster c("Low", "High")
