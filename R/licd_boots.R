local_jcm_BW <- function(res, tab, wc) {
    ldiag <- res[2,1] + res[1,2]
    ntab <- as.numeric(as.vector(tab))
    Ejc <- (wc$S0*(ntab*(ntab-1))) / (2*wc$nwcn1)
    Vjc <- (wc$S1*(ntab*(ntab-1))) / (wc$nwcn1)
    Vjc <- Vjc + (((wc$S2 - 2*wc$S1)*ntab*(ntab-1)*(ntab-2)) /
        (wc$nwcn2))
    Vjc <- Vjc + (((wc$S02 + wc$S1 - wc$S2)*ntab*(ntab-1)*(ntab-2)*
        (ntab-3)) / (wc$nwcn2*wc$n3))
    Vjc <- (0.25 * Vjc) - Ejc^2
    nrns <- ntab[2]*ntab[1]
    pntab <- (ntab*(ntab-1))
    nrns1 <- pntab[2]*pntab[1]
    Exp <- (wc$S0*(nrns)) / (wc$nwcn1)
    Var <- (2*wc$S1*nrns)/(wc$nwcn1)
    Var <- Var + (((wc$S2 - 2*wc$S1) * nrns * (ntab[2]+ntab[1]-2))/ 
        (wc$nwcn2))
    Var <- Var + ((4*(wc$S02 + wc$S1 - wc$S2)*nrns1) / (wc$nwcn2*wc$n3))
    Var <- (0.25 * Var) - Exp^2
    JtotExp <- sum(Exp)
    O_jc_BW <- c(diag(res), ldiag)
    E_jc_BW <- c(Ejc, Exp)
    V_jc_BW <- c(Vjc, Var)
    jc_conf_chi_BW <- local_chi(O_jc_BW, E_jc_BW)
    list(jc_conf_chi_BW, O_jc_BW, E_jc_BW, V_jc_BW)
}
local_chi <- function(O, E) {
    sum((O - E)^2 / E, na.rm=TRUE)
}
local_anscombe <- function(s, n) {
    asin(sqrt((s+(3/8))/(n+(3/4))))
}

licd_multi <- function(fx, listw, zero.policy=attr(listw, "zero.policy"),
    adjust.n=TRUE, nsim = 0L, iseed = NULL, no_repeat_in_row=FALSE,
    control=list()) {
    timings <- list()
    .ptime_start <- proc.time()
    con <- list(comp_binary=TRUE, binomial_punif_alternative="greater",
        jcm_same_punif_alternative="less",
        jcm_diff_punif_alternative="greater",
        uni_jc_same_punif_alternative="two.sided",
        rank_ties.method="min", unique_ceiling=1/3, check_reps=FALSE,
        pysal_rank=FALSE, pysal_sim_obs="GT", na.last="keep", xtras=FALSE)
    nmsC <- names(con)
    con[(namc <- names(control))] <- control
    if (length(noNms <- namc[!namc %in% nmsC])) 
        warning("unknown names in control: ", paste(noNms, collapse = ", "))
    if (length(con$rank_ties.method) != 1L) stop("control rank_ties.method invalid")
    if (!is.character(con$rank_ties.method))
        stop("control rank_ties.method invalid")
    if (is.na(match(con$rank_ties.method, 
         c("average", "first", "last", "random", "max", "min")))) {
         stop("control rank_ties.method string invalid: ", con$rank_ties.method)
    }
    if (length(con$na.last) != 1L) stop("control na.last invalid")
    if (!is.na(con$na.last)) {
        if (is.character(con$na.last)) {
            if (is.na(match(con$na.last, "keep")))
                stop("control na.last string invalid: ", con$na.last)
        } else if (!is.logical(con$na.last)) {
            stop("control na.last invalid")
        }
    }
    if(!inherits(listw, "listw")) stop(paste(deparse(substitute(listw)),
        "is not a listw object"))
    if(!is.factor(fx)) stop(paste(deparse(substitute(fx)),
        "is not a factor"))
    if (any(is.na(fx))) stop("NA in factor")
    if (inherits(fx, "ordered")) warning("use of joincount tests on ordinal variables is not well understood")
    if (is.null(zero.policy))
        zero.policy <- get.ZeroPolicyOption()
    stopifnot(is.logical(zero.policy))
    n <- length(listw$neighbours)
    if (n != length(fx)) stop("objects of different length")
    nb <- listw$neighbours
    crd <- card(nb)
    if (any(crd == 0L) && !zero.policy)
        stop("regions with no neighbours found - zero.policy not supported")
    uni_jcs <- apply(model.matrix(~ fx - 1), 2, function(x) {
        x * lag.listw(listw, x, zero.policy=zero.policy)})
    flat_uni_jcs <- apply(uni_jcs, 1, sum)
    ifx <- as.integer(fx)
    k <- length(levels(fx))
    tifx <- table(ifx)/n
    if (k < 2) stop("must be at least two levels in factor")
    wt <- listw$weights
    almost.1 <- 1 - 64 * .Machine$double.eps


    if (con$pysal_rank) {
        if (con$pysal_sim_obs == "GE") {
# esda/crand.py line 328 modified to handle sampling exceptions
# now matches ties.method="min", na.last="keep"
            local_rank <- function(sims, obs, ties.method, nsim) {
                if (!is.finite(obs)) ret <- NA_real_
                else ret <- sum(sims <= obs, na.rm=TRUE)+1
                ret
            }
        } else if (con$pysal_sim_obs == "GT") {
            local_rank <- function(sims, obs, ties.method, nsim) {
                if (!is.finite(obs)) ret <- NA_real_
                else ret <- sum(sims < obs, na.rm=TRUE)+1
                ret
            }
        } else stop("pysal_sim_obs given", con$pysal_sim_obs, "must be GE or GT")
    } else {
        if (!is.na(con$na.last) && is.logical(con$na.last) && con$na.last) {
            local_rank <- function(sims, obs, ties.method, nsim) {
                as.numeric(rank(c(sims, obs), 
                na.last=TRUE, ties.method=ties.method)[(nsim + 1L)])
            }
        } else if (!is.na(con$na.last) && is.logical(con$na.last) && 
            !con$na.last) {
            local_rank <- function(sims, obs, ties.method, nsim) {
                as.numeric(rank(c(sims, obs), 
                na.last=FALSE, ties.method=ties.method)[(nsim + 1L)])
            }
        } else if (is.na(con$na.last)) {
            local_rank <- function(sims, obs, ties.method, nsim) {
                as.numeric(rank(c(sims, obs), 
                na.last=NA, ties.method=ties.method)[(nsim + 1L)])
            }
        } else if (!is.na(con$na.last) && is.character(con$na.last) && 
            con$na.last == "keep") {
            local_rank <- function(sims, obs, ties.method, nsim) {
                as.numeric(rank(c(sims, obs), 
                na.last="keep", ties.method=ties.method)[(nsim + 1L)])
            }
        } else stop('na.last given', con$na.last, 'must be valid')
    }

    env <- new.env()
    assign("crd", card(nb), envir = env) # cardinality
    assign("nb", nb, envir = env) # nbs
    assign("lww", wt, envir = env) # weights
    assign("nsim", nsim, envir=env) # sims
    assign("fxi", fx, envir = env) # x col
    assign("xi", ifx, envir = env) # x col
    assign("tifx", tifx, envir = env) # ifx props
    assign("n", n, envir=env)
    assign("k", k, envir=env)
    assign("lk", levels(fx), envir=env)
    assign("listw", listw, envir=env)
    assign("no_repeat_in_row", no_repeat_in_row, envir=env)
    assign("almost.1", almost.1, envir=env)
    assign("comp_binary", con$comp_binary, envir=env)
    assign("ties.method", con$rank_ties.method, envir=env)
    assign("unique_ceiling", nsim*con$unique_ceiling, envir=env)
    assign("check_reps", con$check_reps, envir=env)
    assign("xtras", con$xtras, envir=env)
    assign("local_rank", local_rank, envir=env)
    assign("obs_uni", flat_uni_jcs, envir=env)
#    assign("local_jcm_BW", local_jcm_BW, envir=env)
#    assign("local_chi", local_chi, envir=env)
#    assign("local_anscombe", local_anscombe, envir=env)
    varlist = ls(envir = env)

    permLICD_int <- function(i, env) {

        timingsi <- list()
        .ptime_starti <- proc.time()
        crdi <- get("crd", envir = env)[i] # get the ith cardinality
        x <- get("xi", envir = env) # get xs
        fx <- get("fxi", envir = env) # get xs
        x_i <- x[-i] # remove the observed x value for cond. permutations
        xi <- x[i]
        nb_i <- get("nb", envir = env)[[i]] # nbs for ith element
        w_i <- get("lww", envir = env)[[i]] # weights for ith element
        nsim <- get("nsim", envir = env) # no. simulations
        permutation <- nsim > 0
        n_i <- get("n", envir=env) - 1L
        tifx <- get("tifx", envir=env)
        listw <- get("listw", envir=env)
        almost.1 <- get("almost.1", envir=env)
        k <- get("k", envir=env)
        lk <- get("lk", envir=env)
        comp_binary <- get("comp_binary", envir=env)
        ties.method <- get("ties.method", envir=env)
        unique_ceiling <- get("unique_ceiling", envir=env)
        check_reps <- get("check_reps", envir=env)
        xtras <- get("xtras", envir=env)
        local_rank <- get("local_rank", envir=env)
        obs_uni_i <- get("obs_uni", envir=env)[i]
#        local_jcm_BW <- get("local_jcm_BW", envir=env)
#        local_chi <- get("local_chi", envir=env)
#        local_anscombe <- get("local_anscombe", envir=env)

        if (crdi == 0L) return(rep(NA_real_, 36))

        x_nb_i <- c(xi, x[nb_i])
        x_nb_i_xi <- x_nb_i == xi
        w_i_i <- c(0, w_i)
        if (comp_binary) c1_comp_obs_i <- sum(x_nb_i_xi)
        else c1_comp_obs_i <- sum(w_i_i * x_nb_i_xi) + 1
        c2_prop_level_i <- tifx[xi]
        if (comp_binary) c3_crdip1 <- crdi + 1
        else c3_crdip1 <- sum(w_i_i) + 1
        c4_comp_bin_BW_i <- pbinom(c1_comp_obs_i, c3_crdip1,
            c2_prop_level_i, lower.tail=TRUE)
        c5_comp_bin_BW_i <- pbinom(c1_comp_obs_i, c3_crdip1,
            c2_prop_level_i, lower.tail=FALSE)
        c5a_comp_bin_BW_i <- pbinom(c1_comp_obs_i-1, c3_crdip1,
            c2_prop_level_i, lower.tail=FALSE)
        if (xtras) {
            O_BW <- c(c1_comp_obs_i, (c3_crdip1-c1_comp_obs_i))
            E_BW <- c3_crdip1 * c(c2_prop_level_i, 1-c2_prop_level_i)
            c6_comp_chi_BW_i <- local_chi(O_BW, E_BW)
            O_k <- table(factor(x_nb_i, levels=1:k))
            E_k <- c3_crdip1 * tifx
            c7_comp_chi_K_i <- local_chi(O_k, E_k)
            c8_comp_ans_BW_i <- local_anscombe(c1_comp_obs_i, c3_crdip1)
        } else {
            c6_comp_chi_BW_i <- c7_comp_chi_K_i <- c8_comp_ans_BW_i <- NA_real_
        }
        timingsi[["compi"]] <- proc.time() - .ptime_starti
        .ptime_starti <- proc.time()

        sub_i <- sort(c(i, nb_i))
        sub_i <- 1:n %in% sub_i
        sub_xi0 <- x[sub_i] == xi
        i_here <- which(which(sub_i) %in% i)
        listw_subi <- spdep::subset.listw(listw, sub_i)
        only_i_jc <- as.numeric(sub_xi0)[i_here] * spdep::lag.listw(listw_subi,
                as.numeric(sub_xi0), zero.policy=zero.policy)[i_here]
 	sn <- spdep::listw2sn(listw_subi)
	wc <- spdep::spweights.constants(listw_subi, zero.policy=zero.policy, 
		adjust.n=adjust.n)
        wc$S02 <- wc$S0*wc$S0
        wc$nwcn1 <- wc$n*wc$n1
        wc$nwcn2 <- wc$nwcn1*wc$n2
        sub_xi <- factor(as.numeric(!sub_xi0)+1, levels=c(1L, 2L),
            labels=c("1", "2"))
        tab <- table(sub_xi)
        kBW <- length(tab)
        y <- factor(paste(sub_xi[sn[,1]], ":", sub_xi[sn[,2]], sep=""),
            levels=as.vector(outer(1:kBW, 1:kBW, 
            FUN=function(X,Y) paste(X,Y,sep=":"))))
        res <- matrix(tapply(sn[,3], y, sum), ncol=kBW)/2
        res[is.na(res)] <- 0
        if (all(dim(res) == c(1L, 1L))) res <- rbind(cbind(res, 0), c(0, 0))

        local_jcm_obs <- local_jcm_BW(res, tab, wc)
        local_jcm_obs_chi_BW <- local_jcm_obs[[1]]
        local_jcm_obs_count_BB <- local_jcm_obs[[2]][1]
        local_jcm_obs_count_WW <- local_jcm_obs[[2]][2]
        local_jcm_obs_count_BW <- local_jcm_obs[[2]][3]
        vs <- local_jcm_obs[[4]]
        vs_non_pos <- any(vs <=  0)
        Oi_Ei_zero <- all(local_jcm_obs[[2]]-local_jcm_obs[[3]] == 0)
        zs <- (local_jcm_obs[[2]] - local_jcm_obs[[3]]) / sqrt(vs)
        local_jcm_obs_z_BB <- zs[1]
        local_jcm_obs_z_WW <- zs[2]
        local_jcm_obs_z_BW <- zs[3]
        local_jcm_all_BB <- all(sub_xi0)
        timingsi[["configi"]] <- proc.time() - .ptime_starti

        if (permutation) {

            no_repeat_in_row <- get("no_repeat_in_row", envir=env)
        # create matrix of replicates for composition
            if (no_repeat_in_row) {
                samples <- .Call("perm_no_replace", as.integer(nsim),
                  as.integer(n_i), as.integer(crdi), PACKAGE="spdep")
                sx_i <- matrix(x_i[samples], ncol=crdi, nrow=nsim)
            } else {
                sx_i <- matrix(sample(x_i, size = crdi * nsim, replace = TRUE),
                    ncol = crdi, nrow = nsim)
            }

            if (xtras) {
                sx_i_bin <- (sx_i == xi)
                uni_sims <- c(sx_i_bin %*% w_i)
                uni_sims_obs_rank <- local_rank(uni_sims, obs_uni_i,
                    ties.method=ties.method, nsim)
            } else {
               uni_sims_obs_rank <- NA_real_
            }


            sx_i <- cbind(rep(xi, times=nsim), sx_i)
            do_reps <- FALSE
            if (check_reps) {
                u_sx_i <- unique(sx_i)
                len_reps <- nrow(u_sx_i)
                if (len_reps < unique_ceiling) do_reps <- TRUE
                if (do_reps) {
                    reps <- integer(len_reps)
                    for (i in seq(along=reps))
                    reps[i] <- sum(apply(sx_i, 1, function(x)
                        all(x == u_sx_i[i,])))
                }
            } else {
                len_reps <- NA_real_
            }


            if (do_reps) {
                u_c1_comp_sim_i <- apply(u_sx_i, 1,
                    function(y) ifelse(comp_binary, sum(y == xi),
                        sum(w_i_i * (y == xi)) + 1))
                c1_comp_sim_i <- rep(u_c1_comp_sim_i, times=reps)
            } else {
                c1_comp_sim_i <- apply(sx_i, 1,
                    function(y) ifelse(comp_binary, sum(y == xi),
                        sum(w_i_i * (y == xi)) + 1))
            }
            c1_comp_sim_i_rank <- local_rank(c1_comp_sim_i,
                c1_comp_obs_i, ties.method=ties.method, nsim=nsim)

            if (do_reps) {
                u_c4_comp_bin_BW_sim_i <- sapply(u_c1_comp_sim_i, function(y) {
                    pbinom(y, c3_crdip1, c2_prop_level_i, lower.tail=TRUE)
                })
                c4_comp_bin_BW_sim_i <- rep(u_c4_comp_bin_BW_sim_i, times=reps)
            } else {
                c4_comp_bin_BW_sim_i <- sapply(c1_comp_sim_i, function(y) {
                    pbinom(y, c3_crdip1, c2_prop_level_i, lower.tail=TRUE)
                })
            }
            c4_comp_bin_BW_sim_i_rank <- local_rank(c4_comp_bin_BW_sim_i,
                c4_comp_bin_BW_i, ties.method=ties.method, nsim=nsim)

            if (do_reps) {
                u_c5_comp_bin_BW_sim_i <- sapply(u_c1_comp_sim_i, function(y) {
                    pbinom(y, c3_crdip1, c2_prop_level_i, lower.tail=FALSE)
                })
                c5_comp_bin_BW_sim_i <- rep(u_c5_comp_bin_BW_sim_i, times=reps)
            } else {
                c5_comp_bin_BW_sim_i <- sapply(c1_comp_sim_i, function(y) {
                    pbinom(y, c3_crdip1, c2_prop_level_i, lower.tail=FALSE)
                })
            }
            c5_comp_bin_BW_sim_i_rank <- local_rank(c5_comp_bin_BW_sim_i,
                c5_comp_bin_BW_i, ties.method=ties.method, nsim=nsim)

            if (do_reps) {
                u_c5a_comp_bin_BW_sim_i <- sapply(u_c1_comp_sim_i, function(y) {
                    pbinom(y-1, c3_crdip1, c2_prop_level_i, lower.tail=FALSE)
                })
                c5a_comp_bin_BW_sim_i <- rep(u_c5a_comp_bin_BW_sim_i, times=reps)
            } else {
                c5a_comp_bin_BW_sim_i <- sapply(c1_comp_sim_i, function(y) {
                    pbinom(y-1, c3_crdip1, c2_prop_level_i, lower.tail=FALSE)
                })
            }
            c5a_comp_bin_BW_sim_i_rank <- local_rank(c5a_comp_bin_BW_sim_i,
                c5a_comp_bin_BW_i, ties.method=ties.method, nsim=nsim)

            if (xtras) {
                if (do_reps) {
                    u_c6_comp_chi_BW_sim_i <- sapply(u_c1_comp_sim_i, function(y) {
                        local_chi(c(y, c3_crdip1-y), E_BW)
                    })
                    c6_comp_chi_BW_sim_i <- rep(u_c6_comp_chi_BW_sim_i, times=reps)
                } else {
                    c6_comp_chi_BW_sim_i <- sapply(c1_comp_sim_i, function(y) {
                        local_chi(c(y, c3_crdip1-y), E_BW)
                    })
                }
                c6_comp_chi_BW_sim_i_rank <- local_rank(c6_comp_chi_BW_sim_i,
                    c6_comp_chi_BW_i, ties.method=ties.method, nsim=nsim)

                if (do_reps) {
                    u_c7_comp_chi_K_sim_i <- apply(u_sx_i, 1, function(y) {
                        O_k <- table(factor(y, levels=1:k))
                        local_chi(O_k, E_k)
                    })
                    c7_comp_chi_K_sim_i <- rep(u_c7_comp_chi_K_sim_i, times=reps)
                } else {
                    c7_comp_chi_K_sim_i <- apply(sx_i, 1, function(y) {
                        O_k <- table(factor(y, levels=1:k))
                        local_chi(O_k, E_k)
                    })
                }
                c7_comp_chi_K_sim_i_rank <- local_rank(c7_comp_chi_K_sim_i,
                    c7_comp_chi_K_i, ties.method=ties.method, nsim=nsim)

                if (do_reps) {
                    u_c8_comp_ans_BW_sim_i <- sapply(u_c1_comp_sim_i, function(s) {
                        local_anscombe(s, c3_crdip1)
                    })
                    c8_comp_ans_BW_sim_i <- rep(u_c8_comp_ans_BW_sim_i, times=reps)
                } else {
                    c8_comp_ans_BW_sim_i <- sapply(c1_comp_sim_i, function(s) {
                        local_anscombe(s, c3_crdip1)
                    })
                }
                c8_comp_ans_BW_sim_i_rank <- local_rank(c8_comp_ans_BW_sim_i,
                    c8_comp_ans_BW_i, ties.method=ties.method, nsim=nsim)
            } else {
                c6_comp_chi_BW_sim_i_rank <- c7_comp_chi_K_sim_i_rank <- 
                    c8_comp_ans_BW_sim_i_rank <- NA_real_
            }

        # create matrix of replicates for configuration
            x_nb_iBW <- x[nb_i] == xi
            if (no_repeat_in_row) {
                samples <- .Call("perm_no_replace", as.integer(nsim),
                    as.integer(crdi), as.integer(crdi), PACKAGE="spdep")
                sx_i_i <- matrix(x_nb_iBW[samples], ncol=crdi, nrow=nsim)
            } else {
                sx_i_i <- matrix(sample(x_nb_iBW, size = crdi * nsim,
                    replace = TRUE), ncol = crdi, nrow = nsim)
            }

            xi1 <- rep(TRUE, times=nsim)
            if (i_here == 1L) sx_i_i <- cbind(xi1, sx_i_i)
            else if (i_here == (crdi+1)) sx_i_i <- cbind(sx_i_i, xi1)
            else sx_i_i <- cbind(sx_i_i[, 1:(i_here-1)], xi1, sx_i_i[, 
                i_here:crdi])

            do_reps_i <- FALSE
            if (check_reps) {
                u_sx_i_i <- unique(sx_i_i)
                len_reps_i <- nrow(u_sx_i_i)
                if (len_reps_i < unique_ceiling) do_reps_i <- TRUE
                if (do_reps_i) {
                    reps_i <- integer(len_reps_i)
                    len_reps_i <- length(reps_i)
                    for (i in seq(along=reps_i)) reps_i[i] <- sum(apply(
                        sx_i_i, 1, function(x) all(x == u_sx_i_i[i,])))
                }
            } else {
                len_reps_i <- NA_real_
            }


            if (do_reps_i) {
                u_local_jcm_sim <- apply(u_sx_i_i, 1, function(y) {
                    sub_xi <- factor(as.numeric(!y)+1, levels=c(1L, 2L),
                        labels=c("1", "2"))
                    tab <- table(sub_xi)
                    kBW <- length(tab)
                    yy <- factor(paste(sub_xi[sn[,1]], ":", sub_xi[sn[,2]],
                        sep=""), levels=as.vector(outer(1:kBW, 1:kBW, 
                        FUN=function(X,Y) paste(X,Y,sep=":"))))
                    res <- matrix(tapply(sn[,3], yy, sum), ncol=kBW)/2
                    res[is.na(res)] <- 0
                    if (all(dim(res) == c(1L, 1L))) res <- rbind(cbind(res, 0),
                        c(0, 0))
                    local_jcm_BW(res, tab, wc)
                }, simplify=FALSE)
                u_local_jcm_chi_BW_sim <- sapply(u_local_jcm_sim, "[[", 1)
                local_jcm_chi_BW_sim <- rep(u_local_jcm_chi_BW_sim, times=reps_i)
            } else {
                local_jcm_sim <- apply(sx_i_i, 1, function(y) {
                    sub_xi <- factor(as.numeric(!y)+1, levels=c(1L, 2L),
                        labels=c("1", "2"))
                    tab <- table(sub_xi)
                    kBW <- length(tab)
                    yy <- factor(paste(sub_xi[sn[,1]], ":", sub_xi[sn[,2]],
                        sep=""), levels=as.vector(outer(1:kBW, 1:kBW, 
                        FUN=function(X,Y) paste(X,Y,sep=":"))))
                    res <- matrix(tapply(sn[,3], yy, sum), ncol=kBW)/2
                    res[is.na(res)] <- 0
                    if (all(dim(res) == c(1L, 1L))) res <- rbind(cbind(res, 0),
                        c(0, 0))
                    local_jcm_BW(res, tab, wc)
                }, simplify=FALSE)
                local_jcm_chi_BW_sim <- sapply(local_jcm_sim, "[[", 1)
            }
            local_jcm_chi_BW_sim_rank <- local_rank(local_jcm_chi_BW_sim,
                local_jcm_obs_chi_BW, ties.method=ties.method, nsim=nsim)

            if (do_reps_i) {
                u_zs <- t(sapply(u_local_jcm_sim, function(y) {
                    (y[[2]] - y[[3]]) / sqrt(y[[4]])
                }))
                zs1 <- rep(u_zs[,1], times=reps_i)
                zs2 <- rep(u_zs[,2], times=reps_i)
                zs3 <- rep(u_zs[,3], times=reps_i)
            } else {
                zs <- t(sapply(local_jcm_sim, function(y) {
                    (y[[2]] - y[[3]]) / sqrt(y[[4]])
                }))
                zs1 <- zs[,1]
                zs2 <- zs[,2]
                zs3 <- zs[,3]
            }
            local_jcm_z_BB_sim_rank <- local_rank(zs1,
                local_jcm_obs_z_BB, ties.method=ties.method, nsim=nsim)
            local_jcm_z_WW_sim_rank <- local_rank(zs2,
                local_jcm_obs_z_WW, ties.method=ties.method, nsim=nsim)
            local_jcm_z_BW_sim_rank <- local_rank(zs3,
                local_jcm_obs_z_BW, ties.method=ties.method, nsim=nsim)
        } else {
            c1_comp_sim_i_rank <- c4_comp_bin_BW_sim_i_rank <-
            c5_comp_bin_BW_sim_i_rank <- c5a_comp_bin_BW_sim_i_rank <- 
            c6_comp_chi_BW_sim_i_rank <- c7_comp_chi_K_sim_i_rank <- 
            c8_comp_ans_BW_sim_i_rank <- local_jcm_chi_BW_sim_rank <- 
            local_jcm_z_BB_sim_rank <- local_jcm_z_BW_sim_rank <- 
            local_jcm_z_WW_sim_rank <- len_reps <- 
            len_reps_i <- vs_non_pos <- Oi_Ei_zero <-
            uni_sims_obs_rank <- NA_real_
        }
        tcompi <- as.numeric(timingsi[["compi"]][1])
        tconfigi <- as.numeric(timingsi[["configi"]][1])

       
        res_i <- c(xi, c1_comp_obs_i, unname(c2_prop_level_i), c3_crdip1,
            c4_comp_bin_BW_i, c5_comp_bin_BW_i, c5a_comp_bin_BW_i,
            c6_comp_chi_BW_i, c7_comp_chi_K_i, c8_comp_ans_BW_i,
            local_jcm_obs_chi_BW, local_jcm_obs_count_BB,
            local_jcm_obs_count_BW, local_jcm_obs_count_WW, 
            local_jcm_obs_z_BB, local_jcm_obs_z_WW, local_jcm_obs_z_BW,
            c1_comp_sim_i_rank, c4_comp_bin_BW_sim_i_rank,
            c5_comp_bin_BW_sim_i_rank, c5a_comp_bin_BW_sim_i_rank, 
            c6_comp_chi_BW_sim_i_rank, c7_comp_chi_K_sim_i_rank,
            c8_comp_ans_BW_sim_i_rank, local_jcm_chi_BW_sim_rank,
            local_jcm_z_BB_sim_rank, local_jcm_z_BW_sim_rank,
            local_jcm_z_WW_sim_rank, local_jcm_all_BB, as.numeric(len_reps), 
            as.numeric(len_reps_i), tcompi, tconfigi, as.numeric(vs_non_pos),
            as.numeric(Oi_Ei_zero), only_i_jc, uni_sims_obs_rank)
        res_i
    }
    timings[["set_up"]] <- proc.time() - .ptime_start
    .ptime_start <- proc.time()
    inputSGO <- get.SubgraphOption()
    invisible(set.SubgraphOption(FALSE))
    out <- run_perm(fun=permLICD_int, idx=1:n, env=env, iseed=iseed,
        varlist=varlist)
    invisible(set.SubgraphOption(inputSGO))

    colnames(out) <- c("category_i", "count_like_i", "prop_i", "count_nbs_i",
        "bin_like_BW", "bin_unlike_BW", "bin_unlike_BW_alt", "chi_BW_i",
        "chi_K_i", "anscombe_BW", "jcm_chi_obs", "jcm_count_BB_obs", 
        "jcm_count_BW_obs", "jcm_count_WW_obs", "local_jcm_obs_z_BB",
        "local_jcm_obs_z_WW", "local_jcm_obs_z_BW",
        "rank_sim_count_like_i", "rank_sim_bin_like_BW",
        "rank_sim_bin_unlike_BW", "rank_sim_bin_unlike_BW_alt",
        "rank_sim_chi_BW", "rank_sim_chi_K", "rank_sim_anscombe_BW",
        "jcm_chi_sim_rank", "jcm_z_BB_sim_rank", "jcm_z_BW_sim_rank",
        "jcm_z_WW_sim_rank", "local_jcm_all_BB", "len_reps", "len_reps_i",
        "tcompi", "tconfigi", "vs_non_pos", "Oi_Ei_zero", "only_i_jc",
        "uni_sims_obs_rank")

    timings[["processing"]] <- proc.time() - .ptime_start
    .ptime_start <- proc.time()

    local_comp <- data.frame(ID=1:n, category_i=out[,1], count_like_i=out[,2],
        prop_i=out[,3], count_nbs_i=out[,4], pbinom_like_BW=out[,5],
        pbinom_unlike_BW=out[,6], pbinom_unlike_BW_alt=out[,7],
        chi_BW_i=out[,8], chi_K_i=out[,9], anscombe_BW=out[,10])

    pval_jcm_obs_BB <- pnorm(out[,15], lower.tail=FALSE)
    pval_jcm_obs_WW <- pnorm(out[,16], lower.tail=FALSE)
    pval_jcm_obs_BW <- pnorm(out[,17], lower.tail=TRUE)
    sameB <- as.logical(out[,29])
    if (any(sameB)) {
        pval_jcm_obs_BB[sameB] <- 0
        pval_jcm_obs_WW[sameB] <- 1
        pval_jcm_obs_BW[sameB] <- 1
    }
    pval_jcm_obs_BB[!is.finite(pval_jcm_obs_BB)] <- 1
    pval_jcm_obs_WW[!is.finite(pval_jcm_obs_WW)] <- 1
    pval_jcm_obs_BW[!is.finite(pval_jcm_obs_BW)] <- 1
    local_config <- data.frame(ID=1:n, jcm_chi_obs=out[,11],
        jcm_count_BB_obs=out[,12], jcm_count_BW_obs=out[,13],
        jcm_count_WW_obs=out[,14], pval_jcm_obs_BB=pval_jcm_obs_BB,
        pval_jcm_obs_WW=pval_jcm_obs_WW, pval_jcm_obs_BW=pval_jcm_obs_BW,
        only_i_jc=out[,36])
    local_comp_sim <- local_config_sim <- local_uni_sim <- NULL
    if (nsim > 0) {
        pr_bpnsim <- probs_lut("pbinom_like_BW", nsim,
            alternative=con$binomial_punif_alternative)
        local_comp_sim <- data.frame(ID=1:n, rank_sim_count_like_i=out[,18],
            pbinom_like_BW=pr_bpnsim[out[,19]],
            pbinom_unlike_BW=pr_bpnsim[out[,20]],
            pbinom_unlike_BW_alt=pr_bpnsim[out[,21]],
            rank_sim_chi_BW=out[,22], rank_sim_chi_K=out[,23],
            rank_sim_anscombe_BW=out[,24])
        pr_jcmnsim <- probs_lut("jcm_same", nsim,
            alternative=con$jcm_same_punif_alternative)
        pr_jcmnsim1 <- probs_lut("jcm_same", nsim,
            alternative=con$jcm_diff_punif_alternative)
        rnk <- out[,26]
        pval_jcm_obs_BB <- ifelse(is.finite(rnk), pr_jcmnsim[rnk], NA_real_)
        rnk <- out[,27]
        pval_jcm_obs_BW <- ifelse(is.finite(rnk), pr_jcmnsim1[rnk], NA_real_)
        rnk <- out[,28]
        pval_jcm_obs_WW <- ifelse(is.finite(rnk), pr_jcmnsim[rnk], NA_real_)
        if (any(sameB)) {
            pval_jcm_obs_BB[sameB] <- 0
            pval_jcm_obs_WW[sameB] <- 1
            pval_jcm_obs_BW[sameB] <- 1
        }
        pval_jcm_obs_BB[is.nan(pval_jcm_obs_BB)] <- 1
        pval_jcm_obs_WW[is.nan(pval_jcm_obs_WW)] <- 1
        pval_jcm_obs_BW[is.nan(pval_jcm_obs_BW)] <- 1
        local_config_sim <- data.frame(ID=1:n, jcm_chi_sim_rank=out[,25],
        pval_jcm_obs_BB=pval_jcm_obs_BB,
        pval_jcm_obs_BW=pval_jcm_obs_BW,
        pval_jcm_obs_WW=pval_jcm_obs_WW)
        if (con$xtras) {
            stopifnot("uni_sims_obs_rank" %in% colnames(out))
            uni_sim_BB <- out[,37]
            pr_uni_sim <- probs_lut("uni_same", nsim,
                alternative=con$uni_jc_same_punif_alternative)
            pr_uni_sim_BB <- pr_uni_sim[floor(uni_sim_BB)]
            local_uni_sim <- data.frame(1:n, flat_uni_jcs, uni_sim_BB,
                pr_uni_sim_BB, fx)
            names(local_uni_sim) <- c("ID", "obs", "rank",
                attr(pr_uni_sim, "Prname"), "fx")
        }
    }

    timings[["postprocessing"]] <- proc.time() - .ptime_start
    res <- list(local_comp=local_comp, local_config=local_config, local_comp_sim=local_comp_sim, local_config_sim=local_config_sim, local_uni_sim=local_uni_sim)
    attr(res, "timings") <- timings
    attr(res, "out") <- out
    attr(res, "ncpus") <- attr(out, "ncpus")
    attr(res, "nsim") <- nsim
    attr(res, "con") <- con
    attr(res, "uni_jcs") <- uni_jcs
    class(res) <- c("licd", class(res))
    res
}

hotspot.licd <- function(obj, type="both", cutoff=0.05, p.adjust="none", 
    droplevels=TRUE, control=list(), ...) {
    con <- list(binomial_sidak=2, binomial_overlap=TRUE, jcm_sidak=3)
    nmsC <- names(con)
    con[(namc <- names(control))] <- control
    if (length(noNms <- namc[!namc %in% nmsC])) 
        warning("unknown names in control: ", paste(noNms, collapse = ", "))
    comp <- config <- uni <- FALSE
    type <- match.arg(type, c("both", "comp", "config", "uni"))
    if (type == "both") comp <- config <- TRUE
    else if (type == "comp") comp <- TRUE
    else if (type == "config") config <- TRUE
    else uni <- TRUE
    if (uni) {
        if (is.null(obj$local_uni_sim))
            stop("no local joincount uni output in licd object")
    }
    binom_sc <- 1-(1-cutoff)^(1/con$binomial_sidak)
    jcm_sc <- 1-(1-cutoff)^(1/con$jcm_sidak)
    local_comp <- NULL
    local_comp_sim <- NULL
    if (comp) {
        if (!con$binomial_overlap) unlike <- obj$local_comp$pbinom_unlike_BW
        else unlike <- obj$local_comp$pbinom_unlike_BW_alt
        like <- obj$local_comp$pbinom_like_BW
        local_comp <- factor(
            ifelse(p.adjust(unlike, p.adjust) < binom_sc, "Cluster", 
            ifelse(p.adjust(like, p.adjust) < binom_sc, "Outlier", "Dispersed")))
        local_comp_sim <- NULL
        if (attr(obj, "nsim") > 0) {
            if (!con$binomial_overlap)
                unlike <- obj$local_comp_sim$pbinom_unlike_BW
            else unlike <- obj$local_comp_sim$pbinom_unlike_BW_alt
            like <- obj$local_comp_sim$pbinom_like_BW
            local_comp_sim <- factor(
            ifelse(p.adjust(unlike, p.adjust) < binom_sc, "Cluster", 
            ifelse(p.adjust(like, p.adjust) < binom_sc, "Outlier", "Dispersed")))
        }
    }
    local_config <- NULL
    local_config_sim <- NULL
    if (config) {
        BB <- obj$local_config$pval_jcm_obs_BB
        WW <- obj$local_config$pval_jcm_obs_WW
        BW <- obj$local_config$pval_jcm_obs_BW
        BB <- p.adjust(BB, p.adjust)
        WW <- p.adjust(WW, p.adjust)
        BW <- p.adjust(BW, p.adjust)
        JC.pvalue <- cbind(BB, WW, BW)
        crit <- suppressWarnings(apply(JC.pvalue, 1,
            function(y) min(y, na.rm=TRUE) < jcm_sc))
        wh13 <- apply(JC.pvalue, 1, function(y) {
            o <- which.min(y)
            if (length(o) == 0L) o <- 2L
            o
        })
        local_config <- rep("No cluster", nrow(JC.pvalue))
        local_config[crit & wh13 == 3] <- "Dispersed"
        local_config[crit & wh13 == 1] <- "Cluster"
        local_config <- factor(local_config)
        local_config_sim <- NULL
        if (attr(obj, "nsim") > 0) {
            BB <- obj$local_config_sim$pval_jcm_obs_BB
            WW <- obj$local_config_sim$pval_jcm_obs_WW
            BW <- obj$local_config_sim$pval_jcm_obs_BW
            BB <- p.adjust(BB, p.adjust)
            WW <- p.adjust(WW, p.adjust)
            BW <- p.adjust(BW, p.adjust)
            JC.pvalue <- cbind(BB, WW, BW)
            crit <- suppressWarnings(apply(JC.pvalue, 1,
                function(y) min(y, na.rm=TRUE) < jcm_sc))
            wh13 <- apply(JC.pvalue, 1, function(y) {
                o <- which.min(y)
                if (length(o) == 0L) o <- 2L
                o
            })
            local_config_sim <- rep("No cluster", nrow(JC.pvalue))
            local_config_sim[crit & wh13 == 3] <- "Dispersed"
            local_config_sim[crit & wh13 == 1] <- "Cluster"
            local_config_sim <- factor(local_config_sim)
        }
    }
    local_uni_sim <- NULL
    if (uni) {
        pv <- obj$local_uni_sim[,4]
        pv <- p.adjust(pv, p.adjust)
        local_uni_sim <- rep(as.character(NA), length(pv))
        local_uni_sim[pv <= cutoff] <- paste("Cluster:",
            obj$local_uni_sim$fx[pv <= cutoff], ":",
            as.character(cutoff), sep="")
        local_uni_sim <- factor(local_uni_sim)
        
    }
    both <- NULL
    both_sim <- NULL
    both_recode <- NULL
    both_recode_sim <- NULL
    if (type == "both") {
        both <- interaction(local_comp, local_config)
        both_recode <- rep("No cluster", length(local_comp))
        both_recode[both == "Cluster.Cluster"] <- "Cluster"
        both_recode[both == "Cluster.No cluster"] <- "Clump"
        both_recode[both == "Outlier.No cluster"] <- "Outlier"
        both_recode[both == "Dispersed.Dispersed"] <- "Dispersed"
        both_recode[both == "Outlier.Dispersed"] <- "Outlier in dispersion area"
        both_recode <- factor(both_recode)
        if (attr(obj, "nsim") > 0) {
            both_sim <- interaction(local_comp_sim, local_config_sim)
            both_recode_sim <- rep("No cluster", length(local_comp))
            both_recode_sim[both_sim == "Cluster.Cluster"] <- "Cluster"
            both_recode_sim[both_sim == "Cluster.No cluster"] <- "Clump"
            both_recode_sim[both_sim == "Outlier.No cluster"] <- "Outlier"
            both_recode_sim[both_sim == "Dispersed.Dispersed"] <- "Dispersed"
            both_recode_sim[both_sim == "Outlier.Dispersed"] <- "Outlier in dispersion area"
            both_recode_sim <- factor(both_recode_sim)
        }
    }
    list(ID=obj$local_comp$ID, local_comp=local_comp,
        local_comp_sim=local_comp_sim, local_config=local_config,
        local_config_sim=local_config_sim, local_uni_sim=local_uni_sim, 
        both=both, both_sim=both_sim, both_recode=both_recode,
        both_recode_sim=both_recode_sim)
}
