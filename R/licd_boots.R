licd_multi <- function(fx, listw, zero.policy=attr(listw, "zero.policy"),
    adjust.n=TRUE, alternative = "two.sided", nsim = 199, iseed = NULL,             no_repeat_in_row=FALSE) {
    if(!inherits(listw, "listw")) stop(paste(deparse(substitute(listw)),
        "is not a listw object"))
    if(!is.factor(fx)) stop(paste(deparse(substitute(fx)),
        "is not a factor"))
    if (any(is.na(fx))) stop("NA in factor")
    if (is.null(zero.policy))
        zero.policy <- get("zeroPolicy", envir = .spdepOptions)
    stopifnot(is.logical(zero.policy))
    n <- length(listw$neighbours)
    if (n != length(fx)) stop("objects of different length")
    nb <- listw$neighbours
    crd <- card(nb)
    if (any(crd == 0L))
        stop("regions with no neighbours found - zero.policy not supported")
    ifx <- as.integer(fx)
    k <- length(levels(fx))
    tifx <- table(ifx)/n
    if (k < 2) stop("must be at least two levels in factor")
    wt <- listw$weights
    almost.1 <- 1 - 64 * .Machine$double.eps

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
    varlist = ls(envir = env)

    permLICD_int <- function(i, env) {

        crdi <- get("crd", envir = env)[i] # get the ith cardinality
        x <- get("xi", envir = env) # get xs
        fx <- get("fxi", envir = env) # get xs
        x_i <- x[-i] # remove the observed x value for cond. permutations
        xi <- x[i]
        nb_i <- get("nb", envir = env)[[i]] # nbs for ith element
        w_i <- get("lww", envir = env)[[i]] # weights for ith element
        nsim <- get("nsim", envir = env) # no. simulations
        n_i <- get("n", envir=env) - 1L
        tifx <- get("tifx", envir=env)
        listw <- get("listw", envir=env)
        almost.1 <- get("almost.1", envir=env)
        k <- get("k", envir=env)
        lk <- get("lk", envir=env)

        x_nb_i <- c(xi, x[nb_i])
        x_nb_i_xi <- x_nb_i == xi
        c1_comp_obs_i <- sum(x_nb_i_xi)
        c2_prop_level_i <- tifx[xi]
        c3_crdip1 <- crdi + 1
        c4_comp_bin_BW_i <- sum(dbinom(0:c1_comp_obs_i, c3_crdip1,
            c2_prop_level_i))
        c5_comp_bin_BW_i <- 1 - c4_comp_bin_BW_i
        local_chi <- function(O, E) {
            sum((O - E)^2 / E, na.rm=TRUE)
        }
        local_anscombe <- function(s, n) {
            asin(sqrt((s+(3/8))/(n+(3/4))))
        }
        O_BW <- c(c1_comp_obs_i, (c3_crdip1-c1_comp_obs_i))
        E_BW <- c3_crdip1 * c(c2_prop_level_i, 1-c2_prop_level_i)
        c6_comp_chi_BW_i <- local_chi(O_BW, E_BW)
        O_k <- table(factor(x_nb_i, levels=1:k))
        E_k <- c3_crdip1 * tifx
        c7_comp_chi_K_i <- local_chi(O_k, E_k)
        c8_comp_ans_BW_i <- local_anscombe(c1_comp_obs_i, c3_crdip1)

        sub_i <- sort(c(i, nb_i))
        sub_i <- 1:n %in% sub_i
        sub_xi <- x[sub_i] == xi
        i_here <- which(i %in% which(sub_i))
        listw_subi <- subset.listw(listw, sub_i)
 	sn <- listw2sn(listw_subi)
	wc <- spweights.constants(listw_subi, zero.policy=zero.policy, 
		adjust.n=adjust.n)

        local_jcm_BW <- function(xi, sn, wc) {
            xi <- factor(as.numeric(!xi)+1, levels=c(1L, 2L),
                labels=c("1", "2"))
	    tab <- table(xi)
            kBW <- length(tab)
	    y <- factor(paste(xi[sn[,1]], ":",
                xi[sn[,2]], sep=""),
                levels=as.vector(outer(1:kBW, 1:kBW, 
	        FUN=function(X,Y) paste(X,Y,sep=":"))))
	    res <- matrix(tapply(sn[,3], y, sum), ncol=kBW)/2
            res[is.na(res)] <- 0
            if (all(dim(res) == c(1L, 1L))) res <- rbind(cbind(res, 0), c(0, 0))
            BWlk <- c("B", "W")
            ldiag <- res[2,1] + res[1,2]
	    ntab <- as.numeric(as.vector(tab))
	    Ejc <- (wc$S0*(ntab*(ntab-1))) / (2*wc$n*wc$n1)
	    Exp <- (wc$S0*(ntab[2]*ntab[1])) / (wc$n*wc$n1)
	    JtotExp <- sum(Exp)
            O_jc_BW <- c(diag(res), ldiag)
            E_jc_BW <- c(Ejc, Exp)
            jc_conf_chi_BW <- local_chi(O_jc_BW, E_jc_BW)
            list(jc_conf_chi_BW, O_jc_BW, E_jc_BW)
        }

        local_jcm_obs <- local_jcm_BW(sub_xi, sn, wc)
        local_jcm_obs_chi_BW <- local_jcm_obs[[1]]
        local_jcm_obs_count_BB <- local_jcm_obs[[2]][1]
        local_jcm_obs_count_WW <- local_jcm_obs[[2]][2]
        local_jcm_obs_count_BW <- local_jcm_obs[[2]][3]

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

        sx_i <- cbind(rep(xi, times=nsim), sx_i)
        
        c1_comp_sim_i <- apply(sx_i, 1, function(y) sum(y == xi))
        c1_comp_sim_i_rank <- rank(c(c1_comp_sim_i,
            c1_comp_obs_i))[(nsim + 1L)]
        c4_comp_bin_BW_sim_i <- sapply(c1_comp_sim_i, function(y) {
            sum(dbinom(0:y, c3_crdip1, c2_prop_level_i))
        })
        c4_comp_bin_BW_sim_i_rank <- rank(c(c4_comp_bin_BW_sim_i,
            almost.1 * c4_comp_bin_BW_i))[(nsim + 1L)]
        c5_comp_bin_BW_sim_i <- 1 - c4_comp_bin_BW_sim_i
        c5_comp_bin_BW_sim_i_rank <- rank(c(c5_comp_bin_BW_sim_i,
            almost.1 * c5_comp_bin_BW_i))[(nsim + 1L)]
        c6_comp_chi_BW_sim_i <- sapply(c1_comp_sim_i, function(y) {
            local_chi(c(y, c3_crdip1-y), E_BW)
        })
        c6_comp_chi_BW_sim_i_rank <- rank(c(c6_comp_chi_BW_sim_i,
            almost.1 * c6_comp_chi_BW_i))[(nsim + 1L)]
        c7_comp_chi_K_sim_i <- apply(sx_i, 1, function(y) {
            O_k <- table(factor(y, levels=1:k))
            local_chi(O_k, E_k)
        })
        c7_comp_chi_K_sim_i_rank <- rank(c(c7_comp_chi_K_sim_i,
            almost.1 * c7_comp_chi_K_i))[(nsim + 1L)]
        c8_comp_ans_BW_sim_i <- sapply(c1_comp_sim_i, function(s) {
            local_anscombe(s, c3_crdip1)
        })
        c8_comp_ans_BW_sim_i_rank <- rank(c(c8_comp_ans_BW_sim_i,
            almost.1 * c8_comp_ans_BW_i))[(nsim + 1L)]

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

        xi1 <- rep(1L, times=nsim)
        if (i_here == 1L) sx_i <- cbind(xi1, sx_i)
        else if (i_here == (crdi+1)) sx_i <- cbind(sx_i, xi1)
        else sx_i <- cbind(sx_i[, 1:(i_here-1)], xi1, sx_i[, i_here:crdi])

        local_jcm_sim <- apply(sx_i, 1, function(y) {
            local_jcm_BW(y, sn, wc)
        }, simplify=FALSE)

        local_jcm_chi_BW_sim <- sapply(local_jcm_sim, "[[", 1)
        local_jcm_chi_BW_sim_rank <- rank(c(local_jcm_chi_BW_sim,
            almost.1 * local_jcm_obs_chi_BW))[(nsim + 1L)]
        local_jcm_count_BB_sim <- sapply(local_jcm_sim, function(y) y[[2]][1])
        local_jcm_count_BB_sim_rank <- rank(c(local_jcm_count_BB_sim,
            local_jcm_obs_count_BB))[(nsim + 1L)]
        local_jcm_count_WW_sim <- sapply(local_jcm_sim, function(y) y[[2]][2])
        local_jcm_count_WW_sim_rank <- rank(c(local_jcm_count_WW_sim,
            local_jcm_obs_count_WW))[(nsim + 1L)]
        local_jcm_count_BW_sim <- sapply(local_jcm_sim, function(y) y[[2]][3])
        local_jcm_count_BW_sim_rank <- rank(c(local_jcm_count_BW_sim,
            local_jcm_obs_count_BW))[(nsim + 1L)]
       
        res_i <- c(xi, c1_comp_obs_i, unname(c2_prop_level_i), c3_crdip1,
            c4_comp_bin_BW_i, c5_comp_bin_BW_i,
            c6_comp_chi_BW_i, c7_comp_chi_K_i, c8_comp_ans_BW_i,
            local_jcm_obs_chi_BW, local_jcm_obs_count_BB,
            local_jcm_obs_count_BW, local_jcm_obs_count_WW, 
            c1_comp_sim_i_rank, c4_comp_bin_BW_sim_i_rank,
            c5_comp_bin_BW_sim_i_rank, c6_comp_chi_BW_sim_i_rank,
            c7_comp_chi_K_sim_i_rank, c8_comp_ans_BW_sim_i_rank,
            local_jcm_chi_BW_sim_rank, local_jcm_count_BB_sim_rank,
            local_jcm_count_BW_sim_rank, local_jcm_count_WW_sim_rank)
        res_i
    }
    out <- spdep:::run_perm(fun=permLICD_int, idx=1:n, env=env, iseed=iseed, varlist=varlist)
#    ncpus <- attr(out, "ncpus")
    colnames(out) <- c("category_i", "count_like_i", "prop_i", "count_nbs_i",
        "bin_like_BW", "bin_unlike_BW", "chi_BW_i", "chi_K_i", 
        "anscombe_BW", "jcm_chi_obs", "jcm_count_BB_obs", "jcm_count_BW_obs",
        "jcm_count_WW_obs", "rank_sim_count_like_i", "rank_sim_bin_like_BW",
        "rank_sim_bin_unlike_BW", "rank_sim_chi_BW", "rank_sim_chi_K",
        "rank_sim_anscombe_BW", "jcm_chi_sim_rank", "jcm_count_BB_sim_rank",
        "jcm_count_BW_sim_rank", "jcm_count_WW_sim_rank")
    out
}
