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
    if (!zero.policy && any(cards == 0))
        stop("regions with no neighbours found")
    ifx <- as.integer(fx)
    k <- length(levels(fx))
    tifx <- table(ifx)/n
    if (k < 2) stop("must be at least two levels in factor")
    wt <- listw$weights

    env <- new.env()
    assign("crd", card(nb), envir = env) # cardinality
    assign("nb", nb, envir = env) # nbs
    assign("lww", wt, envir = env) # weights
    assign("nsim", nsim, envir=env) # sims
    assign("xi", ifx, envir = env) # x col
    assign("tifx", tifx, envir = env) # ifx props
    assign("n", n, envir=env)
    assign("no_repeat_in_row", no_repeat_in_row, envir=env)
    varlist = ls(envir = env)

    permBB_int <- function(i, env) {

        crdi <- get("crd", envir = env)[i] # get the ith cardinality
        x <- get("xi", envir = env) # get xs
        x_i <- x[-i] # remove the observed x value for cond. permutations
        xi <- x[i]
        nb_i <- get("nb", envir = env)[[i]] # nbs for ith element
        w_i <- get("lww", envir = env)[[i]] # weights for ith element
        nsim <- get("nsim", envir = env) # no. simulations
        n_i <- get("n", envir=env) - 1L
        c1_comp_obs_i <- sum(x[nb_i] == xi) + 1
        c2_prop_level_i <- tifx[xi]
        c3_crdip1 <- crdi + 1
        
        no_repeat_in_row <- get("no_repeat_in_row", envir=env)
        # create matrix of replicates
        if (no_repeat_in_row) {
            samples <- .Call("perm_no_replace", as.integer(nsim),
              as.integer(n_i), as.integer(crdi), PACKAGE="spdep")
            sx_i <- matrix(x_i[samples], ncol=crdi, nrow=nsim)
        } else {
            sx_i <- matrix(sample(x_i, size = crdi * nsim, replace = TRUE),
                ncol = crdi, nrow = nsim)
        }

        

    }

}
