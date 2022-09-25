# Steps for calculating local join count:
#
# 1. Find xj
# 2. Find valid p-value indexes
# 3. Calculate observed jc
# 4. Create replications for p-value
# 5. Calculate p-values from reps
# 6. Put table together
# p-values are only reported for those where xi = 1L and have at least 1 neighbor with one xj == 1L value
local_joincount_uni <- function(fx, chosen, listw,
                                alternative = "two.sided",
                                nsim = 199,
                                iseed = NULL) {

  # check that fx is a factor with 2 levels
  stopifnot(is.factor(fx))
  stopifnot(
    "`fx` must have exactly two levels" = length(levels(fx)) == 2
  )

  # check chosen value
  stopifnot(is.character(chosen))
  stopifnot("`chosen` value is not a level of `fx`" = chosen %in% levels(fx))


  # retrieve weights and neighbors
  nb <- listw[["neighbours"]]
  wt <- listw[["weights"]]

  # check for binary weights
  if (is.null(attr(wt, "B"))) {
    stop("weights must be binary.")
  }

  # cast `fx` to integer
  x <- as.integer(fx == chosen)

  # warn if == 1 more than half of the time
  prop_occurence <- (sum(x == 1) / length(x))
  if (!prop_occurence <= 1/2) {
    warning(
      "Chosen level of `fx` should occur less than half the time. \nOccurs in ",
            round(prop_occurence * 100, 1), "% of observations.")
  }

  # find xj values for identifying where p-vals are to be calculated
  xj <-  lapply(nb, FUN = function(nbs_i) x[nbs_i])
  obs <- x * lag.listw(listw, x)

  # identify which observations should have a p-value
  x_index <- which(x == 1L)
  xj_index <- which(unlist(lapply(xj, function(x) any(x == 1L))) == TRUE)
  index <- intersect(xj_index, x_index)

  obs <- x * lag.listw(listw, x)


  # create new environment to pass into parallel / permBB_int function
  env <- new.env()
  assign("crd", card(nb), envir = env) # cardinality
  assign("lww", wt, envir = env) # weights
  assign("nsim", nsim, envir=env) # weights
  assign("xi", x, envir = env) # x col
  assign("obs", obs, envir = env) # observed values
  varlist = ls(envir = env)

  permBB_int <- function(i, env) {

    crdi <- get("crd", envir = env)[i] # get the ith cardinality
    x <- get("xi", envir = env) # get xs
    x_i <- x[-i] # remove the observed x value for cond. permutations
    w_i <- get("lww", envir = env)[[i]] # weights for ith element
    nsim <- get("nsim", envir = env) # no. simulations
    obs <- get("obs", envir = env) # observed values
    # create matrix of replicates
    sx_i <- matrix(sample(x_i,
                          size = crdi * nsim,
                          replace = TRUE),
                   ncol = crdi,
                   nrow = nsim)

    # calculate join counts for replicates
    res_i <- x[i] * (sx_i %*% w_i)

    # return p-value ranks these will calculate p-value
    # uses look up table approach rather than counting
    # no. observations in each tail
    rank(c(res_i, obs[i]))[(nsim + 1)]

  }

  # create look up table of probabilities
  probs <- probs_lut("BB", nsim, alternative)

  # correct for the two-sided case
  # commenting out because this should be handled in `probs_lut()`
  # if (alternative == "two.sided") probs <- probs / 2
  p_ranks <- run_perm(permBB_int, index, env, iseed, varlist)

  p_res <- rep(NA_real_, length(x))
  p_res[index] <- probs[floor(p_ranks)]

  res <- data.frame(obs, p_res)
  colnames(res) <- c("BB", attr(probs, "Prname"))
  res
}



