# Steps for calculating local join count:
#
# 1. Find xj
# 2. Find valid p-value indexes
# 3. Calculate observed jc
# 4. Create replications for p-value
# 5. Calculate p-values from reps
# 6. Put table together
# p-values are only reported for those where xi = 1L and have at least 1 neighbor with one xj == 1L value
local_joincount_uni <- function(x, listw, nsim = 199, alternative = "two.sided") {

  # check that x contains only 0 and 1 or TRUE FALSE
  if (!(is.numeric(x) || is.logical(x))) stop("`x` must be an integer or logical.")

  # make sure values are 0 or 1
  if (!all(as.integer(x) %in% c(0L, 1L))) {
    stop("`x` must only contain 0 and 1 or TRUE and FALSE")
  }

  # Check for x input type  and proportion of 1s vs TRUES
  nb <- listw[["neighbours"]]
  wt <- listw[["weights"]]

  # check for binary weights
  if (is.null(attr(wt, "B"))) {
    stop("weights must be binary.")
  }

  # warn if == 1 more than half of the time
  prop_occurence <- (sum(x == 1) / length(x))
  if (!prop_occurence <= 1/2) {
    warning("`x` should occur less than half the time. \nOccurs in ",
            round(prop_occurence * 100, 1), "% of observations.")
  }

  # find xj values for identifying where p-vals are to be calculated
  xj <-  lapply(nb, FUN = function(nbs_i) x[nbs_i])
  obs <- x * lag.listw(listw, x)

  # identify which observations should have a p-value
  x_index <- which(x == 1L)
  xj_index <- which(unlist(lapply(xj, function(x) any(x == 1L))) == TRUE)
  index <- intersect(xj_index, x_index)

  # create replicates
  reps <- replicate(
    nsim,
    (x * lag.listw(permute_listw(listw), x))[index]
  )


  g <- (rowSums(reps <= obs[index]) + 1) / (nsim + 1)
  l <- (rowSums(reps >= obs[index]) + 1)/ (nsim + 1)
  p_value <- switch(
    alternative,
    less = l,
    greater = g,
    two.sided = pmin(l, 1 - l)
  )

  p_values <- rep(NA_real_, length(x))
  p_values[index] <- p_value

  data.frame(
    join_count = obs,
    p_sim = p_values
  )
}

