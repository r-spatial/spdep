# small functions to calculate join count
jc_clc_calc <- function(x, z, listw) lag.listw(listw, x * z) * (x * z)
jc_bjc_calc <- function(x, z, listw) (x * (1 - z)) * lag.listw(listw, z * (1-x))


local_joincount_bv <- function(x, z, listw, nsim = 199, alternative = "two.sided") {

  # check that x & z contains only 0 and 1 or TRUE FALSE
  if (!(is.numeric(x) || is.logical(x))) stop("`x` must be an integer or logical.")
  if (!(is.numeric(z) || is.logical(z))) stop("`z` must be an integer or logical.")

  # make sure values are 0 or 1
  if (!all(as.integer(x) %in% c(0L, 1L))) {
    stop("`x` must only contain 0 and 1 or TRUE and FALSE")
  }

  # make sure values are 0 or 1
  if (!all(as.integer(z) %in% c(0L, 1L))) {
    stop("`x` must only contain 0 and 1 or TRUE and FALSE")
  }

  # ensure that x, z, and listw elements are same length
  len_check <- all.equal(
    length(z),
    length(x),
    length(listw[["neighbours"]]),
    length(listw[["weights"]])
    )

  if (!len_check) stop("`x`, `z`, and elements of `listw` must be the same length.")

  # check for binary weights
  if (is.null(attr(listw[["weights"]], "B"))) {
    stop("weights must be binary.")
  }

  # determine if CLC or BJC case
  if (any(x == 1 & z == 1)) {
    case <- "CLC"
  } else if (all((z + x) == 1)) {
    # This ensure that BJC case is met
    case <- "BJC"
  } else if (any(x == 0 & z == 0)) {
    # if at this point neither BJC or CLC have been met
    # then check for presence of 0 == 0
    message(
      "Co-location present for non-observed events.\n  e.g. xi == 0 and zi == 0"
    )
    stop("Neither CLC or BJC cases identified.", call. = FALSE)
  } else {
    stop("Neither CLC or BJC cases identified.", call. = FALSE)
  }

  if (case == "BJC") {
    obs <- jc_bjc_calc(x, z, listw)
    index <- which(x == 1L & obs != 0)
    reps <- replicate(nsim, jc_bjc_calc(x, z, permute_listw(listw)))
  } else if (case == "CLC") {
    # matches Pysal Join_Counts_Local_BV
    obs <- jc_clc_calc(x, z, listw)
    index <- which(obs > 0)
    reps <- replicate(nsim, jc_clc_calc(x, z, permute_listw(listw)))
  }


  # greater
  g <- (rowSums(reps <= obs) + 1) / (nsim + 1)
  # less
  l <- (rowSums(reps >= obs) + 1)/ (nsim + 1)
  # two-sided
  ts <- pmin(l, 1 - l)

  p_value <- switch(
    alternative,
    less = l,
    greater = g,
    two.sided = ts
  )

  # assign correct p-values and ensure missing are correct
  p_vals <- rep(NA_real_, length(x))
  p_vals[index] <- p_value[index]

  data.frame(
    join_count = obs,
    p_sim = p_vals
  )

}
