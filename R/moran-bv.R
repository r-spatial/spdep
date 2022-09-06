moran_bv <- function(x, y, listw, nsim = 99, scale = TRUE) {
  # variables should always be centered and scaled
  if (scale) {
    x <- scale(x)
    y <- scale(y)
  }

  obs <- sum(lag.listw(listw, y) * x) / sum(x ^ 2)

  # create simulated distribution using permuted listw object
  reps <- replicate(
    nsim,
    sum(lag.listw(permute_listw(listw), y) * x) / sum(x ^ 2)
  )

  # calculate p-value from replicates
  p_sim <- (sum(obs <= reps) + 1 ) / (nsim + 1)

  list("Ib" = obs,
       # folded p-sim
       p_sim = pmin(p_sim, 1 - p_sim)
  )
}
