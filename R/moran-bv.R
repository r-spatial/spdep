# calculate the bivariate moran
moran_bv_calc <- function(x, yj, wt) {
  numerator <- sum(mapply(function(wij, yj, xi) sum(wij * yj * xi),
                          wt, yj, x))
  denominator <- sum(x^2)

  numerator / denominator
}


# a listw implementation to make permutation easier
moran_bv_impl <- function(x, y, listw) {

  nb <- listw[["neighbours"]]
  wt <- listw[["weights"]]

  yj <- find_xj(y, nb)

  moran_bv_calc(x, yj, wt)

}


moran_bv <- function(x, y, listw, nsim = 99, scale = TRUE) {

  # variables should always be centered and scaled
  if (scale) {
    x <- scale(x)
    y <- scale(y)
  }

  # pull out weights and nbs for observed calc
  nb <- listw[["neighbours"]]
  wt <- listw[["weights"]]

  # calculate observed bivariate moran
  obs <- moran_bv_calc(x, find_xj(y, nb),  wt)

  # create simulated distribution using permuted listw object
  reps <- replicate(
    nsim,
    moran_bv_impl(x, y, permute_listw(listw))
    )

  # calculate p-value from replicates
  p_sim <- (sum(obs <= reps) + 1 )/ (nsim + 1)

  list("Ib" = obs,
       # folded p-sim
       p_sim = pmin(p_sim, 1 - p_sim)
  )
}
