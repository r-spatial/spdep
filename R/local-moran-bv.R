# A function to calculate local bv moran
# used in replicate internal to local_moran_bv
local_moran_bv_calc <- function(x, y, listw) {
  x * lag.listw(listw, y)
}


localmoran_bv <- function(x, y, listw, nsim = 199, scale = TRUE) {

  # the variables should be scaled and are by default
  if (scale) {
    x <- as.numeric(scale(x))
    y <- as.numeric(scale(y))
  }

  obs <- local_moran_bv_calc(x, y, listw)
  # create replicates as reference distribution
  reps <- replicate(
    nsim,
    local_moran_bv_calc(x, y, permute_listw(listw))
    )

  # calculate folded p-value
  p_sim <- (rowSums(obs <= reps) + 1 )/ (nsim + 1)

  # return results as data.frame
  data.frame("Ib" = obs,
             p_sim = pmin(p_sim, 1 - p_sim)
  )
}

