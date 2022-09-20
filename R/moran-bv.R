moran_bv <- function(x, y, listw, nsim = 499, scale = TRUE) {
  stopifnot(length(x) == length(y))
  if(!inherits(listw, "listw")) stop(paste(deparse(substitute(listw)),
    "is not a listw object"))
# FIXME is listw assumed to be row-standardized?
  n <- length(listw$neighbours)
  stopifnot(n == length(x))
  if(!is.numeric(x)) stop(paste(deparse(substitute(x)),
    "is not a numeric vector"))
  if(!is.numeric(y)) stop(paste(deparse(substitute(y)),
    "is not a numeric vector"))
# FIXME no na handling for x or y
  if (missing(nsim)) stop("nsim must be given")
  stopifnot(all(!is.na(x)))
  stopifnot(all(!is.na(y)))

  # variables should always be centered and scaled
  if (scale) {
    x <- scale(x)
    y <- scale(y)
  }
  cards <- card(listw$neighbours)
  stopifnot(all(cards > 0L))
# FIXME no zero.policy handling
  if (nsim < 1) stop("nsim too small")

  bvm <- function(x, y, listw) {
    sum(lag.listw(listw, y) * x) / sum(x ^ 2)
  }
  obs <- bvm(x, y, listw)

  xx <- data.frame(x, y)
  bvm_boot <- function(var, i, ...) {
    return(bvm(x=var[i,1], y=var[i,2], ...)) # as lee.mc
  }
  p_setup <- parallel_setup(NULL)
  parallel <- p_setup$parallel
  ncpus <- p_setup$ncpus
  cl <- p_setup$cl
  res <- boot(xx, statistic=bvm_boot, R=nsim,
    sim="permutation", listw=listw, parallel=parallel, ncpus=ncpus, cl=cl)
  return(res)
}
