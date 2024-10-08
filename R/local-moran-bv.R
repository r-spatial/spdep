# A function to calculate local bv moran
# used in replicate internal to local_moran_bv
local_moran_bv_calc <- function(x, y, listw) {
  x * lag.listw(listw, y)
}


localmoran_bv <- function(x, y, listw, nsim = 199, scale = TRUE,
  alternative="two.sided", iseed=1L, no_repeat_in_row=FALSE,
  zero.policy=attr(listw, "zero.policy")) {
  stopifnot(length(x) == length(y))
  if(!inherits(listw, "listw")) stop(paste(deparse(substitute(listw)),
    "is not a listw object"))
  if (is.null(zero.policy))
    zero.policy <- get.ZeroPolicyOption()
  stopifnot(is.logical(zero.policy))
  n <- length(listw$neighbours)
  if (n != length(x)) stop("Different numbers of observations")
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
  
  xx <- mean(x)
  ly <- lag.listw(listw, y, zero.policy=zero.policy)
  lyy <- mean(ly)
  lbs <- c("Low", "High")
  quadr <- interaction(cut(x, c(-Inf, xx, Inf), labels=lbs), 
    cut(ly, c(-Inf, lyy, Inf), labels=lbs), sep="-")
  xmed <- median(x)
  lymed <- median(ly)
  quadr_med <- interaction(cut(x, c(-Inf, xmed, Inf), labels=lbs),
    cut(ly, c(-Inf, lymed, Inf), labels=lbs), sep="-")

  # the variables should be scaled and are by default
  if (scale) {
    x <- as.numeric(scale(x))
    y <- as.numeric(scale(y))
  }

  ly <- lag.listw(listw, y, zero.policy=zero.policy)
  quadr_ps <- interaction(cut(x, c(-Inf, 0, Inf), labels=lbs), 
    cut(ly, c(-Inf, 0, Inf), labels=lbs), sep="-")

  cards <- card(listw$neighbours)
  stopifnot(all(cards > 0L))
# FIXME no zero.policy handling
  if (nsim < 1) stop("nsim too small")

  obs <- local_moran_bv_calc(x, y, listw)

  crd <- card(listw$neighbours)
  lww <- listw$weights
  n <- length(lww)

  env <- new.env()
  assign("x", x, envir=env)
  assign("y", y, envir=env)
  assign("crd", crd, envir=env)
  assign("lww", lww, envir=env)
  assign("nsim", nsim, envir=env)
  assign("obs", obs, envir=env)
  assign("n", n, envir=env)
  assign("no_repeat_in_row", no_repeat_in_row, envir=env)
  varlist <- ls(envir = env)

  permI_bv_int <- function(i, env) {
    res_i <- rep(as.numeric(NA), 4) # initialize output
    crdi <- get("crd", envir=env)[i]
    if (crdi > 0) { # if i has neighbours
      nsim <- get("nsim", envir=env)
      xi <- get("x", envir=env)[i]
      y_i <- get("y", envir=env)[-i]
      n_i <- get("n", envir=env) - 1L
      no_repeat_in_row <- get("no_repeat_in_row", envir=env)
      if (no_repeat_in_row) {
        samples <- .Call("perm_no_replace", as.integer(nsim),
          as.integer(n_i), as.integer(crdi), PACKAGE="spdep")
        sy_i <- matrix(y_i[samples], ncol=crdi, nrow=nsim)
      } else {
        sy_i <- matrix(sample(c(y_i), size=crdi*nsim, replace=TRUE),
        ncol=crdi, nrow=nsim) # permute nsim*#neighbours from y[-i]
      }
      wtsi <- get("lww", envir=env)[[i]]
      res_p <- xi * sy_i %*% wtsi
      # res_p length nsim for obs i conditional draws
      res_i[1] <- mean(res_p)
      res_i[2] <- var(res_p)
      obsi <- get("obs", envir=env)[i]
      res_i[3] <- rank(c(res_p, obsi))[(nsim + 1L)]
      res_i[4] <- as.integer(sum(res_p >= obsi))
    }
    res_i
  }

  out <- run_perm(fun=permI_bv_int, idx=1:n, env=env, iseed=iseed, varlist=varlist)
  ncpus <- attr(out, "ncpus")

  res <- matrix(nrow=n, ncol=7)
  res[,1] <- obs
  res[,2] <- out[,1]
  res[,3] <- out[,2]
  res[,4] <- (res[,1] - res[,2])/sqrt(res[,3])
  if (alternative == "two.sided") 
    res[,5] <- 2 * pnorm(abs(res[,4]), lower.tail=FALSE)
  else if (alternative == "greater") 
    res[,5] <- pnorm(res[,4], lower.tail=FALSE)
  else res[,5] <- pnorm(res[,4])
# look-up table
  probs <- probs_lut(stat="Ibv", nsim=nsim, alternative=alternative)
  res[,6] <- probs[as.integer(out[,3])]
  rnk0 <- as.integer(out[,4])
  drnk0 <- nsim - rnk0
  rnk <- ifelse(drnk0 < rnk0, drnk0, rnk0)
# folded
  res[,7] <- (rnk + 1.0) / (nsim + 1.0)
  Prname <- attr(probs, "Prname")
  Prname_rank <- paste0(Prname, " Sim")
  Prname_sim <- "Pr(folded) Sim"
  colnames(res) <- c("Ibvi", "E.Ibvi", "Var.Ibvi", "Z.Ibvi", Prname,
    Prname_rank, Prname_sim)
  attr(res, "quadr") <- data.frame(mean=quadr, median=quadr_med, 
    pysal=quadr_ps)
  class(res) <- c("localmoran", class(res))
  attr(res, "ncpus") <- ncpus
  res
}

