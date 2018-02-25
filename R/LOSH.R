LOSH <- function(x, listw, a = 2, var_hi = TRUE, zero.policy = NULL, na.action = na.fail, spChk = NULL) {
  if (is.null(zero.policy)) 
    zero.policy <- get("zeroPolicy", envir = .spdepOptions)
  stopifnot(is.logical(zero.policy))
  n <- length(listw$neighbours)
  #a <- 2    ## "a" could be any other positive value, but the chi-square-based inference is then no longer possible
  if (n != length(x)) 
    stop("Different numbers of observations")
  NAOK <- deparse(substitute(na.action)) == "na.pass"
  if(var_hi) {
    res <- matrix(nrow = n, ncol = 6)
    colnames(res) <- c("Hi", "E.Hi", "Var.Hi", "Z.Hi", "x_bar_i", "ei")
  } else {
    res <- matrix(nrow = n, ncol = 3)
    colnames(res) <- c("Hi", "x_bar_i", "ei")
  }
  ### calculation of the row sums of the spatial weights
  Wi <- vapply(listw$weights, sum, FUN.VALUE = 0.0)
  ### calculation of x_bar_i
  res[,"x_bar_i"] <- lag.listw(listw, x, zero.policy = zero.policy, NAOK = NAOK) / Wi
  ### calculation of ej
  res[,"ei"] <- abs(x - res[,"x_bar_i"])^a
  ### calculation of 1/(h1 * Wi)
  denom_hi <- mean(res[,"ei"]) * Wi
  ### calculation of Hi
  res[,"Hi"] <- lag.listw(listw, res[,"ei"], zero.policy = zero.policy, NAOK = NAOK) / denom_hi
  if(var_hi) {
    ### calculation of the variance of the ei terms (a global value)
    var_ei <- (sum(res[,"ei"]^2) / n) - mean(res[,"ei"])^2
    ### calculation of the variance of Hi
    res[,"Var.Hi"] <- (n-1)^(-1) * denom_hi^(-2) * var_ei * (n * vapply(listw$weights, function(y) sum(y^2), FUN.VALUE = 0.0) - Wi^2)
    ### calculation of Zi
    res[,"Z.Hi"] <- (2 * res[,"Hi"]) / res[,"Var.Hi"]
    ### storing the expectations
    res[,"E.Hi"] <- 1
  }
  res
}
