moran.plot.seismogram <- function(x, listw, locmoran, alpha = 0.05, adjusted_p = NULL, xlab = NULL, ylab = NULL, return_df = TRUE, spChk = NULL, zero.policy=attr(listw, "zero.policy")) {
  if (!inherits(listw, "listw")) 
    stop(paste(deparse(substitute(listw)), "is not a listw object"))
  if (!inherits(locmoran, "localmoran")) 
    stop(paste(deparse(substitute(locmoran)), "is not a localmoran object"))
  stopifnot(is.vector(x))
  stopifnot(is.logical(return_df))
  stopifnot(is.numeric(alpha))
  if (is.null(zero.policy))
    zero.policy <- get.ZeroPolicyOption()
  stopifnot(is.logical(zero.policy))
  xname <- deparse(substitute(x))
  if (!is.numeric(x)) 
    stop(paste(xname, "is not a numeric vector"))
  if (any(is.na(x))) 
    stop("NA in X")
  n <- length(listw$neighbours)
  if (n != length(x)) 
    stop("objects of different length")
  usePadj <- FALSE
  if (!is.null(adjusted_p)) {
    if(n != length(adjusted_p))
      stop("objects of different length")
    if (any(is.na(adjusted_p))) 
      stop("NA in vector of adjusted p values")
    usePadj <- TRUE
  }
  if (is.null(spChk)) 
    spChk <- get.spChkOption()
  if (spChk && !chkIDs(x, listw)) 
    stop("Check of data and weights ID integrity failed")
  if (is.null(xlab)) 
    xlab <- xname
  if (is.null(ylab)) 
    ylab <- paste("spatially lagged", xname)
  Z <- as.vector(scale(x, scale = F))
  WZ <- lag.listw(listw, Z, zero.policy = zero.policy)
  if (anyNA(WZ)) warning("no-neighbour observation(s) found - use zero.policy=TRUE")
  if(usePadj)
    cv <- min(abs(locmoran[(which(adjusted_p <= alpha/2)),4]))
  else
    cv <- qnorm(1 - alpha, lower.tail = FALSE)
  b <- ((cv * sqrt(locmoran[, 3])) + locmoran[, 2]) * var(Z) / (Z - mean(Z))
  b2 <- ((-cv * sqrt(locmoran[, 3])) + locmoran[, 2]) * var(Z) / (Z - mean(Z))
  b[which((Z < 0 & WZ > mean(WZ)) | (Z > 0 & WZ < mean(WZ)))] <- b2[which((Z < 0 & WZ > mean(WZ)) | (Z > 0 & WZ < mean(WZ)))]
  
  x_q1 <- Z[which(Z > 0 & WZ > mean(WZ))]
  y_q1 <- WZ[which(Z > 0 & WZ > mean(WZ))]
  b_q1 <- b[which(Z > 0 & WZ > mean(WZ))]
  x_q2 <- Z[which(Z > 0 & WZ < mean(WZ))]
  y_q2 <- WZ[which(Z > 0 & WZ < mean(WZ))]
  b_q2 <- b[which(Z > 0 & WZ < mean(WZ))]
  x_q3 <- Z[which(Z < 0 & WZ < mean(WZ))]
  y_q3 <- WZ[which(Z < 0 & WZ < mean(WZ))]
  b_q3 <- b[which(Z < 0 & WZ < mean(WZ))]
  x_q4 <- Z[which(Z < 0 & WZ > mean(WZ))]
  y_q4 <- WZ[which(Z < 0 & WZ > mean(WZ))]
  b_q4 <- b[which(Z < 0 & WZ > mean(WZ))]

  lw.lm <- lm(WZ ~ Z)
  plot(Z, WZ, xlab="X_centred", ylab="WX", pch = 20, cex = 0.33, col = "gray70", xlim = c(min(Z),max(Z)), ylim = c(min(WZ, b),max(WZ, b)))
  abline(h = mean(WZ), lty = "dashed", col = "grey30")
  abline(v = 0, lty = "dashed", col = "grey30")
  abline(lw.lm, lty = "dotted", col = "grey40")
    
  df_q1 <- data.frame(x_q1, b_q1)
  df_q1 <- df_q1[order(x_q1),]
  df_q2 <- data.frame(x_q2, b_q2)
  df_q2 <- df_q2[order(x_q2),]
  df_q3 <- data.frame(x_q3, b_q3)
  df_q3 <- df_q3[order(x_q3),]
  df_q4 <- data.frame(x_q4, b_q4)
  df_q4 <- df_q4[order(x_q4),]

  for(i in 1:(length(df_q1$x_q1) - 1))
    lines(c(df_q1$x_q1[i], df_q1$x_q1[i+1]), c(df_q1$b_q1[i], df_q1$b_q1[i+1]), type = "l", col = "firebrick")
  
  for(i in 1:length(df_q2$x_q2)) 
    lines(c(df_q2$x_q2[i], df_q2$x_q2[i+1]), c(df_q2$b_q2[i], df_q2$b_q2[i+1]), type = "l", col = "royalblue")
        
  for(i in 1:length(df_q3$x_q3))
    lines(c(df_q3$x_q3[i], df_q3$x_q3[i+1]), c(df_q3$b_q3[i], df_q3$b_q3[i+1]), type = "l", col = "firebrick")
    
  for(i in 1:length(df_q4$x_q4))
    lines(c(df_q4$x_q4[i], df_q4$x_q4[i+1]), c(df_q4$b_q4[i], df_q4$b_q4[i+1]), type = "l", col = "royalblue")

  if(return_df) {
    res <- data.frame(z = Z, wz = WZ, b = b)
    invisible(res)
  }
}
