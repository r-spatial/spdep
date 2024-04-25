moran.plot.seismogram <- function(x, listw, locmoran, cv = 2.58, plain = FALSE, zero.policy = FALSE, xlab = NULL, ylab = NULL, plot = TRUE, return_df = TRUE, spChk = NULL) {
  if (!inherits(listw, "listw")) 
    stop(paste(deparse(substitute(listw)), "is not a listw object"))
  if (!inherits(locmoran, "localmoran")) 
    stop(paste(deparse(substitute(locmoran)), "is not a localmoran object"))
  stopifnot(is.vector(x))
  stopifnot(is.logical(plain))
  stopifnot(is.logical(return_df))
  if (is.null(zero.policy)) 
    zero.policy <- get("zeroPolicy", envir = .spdepOptions)
  stopifnot(is.logical(zero.policy))
  xname <- deparse(substitute(x))
  if (!is.numeric(x)) 
    stop(paste(xname, "is not a numeric vector"))
  if (any(is.na(x))) 
    stop("NA in X")
  n <- length(listw$neighbours)
  if (n != length(x)) 
    stop("objects of different length")
  if (is.null(spChk)) 
    spChk <- get.spChkOption()
  if (spChk && !chkIDs(x, listw)) 
    stop("Check of data and weights ID integrity failed")
  if (is.null(xlab)) 
    xlab <- xname
  if (is.null(ylab)) 
    ylab <- paste("spatially lagged", xname)
  Z <- as.vector(scale(x))
  ZLXi <- lag.listw(listw, Z, zero.policy = zero.policy)
  ZIi <- locmoran[, 4]
  b <- ((cv * sqrt(locmoran[, 3])) + locmoran[, 2]) * var(Z) / (Z - mean(Z))
  b2 <- ((-cv * sqrt(locmoran[, 3])) + locmoran[, 2]) * var(Z) / (Z - mean(Z))
  b[which((Z < 0 & ZLXi > 0) | (Z > 0 & ZLXi < 0))] <- b2[which((Z < 0 & ZLXi > 0) | (Z > 0 & ZLXi < 0))]
  if(plot) {
    if(!plain) {
      x_q1 <- Z[which(Z > 0 & ZLXi > 0)]
      y_q1 <- ZLXi[which(Z > 0 & ZLXi > 0)]
      b_q1 <- b[which(Z > 0 & ZLXi > 0)]
      x_q2 <- Z[which(Z > 0 & ZLXi < 0)]
      y_q2 <- ZLXi[which(Z > 0 & ZLXi < 0)]
      b_q2 <- b[which(Z > 0 & ZLXi < 0)]
      x_q3 <- Z[which(Z < 0 & ZLXi < 0)]
      y_q3 <- ZLXi[which(Z < 0 & ZLXi < 0)]
      b_q3 <- b[which(Z < 0 & ZLXi < 0)]
      x_q4 <- Z[which(Z < 0 & ZLXi > 0)]
      y_q4 <- ZLXi[which(Z < 0 & ZLXi > 0)]
      b_q4 <- b[which(Z < 0 & ZLXi > 0)]
    }
    lw.lm <- lm(ZLXi ~ Z)
    plot(Z, ZLXi, xlab="Z", ylab="WZ", pch = 20, cex = 0.33, col = "lightgrey", xlim = c(min(Z),max(Z)), ylim = c(min(ZLXi, b),max(ZLXi, b)))
    abline(h = 0, lty = "dashed", col = "grey30")
    abline(v = 0, lty = "dashed", col = "grey30")
    abline(lw.lm, lty = "dotted", col = "grey40")
    if(!plain) {
      df_q1 <- data.frame(x_q1, b_q1)
      df_q1 <- df_q1[order(x_q1),]
      df_q2 <- data.frame(x_q2, b_q2)
      df_q2 <- df_q2[order(x_q2),]
      df_q3 <- data.frame(x_q3, b_q3)
      df_q3 <- df_q3[order(x_q3),]
      df_q4 <- data.frame(x_q4, b_q4)
      df_q4 <- df_q4[order(x_q4),]
        
      for(i in 1:(length(df_q1$x_q1) - 1)) {
        lines(c(df_q1$x_q1[i], df_q1$x_q1[i+1]), c(df_q1$b_q1[i], df_q1$b_q1[i+1]), type = "l", col = "firebrick")
      }
        
      for(i in 1:length(df_q2$x_q2)) {
        lines(c(df_q2$x_q2[i], df_q2$x_q2[i+1]), c(df_q2$b_q2[i], df_q2$b_q2[i+1]), type = "l", col = "royalblue")
      }
        
      for(i in 1:length(df_q3$x_q3)) {
        lines(c(df_q3$x_q3[i], df_q3$x_q3[i+1]), c(df_q3$b_q3[i], df_q3$b_q3[i+1]), type = "l", col = "firebrick")
      }
        
      for(i in 1:length(df_q4$x_q4)) {
        lines(c(df_q4$x_q4[i], df_q4$x_q4[i+1]), c(df_q4$b_q4[i], df_q4$b_q4[i+1]), type = "l", col = "royalblue")
      }
    }
  }
  if(return_df) {
    res <- data.frame(z = Z, wz = ZLXi, b = b)
    invisible(res)
  }
}
