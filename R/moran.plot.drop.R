moran.plot.drop <- function(x, listw, nsim = 999, cv = 2.58, significant = TRUE, plain = FALSE, zero.policy = FALSE, xlab = NULL, ylab = NULL, plot = TRUE, return_df = TRUE, spChk = NULL, labels = NULL) {
  require(spdep)
  if (!inherits(listw, "listw")) 
    stop(paste(deparse(substitute(listw)), "is not a listw object"))
  stopifnot(is.vector(x))
  stopifnot(is.logical(significant))
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
  labs <- TRUE
  if (is.logical(labels) && !labels) 
    labs <- FALSE
  if (is.null(labels) || length(labels) != n) 
    labels <- as.character(attr(listw, "region.id"))
  if (is.null(xlab)) 
    xlab <- xname
  if (is.null(ylab)) 
    ylab <- paste("spatially lagged", xname)
  Z <- as.vector(scale(x))
  ZLXi <- lag.listw(listw, Z, zero.policy = zero.policy)
  locmoran <- localmoran_perm(Z, listw, nsim = nsim)
  ZIi <- locmoran[, 4]
  b <- ((cv * sqrt(locmoran[, 3])) + locmoran[, 2]) * var(Z) / (Z - mean(Z))
  b2 <- ((-cv * sqrt(locmoran[, 3])) + locmoran[, 2]) * var(Z) / (Z - mean(Z))
  b[which((Z < 0 & ZLXi > 0) | (Z > 0 & ZLXi < 0))] <- b2[which((Z < 0 & ZLXi > 0) | (Z > 0 & ZLXi < 0))]
  if(plot) {
    if(!plain) {
      if(significant) {
        x_q1 <- Z[which(Z > 0 & ZLXi > 0 & ZIi >= cv)]
        y_q1 <- ZLXi[which(Z > 0 & ZLXi > 0 & ZIi >= cv)]
        b_q1 <- b[which(Z > 0 & ZLXi > 0 & ZIi >= cv)]
        x_q2 <- Z[which(Z > 0 & ZLXi < 0 & ZIi <= (-1) * cv)]
        y_q2 <- ZLXi[which(Z > 0 & ZLXi < 0 & ZIi <= (-1) * cv)]
        b_q2 <- b[which(Z > 0 & ZLXi < 0 & ZIi <= (-1) * cv)]
        x_q3 <- Z[which(Z < 0 & ZLXi < 0 & ZIi >= cv)]
        y_q3 <- ZLXi[which(Z < 0 & ZLXi < 0 & ZIi >= cv)]
        b_q3 <- b[which(Z < 0 & ZLXi < 0 & ZIi >= cv)]
        x_q4 <- Z[which(Z < 0 & ZLXi > 0 & ZIi <= (-1) * cv)]
        y_q4 <- ZLXi[which(Z < 0 & ZLXi > 0 & ZIi <= (-1) * cv)]
        b_q4 <- b[which(Z < 0 & ZLXi > 0 & ZIi <= (-1) * cv)]
      } else {
        x_q1 <- Z[which(Z > 0 & ZLXi > 0 & ZIi < cv)]
        y_q1 <- ZLXi[which(Z > 0 & ZLXi > 0 & ZIi < cv)]
        b_q1 <- b[which(Z > 0 & ZLXi > 0 & ZIi < cv)]
        x_q2 <- Z[which(Z > 0 & ZLXi < 0 & ZIi > (-1) * cv)]
        y_q2 <- ZLXi[which(Z > 0 & ZLXi < 0 & ZIi > (-1) * cv)]
        b_q2 <- b[which(Z > 0 & ZLXi < 0 & ZIi > (-1) * cv)]
        x_q3 <- Z[which(Z < 0 & ZLXi < 0 & ZIi < cv)]
        y_q3 <- ZLXi[which(Z < 0 & ZLXi < 0 & ZIi < cv)]
        b_q3 <- b[which(Z < 0 & ZLXi < 0 & ZIi < cv)]
        x_q4 <- Z[which(Z < 0 & ZLXi > 0 & ZIi > (-1) * cv)]
        y_q4 <- ZLXi[which(Z < 0 & ZLXi > 0 & ZIi > (-1) * cv)]
        b_q4 <- b[which(Z < 0 & ZLXi > 0 & ZIi > (-1) * cv)]
      }
    }
    lw.lm <- lm(ZLXi ~ Z)
    plot(Z, ZLXi, xlab="Z", ylab="WZ", pch = 20, cex = 0.33, col = "lightgrey", xlim = c(min(Z),max(Z)), ylim = c(min(ZLXi, b),max(ZLXi, b)))
    abline(h = 0, lty = "dashed", col = "grey30")
    abline(v = 0, lty = "dashed", col = "grey30")
    abline(lw.lm, lty = "dotted", col = "grey40")
    if(!plain) {
      for(i in 1:length(x_q1)) {
        lines(c(x_q1[i], x_q1[i]), c(y_q1[i], b_q1[i]), type = "l", col = "lightpink", lwd = 1)
        points(c(x_q1[i]), c(y_q1[i]), pch = 20, cex = 0.5, col = "firebrick")
        if (labs && length(x_q1) > 0) 
          text(x_q1[i], y_q1[i], labels = labels[i], pos = 2, cex = 0.5, col = "firebrick")
      }
      
      for(i in 1:length(x_q2)) {
        lines(c(x_q2[i], x_q2[i]), c(y_q2[i], b_q2[i]), type = "l", col = "lightblue1", lwd = 1)
        points(c(x_q2[i]), c(y_q2[i]), pch = 20, cex = 0.5, col = "royalblue")
        if (labs  && length(x_q2) > 0) 
          text(x_q2[i], y_q2[i], labels = labels[i], pos = 2, cex = 0.5, col = "royalblue")
      }
        
      for(i in 1:length(x_q3)) {
        lines(c(x_q3[i], x_q3[i]), c(y_q3[i], b_q3[i]), type = "l", col = "lightpink", lwd = 1)
        points(c(x_q3[i]), c(y_q3[i]), pch = 20, cex = 0.5, col = "firebrick")
        if (labs  && length(x_q3) > 0) 
          text(x_q3[i], y_q3[i], labels = labels[i], pos = 2, cex = 0.5, col = "firebrick")
      }
      
      for(i in 1:length(x_q4)) {
        lines(c(x_q4[i], x_q4[i]), c(y_q4[i], b_q4[i]), type = "l", col = "lightblue1", lwd = 1)
        points(c(x_q4[i]), c(y_q4[i]), pch = 20, cex = 0.5, col = "royalblue")
        if (labs  && length(x_q4) > 0) 
          text(x_q4[i], y_q4[i], labels = labels[i], pos = 2, cex = 0.5, col = "royalblue")
      }
      if (zero.policy) {
        n0 <- ZLXi == 0
        if (any(n0)) {
          symbols(x[n0], wx[n0], inches = FALSE, circles = rep(diff(range(x))/50, length(which(n0))), bg = "grey", add = TRUE)
        }
      }
    }
  }
  if(return_df) {
    res <- data.frame(z = Z, wz = ZLXi, b = b, line_lengths = abs(ZLXi - b))
    invisible(res)
  }
}