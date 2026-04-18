moran.plot.drop <- function(x, listw, locmoran, alpha = 0.05, adjusted_p = NULL, significant = TRUE, xlab = NULL, ylab = NULL, return_df = TRUE, spChk = NULL, labels = NULL, zero.policy=attr(listw, "zero.policy")) {
  if (!inherits(listw, "listw")) 
    stop(paste(deparse(substitute(listw)), "is not a listw object"))
  if (!inherits(locmoran, "localmoran")) 
    stop(paste(deparse(substitute(locmoran)), "is not a localmoran object"))
  stopifnot(is.vector(x))
  stopifnot(is.logical(significant))
  stopifnot(is.logical(return_df))
  stopifnot(is.numeric(alpha))
  if (is.null(zero.policy))
    zero.policy <- get.ZeroPolicyOption()
  stopifnot(is.logical(zero.policy))
  xname <- deparse(substitute(x))
  if (!is.numeric(x)) 
    stop(paste(xname, "is not a numeric vector"))
  if (anyNA(x)) 
    stop("NA in X")
  n <- length(listw$neighbours)
  if (n != length(x)) 
    stop("objects of different length")
  usePadj <- FALSE
  if (!is.null(adjusted_p)) {
    if(n != length(adjusted_p))
      stop("objects of different length")
    if (anyNA(adjusted_p)) 
      stop("NA in vector of adjusted p values")
    usePadj <- TRUE
  }
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
  
  WX <- lag.listw(listw, scale(x, scale = F), zero.policy = zero.policy)
  if (anyNA(WX)) warning("no-neighbour observation(s) found - use zero.policy=TRUE")
  if(usePadj)
    cv <- min(abs(locmoran[(which(adjusted_p <= alpha/2)), 4]))
  else
    cv <- abs(qnorm(1 - alpha, lower.tail = FALSE))

  b <- ((cv * sqrt(locmoran[, 3])) + locmoran[, 2]) * var(x) / (x - mean(x))
  b2 <- ((-cv * sqrt(locmoran[, 3])) + locmoran[, 2]) * var(x) / (x - mean(x))
  b[which((x < mean(x) & WX > mean(WX)) | (x > mean(x) & WX < mean(WX)))] <- b2[which((x < mean(x) & WX > mean(WX)) | (x > mean(x) & WX < mean(WX)))]
  
  if(significant) {
    x_q1 <- x[which(x > mean(x) & WX > mean(WX) & locmoran[, 4] >= cv)]
    y_q1 <- WX[which(x > mean(x) & WX > mean(WX) & locmoran[, 4] >= cv)]
    b_q1 <- b[which(x > mean(x) & WX > mean(WX) & locmoran[, 4] >= cv)]
    labels_q1 <- labels[which(x > mean(x) & WX > mean(WX) & locmoran[, 4] >= cv)]
    x_q2 <- x[which(x > mean(x) & WX < mean(WX) & locmoran[, 4] <= (-1) * cv)]
    y_q2 <- WX[which(x > mean(x) & WX < mean(WX) & locmoran[, 4] <= (-1) * cv)]
    b_q2 <- b[which(x > mean(x) & WX < mean(WX) & locmoran[, 4] <= (-1) * cv)]
    labels_q2 <- labels[which(x > mean(x) & WX < mean(WX) & locmoran[, 4] <= (-1) * cv)]
    x_q3 <- x[which(x < mean(x) & WX < mean(WX) & locmoran[, 4] >= cv)]
    y_q3 <- WX[which(x < mean(x) & WX < mean(WX) & locmoran[, 4] >= cv)]
    b_q3 <- b[which(x < mean(x) & WX < mean(WX) & locmoran[, 4] >= cv)]
    labels_q3 <- labels[which(x < mean(x) & WX < mean(WX) & locmoran[, 4] >= cv)]
    x_q4 <- x[which(x < mean(x) & WX > mean(WX) & locmoran[, 4] <= (-1) * cv)]
    y_q4 <- WX[which(x < mean(x) & WX > mean(WX) & locmoran[, 4] <= (-1) * cv)]
    b_q4 <- b[which(x < mean(x) & WX > mean(WX) & locmoran[, 4] <= (-1) * cv)]
    labels_q4 <- labels[which(x < mean(x) & WX > mean(WX) & locmoran[, 4] <= (-1) * cv)]
  } else {
    x_q1 <- x[which(x > mean(x) & WX > mean(WX) & locmoran[, 4] < cv)]
    y_q1 <- WX[which(x > mean(x) & WX > mean(WX) & locmoran[, 4] < cv)]
    b_q1 <- b[which(x > mean(x) & WX > mean(WX) & locmoran[, 4] < cv)]
    x_q2 <- x[which(x > mean(x) & WX < mean(WX) & locmoran[, 4] > (-1) * cv)]
    y_q2 <- WX[which(x > mean(x) & WX < mean(WX) & locmoran[, 4] > (-1) * cv)]
    b_q2 <- b[which(x > mean(x) & WX < mean(WX) & locmoran[, 4] > (-1) * cv)]
    x_q3 <- x[which(x < mean(x) & WX < mean(WX) & locmoran[, 4] < cv)]
    y_q3 <- WX[which(x < mean(x) & WX < mean(WX) & locmoran[, 4] < cv)]
    b_q3 <- b[which(x < mean(x) & WX < mean(WX) & locmoran[, 4] < cv)]
    x_q4 <- x[which(x < mean(x) & WX > mean(WX) & locmoran[, 4] > (-1) * cv)]
    y_q4 <- WX[which(x < mean(x) & WX > mean(WX) & locmoran[, 4] > (-1) * cv)]
    b_q4 <- b[which(x < mean(x) & WX > mean(WX) & locmoran[, 4] > (-1) * cv)]
  }
  
  lw.lm <- lm(WX ~ x)
  plot(x, WX, xlab="X", ylab="WX", pch = 20, cex = 0.5, col = "gray70", xlim = c(min(x),max(x)), ylim = c(min(WX, b),max(WX, b)))
  abline(h = mean(WX), lty = "dashed", col = "grey30")
  abline(v = mean(x), lty = "dashed", col = "grey30")
  abline(lw.lm, lty = "dotted", col = "grey40")
  
  for(i in 1:length(x_q1)) {
    lines(c(x_q1[i], x_q1[i]), c(y_q1[i], b_q1[i]), type = "l", col = "lightpink", lwd = 1)
    points(c(x_q1[i]), c(y_q1[i]), pch = 20, cex = 0.5, col = "firebrick")
    if (labs) 
      text(x_q1[i], y_q1[i], labels = labels_q1[i], pos = 2, cex = 0.5, col = "firebrick")
  }

  for(i in 1:length(x_q2)) {
    lines(c(x_q2[i], x_q2[i]), c(y_q2[i], b_q2[i]), type = "l", col = "lightblue1", lwd = 1)
    points(c(x_q2[i]), c(y_q2[i]), pch = 20, cex = 0.5, col = "royalblue")
    if (labs) 
      text(x_q2[i], y_q2[i], labels = labels_q2[i], pos = 2, cex = 0.5, col = "royalblue")
  }
      
  for(i in 1:length(x_q3)) {
    lines(c(x_q3[i], x_q3[i]), c(y_q3[i], b_q3[i]), type = "l", col = "lightpink", lwd = 1)
    points(c(x_q3[i]), c(y_q3[i]), pch = 20, cex = 0.5, col = "firebrick")
    if (labs) 
      text(x_q3[i], y_q3[i], labels = labels_q3[i], pos = 2, cex = 0.5, col = "firebrick")
  }
    
  for(i in 1:length(x_q4)) {
    lines(c(x_q4[i], x_q4[i]), c(y_q4[i], b_q4[i]), type = "l", col = "lightblue1", lwd = 1)
    points(c(x_q4[i]), c(y_q4[i]), pch = 20, cex = 0.5, col = "royalblue")
    if (labs) 
      text(x_q4[i], y_q4[i], labels = labels_q4[i], pos = 2, cex = 0.5, col = "royalblue")
  }
  
  if(return_df) {
    res <- data.frame(region.id = labels, x = x, WX = WX, b = b, line_lengths = abs(WX - b))
    invisible(res)
  }
}
