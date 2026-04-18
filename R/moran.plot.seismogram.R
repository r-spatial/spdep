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

  WX <- lag.listw(listw, scale(x, scale = F), zero.policy = zero.policy)
  if (anyNA(WX)) warning("no-neighbour observation(s) found - use zero.policy=TRUE")
  if(usePadj)
    cv <- min(abs(locmoran[(which(adjusted_p <= alpha/2)),4]))
  else
    cv <- abs(qnorm(1 - alpha, lower.tail = FALSE))
  b <- ((cv * sqrt(locmoran[, 3])) + locmoran[, 2]) * var(x) / (x - mean(x))
  b2 <- ((-cv * sqrt(locmoran[, 3])) + locmoran[, 2]) * var(x) / (x - mean(x))
  b[which((x < mean(x) & WX > mean(WX)) | (x > mean(x) & WX < mean(WX)))] <- b2[which((x < mean(x) & WX > mean(WX)) | (x > mean(x) & WX < mean(WX)))]
  
  x_q1 <- x[which(x > mean(x) & WX > mean(WX))]
  y_q1 <- WX[which(x > mean(x) & WX > mean(WX))]
  b_q1 <- b[which(x > mean(x) & WX > mean(WX))]
  x_q2 <- x[which(x > mean(x) & WX < mean(WX))]
  y_q2 <- WX[which(x > mean(x) & WX < mean(WX))]
  b_q2 <- b[which(x > mean(x) & WX < mean(WX))]
  x_q3 <- x[which(x < mean(x) & WX < mean(WX))]
  y_q3 <- WX[which(x < mean(x) & WX < mean(WX))]
  b_q3 <- b[which(x < mean(x) & WX < mean(WX))]
  x_q4 <- x[which(x < mean(x) & WX > mean(WX))]
  y_q4 <- WX[which(x < mean(x) & WX > mean(WX))]
  b_q4 <- b[which(x < mean(x) & WX > mean(WX))]

  lw.lm <- lm(WX ~ x)
  plot(x, WX, xlab="X", ylab="WX", pch = 20, cex = 0.5, col = "gray70", xlim = c(min(x),max(x)), ylim = c(min(WX, b),max(WX, b)))
  abline(h = mean(WX), lty = "dashed", col = "grey30")
  abline(v = mean(x), lty = "dashed", col = "grey30")
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
    res <- data.frame(x = x, WX = WX, b = b)
    invisible(res)
  }
}
