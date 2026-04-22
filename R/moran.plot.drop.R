moran.plot.drop <- function(x, listw, locmoran=NULL, alpha = 0.05, adjusted_p = NULL, significant = TRUE, xlab = NULL, ylab = NULL, return_df = TRUE, spChk = NULL, labels = NULL, zero.policy=attr(listw, "zero.policy"),
 na.action=na.fail, conditional=TRUE, alternative = "two.sided", mlvar=TRUE,
 adjust.x=FALSE, quadrant.type="mean") {
  if (!inherits(listw, "listw")) 
    stop(paste(deparse(substitute(listw)), "is not a listw object"))
  intlocmoran <- FALSE
  if (is.null(locmoran)) {
    intlocmoran <- TRUE
    locmoran <- localmoran(x=x, listw=listw,
      zero.policy=zero.policy, na.action=na.action, conditional=conditional, 
      alternative=alternative, mlvar=mlvar, adjust.x=adjust.x)
  }
  if (!inherits(locmoran, "localmoran")) 
    stop(paste(deparse(substitute(locmoran)), "is not a localmoran object"))
  if (!intlocmoran) stopifnot(is.vector(x))
  stopifnot(is.logical(significant))
  stopifnot(is.logical(return_df))
  stopifnot(is.numeric(alpha))
  if (!intlocmoran) if (is.null(zero.policy))
    zero.policy <- get.ZeroPolicyOption()
    if (!intlocmoran) stopifnot(is.logical(zero.policy))
  xname <- deparse(substitute(x))
    if (!intlocmoran) if (!is.numeric(x)) 
      stop(paste(xname, "is not a numeric vector"))
  if (!intlocmoran) if (anyNA(x)) 
    stop("NA in X")
  n <- length(listw$neighbours)
  if (!intlocmoran) if (n != length(x)) 
    stop("objects of different length")
  usePadj <- FALSE
  if (!is.null(adjusted_p)) {
    if(n != length(adjusted_p))
      stop("objects of different length")
    if (anyNA(adjusted_p)) 
      stop("NA in vector of adjusted p values")
    usePadj <- TRUE
  }
  if (!intlocmoran) if (is.null(spChk)) 
    spChk <- get.spChkOption()
  if (!intlocmoran) if (spChk && !chkIDs(locmoran, listw)) 
    stop("Check of data and weights ID integrity failed")
  labs <- TRUE
  if (is.logical(labels)) {
    if(!labels)
      labs <- FALSE
    labels <- as.character(attr(listw, "region.id"))
  } else if (!is.logical(labels) && !is.null(labels)) {
    if(length(labels) != n) {
      warning("Length of the labels vector does not match number of regions. region.id from the listw object is used instead.")
      labels <- as.character(attr(listw, "region.id"))
    }
  } else if (is.null(labels)) {
    labs <- FALSE
    labels <- as.character(attr(listw, "region.id"))
  }
  if (is.null(xlab))
    xlab <- xname
  if (is.null(ylab)) 
    ylab <- paste("spatially lagged centred", xname)
  
  if (intlocmoran) {
    WX <- attr(locmoran, "xlz")$lz
  } else {
    WX <- lag.listw(listw, scale(x, scale = FALSE), zero.policy = zero.policy)
  }
  if (anyNA(WX)) warning("no-neighbour observation(s) found - use zero.policy=TRUE")
  if(usePadj)
    cv <- min(abs(locmoran[(which(adjusted_p <= alpha/2)), 4]))
  else
    cv <- abs(qnorm(1 - alpha, lower.tail = FALSE))

  if (intlocmoran) {
    quadrant.type <- match.arg(quadrant.type, c("mean", "median", "pysal"))
    if (is.null(attr(locmoran, "quadr"))) stop("object has no quadr attribute")
    quadr <- attr(locmoran, "quadr")[[quadrant.type]]
    HH <- quadr == "High-High"
    HL <- quadr == "High-Low"
    LH <- quadr == "Low-High"
    LL <- quadr == "Low-Low"
  } else {
    HH <- x > mean(x) & WX > mean(WX)
    HL <- x > mean(x) & WX < mean(WX)
    LL <- x < mean(x) & WX < mean(WX)
    LH <- x < mean(x) & WX > mean(WX)
  }

  b <- ((cv * sqrt(locmoran[, 3])) + locmoran[, 2]) * var(x) / (x - mean(x))
  b2 <- ((-cv * sqrt(locmoran[, 3])) + locmoran[, 2]) * var(x) / (x - mean(x))
  b[which((LH) | (HL))] <- b2[which((LH) | (HL))]
  
  if(significant) {
    x_q1 <- x[which(HH & locmoran[, 4] >= cv)] #HH
    y_q1 <- WX[which(HH & locmoran[, 4] >= cv)]
    b_q1 <- b[which(HH & locmoran[, 4] >= cv)]
    labels_q1 <- labels[which(HH & locmoran[, 4] >= cv)]
    x_q2 <- x[which(HL & locmoran[, 4] <= (-1) * cv)] #HL
    y_q2 <- WX[which(HL & locmoran[, 4] <= (-1) * cv)]
    b_q2 <- b[which(HL & locmoran[, 4] <= (-1) * cv)]
    labels_q2 <- labels[which(HL & locmoran[, 4] <= (-1) * cv)]
    x_q3 <- x[which(LL & locmoran[, 4] >= cv)] #LL
    y_q3 <- WX[which(LL & locmoran[, 4] >= cv)]
    b_q3 <- b[which(LL & locmoran[, 4] >= cv)]
    labels_q3 <- labels[which(LL & locmoran[, 4] >= cv)]
    x_q4 <- x[which(LH & locmoran[, 4] <= (-1) * cv)] #LH
    y_q4 <- WX[which(LH & locmoran[, 4] <= (-1) * cv)]
    b_q4 <- b[which(LH & locmoran[, 4] <= (-1) * cv)]
    labels_q4 <- labels[which(LH & locmoran[, 4] <= (-1) * cv)]
  } else {
    x_q1 <- x[which(HH & locmoran[, 4] < cv)] #HH
    y_q1 <- WX[which(HH & locmoran[, 4] < cv)]
    b_q1 <- b[which(HH & locmoran[, 4] < cv)]
    x_q2 <- x[which(HL & locmoran[, 4] > (-1) * cv)] #HL
    y_q2 <- WX[which(HL & locmoran[, 4] > (-1) * cv)]
    b_q2 <- b[which(HL & locmoran[, 4] > (-1) * cv)]
    x_q3 <- x[which(LL & locmoran[, 4] < cv)] #LL
    y_q3 <- WX[which(LL & locmoran[, 4] < cv)]
    b_q3 <- b[which(LL & locmoran[, 4] < cv)]
    x_q4 <- x[which(LH & locmoran[, 4] > (-1) * cv)] #LH
    y_q4 <- WX[which(LH & locmoran[, 4] > (-1) * cv)]
    b_q4 <- b[which(LH & locmoran[, 4] > (-1) * cv)]
  }
  
  lw.lm <- lm(WX ~ x)
  plot(x, WX, xlab=xlab, ylab=ylab, pch = 20, cex = 0.5, col = "gray70", xlim = c(min(x),max(x)), ylim = c(min(WX, b),max(WX, b)))
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
    res <- data.frame(labels = labels, x = x, WX = WX, b = b, line_lengths = abs(WX - b))
    attr(res, "plot_objs") <- list(df_q1=data.frame(x_q1, y_q1, b_q1), df_q2=data.frame(x_q2, y_q2, b_q2), df_q3=data.frame(x_q3, y_q3, b_q3), df_q4=data.frame(x_q4, y_q4, b_q4))
    invisible(res)
  }
}
