moran.plot.seismogram <- function(x, listw, locmoran=NULL, alpha = 0.05, adjusted_p = NULL, xlab = NULL, ylab = NULL, return_df = TRUE, spChk = NULL, zero.policy=attr(listw, "zero.policy"), na.action=na.fail, conditional=TRUE,
 alternative = "two.sided", mlvar=TRUE, adjust.x=FALSE, quadrant.type="mean") {
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
  stopifnot(is.logical(return_df))
  stopifnot(is.numeric(alpha))
  if (!intlocmoran) if (is.null(zero.policy))
    zero.policy <- get.ZeroPolicyOption()
  if (!intlocmoran) stopifnot(is.logical(zero.policy))
  xname <- deparse(substitute(x))
  if (!intlocmoran) if (!is.numeric(x)) 
    stop(paste(xname, "is not a numeric vector"))
  if (!intlocmoran) if (any(is.na(x))) 
    stop("NA in X")
  n <- length(listw$neighbours)
  if (!intlocmoran) if (n != length(x)) 
    stop("objects of different length")
  usePadj <- FALSE
  if (!is.null(adjusted_p)) {
    if(n != length(adjusted_p))
      stop("objects of different length")
    if (any(is.na(adjusted_p))) 
      stop("NA in vector of adjusted p values")
    usePadj <- TRUE
  }
  if (!intlocmoran) if (is.null(spChk)) 
    spChk <- get.spChkOption()
  if (!intlocmoran) if (spChk && !chkIDs(locmoran, listw)) 
    stop("Check of locmoran and weights ID integrity failed")
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
    cv <- min(abs(locmoran[(which(adjusted_p <= alpha/2)),4]))
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
  
  x_q1 <- x[which(HH)] #HH
  y_q1 <- WX[which(HH)]
  b_q1 <- b[which(HH)]
  x_q2 <- x[which(HL)] #HL
  y_q2 <- WX[which(HL)]
  b_q2 <- b[which(HL)]
  x_q3 <- x[which(LL)] #LL
  y_q3 <- WX[which(LL)]
  b_q3 <- b[which(LL)]
  x_q4 <- x[which(LH)] #LH
  y_q4 <- WX[which(LH)]
  b_q4 <- b[which(LH)]

  lw.lm <- lm(WX ~ x)
  plot(x, WX, xlab=xlab, ylab=ylab, pch = 20, cex = 0.5, col = "gray70", xlim = c(min(x),max(x)), ylim = c(min(WX, b),max(WX, b)))
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
    res <- data.frame(labels=as.character(attr(listw, "region.id")), x=x, wx=WX, b=b)
    attr(res, "plot_objs") <- list(df_q1=df_q1, df_q2=df_q2, df_q3=df_q3, df_q4=df_q4)
    invisible(res)
  }
}
