# Copyright 2015 by Roger S. Bivand, Yongwan Chun and Daniel A. Griffith

l_max <- function(lw, zero.policy=TRUE, control=list()) {
    tol <- control$tol
    if (is.null(tol)) tol <- .Machine$double.eps^(1/2)
    stopifnot(is.numeric(tol))
    stopifnot(length(tol) == 1)
    trace <- control$trace
    if (is.null(trace)) trace <- FALSE
    stopifnot(is.logical(trace))
    stopifnot(length(trace) == 1)
# n number of observations
    n <- as.integer(length(lw$neighbours))
    maxiter <- control$maxiter
    if (is.null(maxiter)) maxiter <- 6L*(n-2L)
    stopifnot(is.integer(maxiter))
    stopifnot(length(maxiter) == 1)
    if (lw$style != "B" && trace)
      cat("l_max: weights style is: ", lw$style, "\n")
    if (!is.symmetric.glist(lw$neighbours, lw$weights) && trace)
        cat("l_max: asymmetric weights\n")
# size of neighbour sets
#    ni <- card(lw$neighbours)
# size of generalised neighbour weights
    ni <- sapply(lw$weights, sum)
# initialize variables
    constant <- sum(ni^2)
    nil <- ni/constant
    denom <- sum(nil)
    lamlag <- denom/n
    constant0 <- constant
# start while loop
    keepgoing <- TRUE; k <- 0L
    while (keepgoing) {
        k <- k + 1L
        nik <- lag.listw(lw, nil, zero.policy=zero.policy)
        sumsqnik <- sum(nik^2)
        numer <- sum(nik)
        constant <- sqrt(sumsqnik)
        if (abs(denom) < .Machine$double.eps^2) {
            msg <- "divide by zero"
            keepgoing <- FALSE
            break
        }
        lambda1 <- numer/denom
        if (trace) cat(k, lambda1, numer, denom, constant, "\n")
#        if (abs(lamlag - lambda1) < tol) {
        if (abs(constant0 - constant) < tol) {
            lambda1 <- constant
            msg <- "converged"
            keepgoing <- FALSE
            break
        }
#        lamlag <- lambda1
        denom <- numer/constant
        nil <- nik/constant
        constant0 <- constant
        if (k > maxiter) {
            msg <- "iteration limit exceeded"
            keepgoing <- FALSE
            break
      }
    }
    attr(lambda1, "k") <- k
    attr(lambda1, "msg") <- msg
    attr(lambda1, "constant") <- constant
    attr(lambda1, "e1") <- nik
    lambda1
}


lextrB <- function(lw, zero.policy=TRUE, control=list()) {
# must be binary listw object
  stopifnot(lw$style == "B")
  if (!is.symmetric.glist(lw$neighbours, lw$weights))
    stop("lextrB: asymmetric weights\n")
  tol <- control$tol
  if (is.null(tol)) tol <- .Machine$double.eps^(1/2)
  stopifnot(is.numeric(tol))
  stopifnot(length(tol) == 1)
  control$tol <- tol
  trace <- control$trace
  if (is.null(trace)) trace <- FALSE
  control$trace <- trace
  stopifnot(is.logical(trace))
  stopifnot(length(trace) == 1)
# n number of observations
  lwcard <- card(lw$neighbours)
  n <- as.integer(length(lwcard))
  stopifnot(attr(lw$weights, "mode") == "binary")
  maxiter <- control$maxiter
  if (is.null(maxiter)) maxiter <- 6L*(n-2L)
  stopifnot(is.integer(maxiter))
  stopifnot(length(maxiter) == 1)
  control$maxiter <- maxiter
  useC <- control$useC
  if (is.null(useC)) useC <- TRUE
  stopifnot(is.logical(useC))
  stopifnot(length(useC) == 1)
  control$useC <- useC
  resl1 <- l_max(lw=lw, zero.policy=zero.policy, control=control)
  if (attr(resl1, "msg") != "converged") warning("lextrB: l_max not converged")
  resln_2.1 <- lminC_2.1(lw=lw, y=attr(resl1, "e1")/c(resl1), crd=lwcard,
    zero.policy=zero.policy, control=control)
  if (attr(resln_2.1, "msg") != "converged") warning("lextrB: 2.1 not converged")
  resln_2.2 <- lminC_2.2(lw, resln_2.1, crd=lwcard, zero.policy=zero.policy,
    control=control)
  lambda.n <- lminC_2.3(lw, resln_2.2, attr(resln_2.1, "sse"), crd=lwcard,
    zero.policy=zero.policy, control=control)
  if (attr(lambda.n, "msg") != "converged") warning("lextrB: 2.3 not converged")
  res <- c(lambda_n=c(lambda.n), lambda_1=c(resl1))
  attr(res, "en1") <- cbind(en=attr(lambda.n, "en"),
    e1=attr(resl1, "e1")/c(resl1))
  res
}

lminC_2.3 <- function(lw, y, sse.new, crd, zero.policy=TRUE,
  control=list(
  trace=TRUE,
  tol=.Machine$double.eps^(1/2),
  maxiter=6*(length(lw$neighbours)-2), useC=FALSE)) {
## 2-3. updaing Ei with predicted values Ei_hat until a new
# sse gets smaller (but it needs to be smaller than 1e-15)
#
# must be binary listw object
  stopifnot(lw$style == "B")
  tol <- control$tol
  if (is.null(tol)) tol <- .Machine$double.eps^(1/2)
  trace <- control$trace
  if (is.null(trace)) trace <- TRUE
# n number of observations
  n <- length(lw$neighbours)
  maxiter <- control$maxiter
  if (is.null(maxiter)) maxiter <- 6*(n-2)
  cy <-  lag.listw(lw, y, zero.policy=zero.policy)
  keepgoing4 <- TRUE
  iter <- 0L
  RV.lm.fit <- paste(R.version$major, R.version$minor, sep=".") > "3.0.3"
  if (!RV.lm.fit) .lm.fit <- function() {}
  while(keepgoing4) {
    iter <- iter + 1L
    if (RV.lm.fit) {
      lm.y <- .lm.fit(x=cbind(1,cy), y=y)#lm(y ~ cy)
    } else {
      lm.y <- lm.fit(x=cbind(1,cy), y=y)
    }
    sse.new <- crossprod(lm.y$residuals)#summary(lm.y)$sigma
    beta <- lm.y$coefficients
 
### lw$neighbours y cy beta
### n.switch4 y
    ttol <- tol
    if (control$useC) {
#      crd <- card(lw$neighbours)
      uCres23 <- .Call("lmin23", lw$neighbours, y, cy, crd, beta, ttol,
        PACKAGE="spdep")
      y <- uCres23[[1]]
      n.switch4 <- uCres23[[2]]
    } else {
      n.switch4 <- 0L
      for (i in 1:n) {
        neis <- lw$neighbours[[i]] 
        if (neis[1] > 0L) {
          if (abs(y[i] - (beta[1] + beta[2] * cy[i])) > ttol) {
            tmp <- y[i]
            y[i] <- beta[1] + beta[2] * cy[i] 
            cy[neis] <- cy[neis] - tmp + y[i]
            n.switch4 = n.switch4 + 1L
          }
        }  
      }
    } 
###  
    y <- y - mean(y)
    y <- y/sqrt(sum(y^2))
    cy <- lag.listw(lw, y, zero.policy=zero.policy)
    if (trace) cat("Phase 2.3:", iter, n.switch4, sse.new, "\n")
    if (sse.new <= tol) {
      msg <- "converged"
      break#keepgoing4 <- FALSE
    }  
    if (iter >= maxiter) {
      msg <- "iteration limit exceeded"
      break#keepgoing4 <- FALSE
    } 
    sse.old <- sse.new
  }
   
  lambda.n <- sum(y * cy)/sum(y^2)
  attr(lambda.n, "iter") <- iter
  attr(lambda.n, "msg") <- msg
  attr(lambda.n, "en") <- y
  lambda.n
}


lminC_2.2 <- function(lw, res_2.1, crd, zero.policy=TRUE,
  control=list(
  trace=TRUE,
  tol=.Machine$double.eps^(1/2),
  maxiter=6*(length(lw$neighbours)-2), useC=FALSE)) {
## 2-2. updaing Ei with predicted values Ei_hat based on
# slightly different t1 and t2 
## Ei_hat changes with an update of Ei and lag.Ei at every iteration.
# So predicted values are calculated using coefficients. 

# must be binary listw object
  stopifnot(lw$style == "B")
  tol <- control$tol
  if (is.null(tol)) tol <- .Machine$double.eps^(1/2)
  trace <- control$trace
  if (is.null(trace)) trace <- TRUE
# n number of observations
  n <- length(lw$neighbours)
  maxiter <- control$maxiter
  if (is.null(maxiter)) maxiter <- 6*(n-2)
  beta <- attr(res_2.1, "lm.y")$coefficients
  y <- c(res_2.1)
  cy <- lag.listw(lw, y, zero.policy=zero.policy)
### n.switch3 lw$neighbours y cy beta
### n.switch3 y
  if (control$useC) {
#    crd <- card(lw$neighbours)
    uCres22 <- .Call("lmin22", lw$neighbours, y, cy, crd, beta,
      PACKAGE="spdep")
    y <- uCres22[[1]]
    n.switch3 <- uCres22[[2]]
  } else {
    n.switch3 <- 0L
    for (i in 1:n) {
      neis <- lw$neighbours[[i]] 
      if (neis[1] > 0L) {
        t1 <- abs(y[i] - cy[i]) + sum(abs(y[neis] - cy[neis]))
        t2 <- abs(beta[1] + beta[2] * cy[i] - cy[i]) + sum(abs(y[neis] - 
          (cy[neis] - y[i] + beta[1] + beta[2] * cy[i])))
        if (t1 <= t2) {
          tmp <- y[i]
          y[i] <- beta[1] + beta[2] * cy[i] 
          cy[neis] <- cy[neis] - tmp + y[i]
          n.switch3 <- n.switch3 + 1L
        }
      }
    }
  }
###
  y <- y - mean(y)
  y <- y/sqrt(sum(y^2))
  if (trace) cat("Phase 2.2:", n.switch3, "\n")
  y
}

lminC_2.1 <- function(lw, y, crd, zero.policy=TRUE,
  control=list(
  trace=TRUE,
  tol=.Machine$double.eps^(1/2),
  maxiter=6*(length(lw$neighbours)-2), useC=FALSE)) {
## 2-1. updaing Ei with lag.Ei based on comparison of test1 (t1)
# and test2 (t2)
#
# must be binary listw object
  stopifnot(lw$style == "B")
  tol <- control$tol
  if (is.null(tol)) tol <- .Machine$double.eps^(1/2)
  trace <- control$trace
  if (is.null(trace)) trace <- TRUE
# n number of observations
  n <- length(lw$neighbours)
  maxiter <- control$maxiter
  if (is.null(maxiter)) maxiter <- 6*(n-2)
  y <- y - mean(y)
  y <- y/(sqrt(sum(y^2)))
  cy <- lag.listw(lw, y, zero.policy=zero.policy)
  iter <- 0L
# 998 while
  keepgoing2 <- TRUE
  sse.old <- n
  iter <- 0L
  RV.lm.fit <- paste(R.version$major, R.version$minor, sep=".") > "3.0.3"
  if (!RV.lm.fit) .lm.fit <- function() {}
  while (keepgoing2) {
    iter <- iter + 1L
### lw$neighbours y cy
### n.switch y
    if (control$useC) {
#      crd <- card(lw$neighbours)
      uCres21 <- .Call("lmin21", lw$neighbours, y, cy, crd,
        PACKAGE="spdep")
      y <- uCres21[[1]]
      n.switch <- uCres21[[2]]
    } else {
      n.switch <- 0L
      for (i in 1:n) {
        neis <- lw$neighbours[[i]] 
        if (neis[1] > 0L) {
          t1 <- abs(y[i] - cy[i]) + sum(abs(y[neis] - cy[neis]))
          t2 <- abs(-2 * cy[i]) + sum(abs(y[neis] -
            (cy[neis] - y[i] - cy[i])))
          if (t1 <= t2) {
            tmp <- y[i]
            y[i] <- cy[i] * -1                     # update Ei with lag.Ei
            cy[neis] <- cy[neis] - tmp + y[i]
# update lag.Ei with replacing old Ei with new Ei
            n.switch <- n.switch + 1L
          }
        }
      }
###
    }
    y <- y - mean(y)
    y <- y/sqrt(sum(y^2))
    cy <- lag.listw(lw, y, zero.policy=zero.policy)

    if (RV.lm.fit) {
      lm.y <- .lm.fit(x=cbind(1,cy), y=y)#lm(y ~ cy)
    } else {
      lm.y <- lm.fit(x=cbind(1,cy), y=y)
    }
#    lm.y <- .lm.fit(x=cbind(1,cy), y=y)#lm(y ~ cy)
    sse.new <- crossprod(lm.y$residuals)#summary(lm.y)$sigma

    if (iter > maxiter) {
      msg <- "iteration limit exceeded"
      break #keepgoing2 <- FALSE
    }

    if (sse.new < sse.old) {            # Dan used sum of squares of errors. Maybe need to consider to replace sigma with SSE if this code produces a different outcome 
      if (trace) cat("Phase 2.1:", iter, n.switch, sse.old, sse.new, "\n")
      sse.old <- sse.new
    } else {
      msg <- "converged"
      break #keepgoing2 <- FALSE
    } 
  }
  attr(y, "iter") <- iter
  attr(y, "msg") <- msg
  attr(y, "sse") <- sse.new
  attr(y, "lm.y") <- lm.y
  y
}

