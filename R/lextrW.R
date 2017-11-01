
lextrW <- function(lw, zero.policy=TRUE, control=list()) {
# must be row-standardized listw object

  stopifnot(lw$style == "W")
  stopifnot(attr(lw$weights, "mode") == "binary")
# n number of observations
  lwcard <- as.integer(attr(lw$weights, "comp")$d) #card(lw$neighbours)
  n <- as.integer(length(lwcard))
  stopifnot(can.be.simmed(lw))
  lw <- similar.listw(lw)

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
  if (attr(resl1, "msg") != "converged") warning("lextrW: l_max not converged")
  lwB <- nb2listw(lw$neighbours, style="B", zero.policy=zero.policy)
  resln_2.1 <- lminC_2.1(lw=lwB, y=attr(resl1, "e1")/c(resl1), crd=lwcard,
    zero.policy=zero.policy, control=control)
  if (attr(resln_2.1, "msg") != "converged") warning("lextrW: 2.1 not converged")
  resln_2.2 <- lminC_2.2(lwB, resln_2.1, crd=lwcard, zero.policy=zero.policy,
    control=control)
  lambda.n <- lminC_2.3(lwB, resln_2.2, attr(resln_2.1, "sse"), crd=lwcard, 
    zero.policy=zero.policy, control=control)
  if (attr(lambda.n, "msg") != "converged") warning("lextrW: 2.3 not converged")

  resln_3 <- lminW_3(lw, ev1=attr(lambda.n, "en"), n.nei=lwcard,
    zero.policy=zero.policy, control=control)
  if (attr(resln_3, "msg") != "converged") warning("lextrW: 3 not converged")

  res <- c(lambda_n=c(resln_3), lambda_1=c(resl1))
  attr(res, "en1") <- cbind(en=attr(resln_3, "en"),
    e1=attr(resl1, "e1")/c(resl1))
  res
}

lminW_3 <- function(lw, ev1, n.nei, zero.policy=TRUE,
  control=list(
  trace=TRUE,
  tol=.Machine$double.eps^(1/2),
  maxiter=6*(length(lw$neighbours)-2), useC=FALSE)) {
  stopifnot(lw$style == "W:sim")
  tol <- control$tol
  if (is.null(tol)) tol <- .Machine$double.eps^(1/2)
  trace <- control$trace
  if (is.null(trace)) trace <- TRUE
# n number of observations
  n <- length(lw$neighbours)
  maxiter <- control$maxiter
  if (is.null(maxiter)) maxiter <- 6*(n-2)

#  n.nei <- card(lw$neighbours)
  n.nei.sq <- sqrt(n.nei)
  nn.nei <- n.nei.sq/sqrt(sum(n.nei.sq^2))
  ortho <- sum(ev1 * nn.nei)
  ev1 <- ev1 - ortho/(n*nn.nei)
  ev1 <- ev1/sqrt(sum(ev1^2))   
  ev1.lag <- lag.listw(lw, ev1, zero.policy=zero.policy)
  eval.min.old <- sum(ev1 * ev1.lag)/sum(ev1^2)
   
  n.regress <- 0L
  keepgoing6 <- TRUE
  RV.lm.fit <- paste(R.version$major, R.version$minor, sep=".") > "3.0.3"
  if (!RV.lm.fit) .lm.fit <- function() {}
  while (keepgoing6) {
    n.regress <- n.regress + 1L
    if (RV.lm.fit) {
      lm.y <- .lm.fit(x=cbind(1,ev1.lag), y=ev1)#lm(y ~ cy)
    } else {
      lm.y <- lm.fit(x=cbind(1,ev1.lag), y=ev1)
    }
#    lm.y <- .lm.fit(x=cbind(1,ev1.lag), y=ev1)
    sse.new <- crossprod(lm.y$residuals)
    beta <- lm.y$coefficients
   
    if (control$useC) {
      uCres3 <- .Call("lmin3", lw$neighbours, ev1, ev1.lag, n.nei, beta, tol,
        PACKAGE="spdep")
      ev1 <- uCres3[[1]]
    } else {
      for (i in 1:n) {
        neis <- lw$neighbours[[i]]
        if (neis[1] > 0L) {
          if (abs(ev1[i] - beta[1] + beta[2] * ev1.lag[i]) >= tol) {
            tmp <- ev1[i]
            ev1[i] <- beta[1] + beta[2] * ev1.lag[i] 
            ev1.lag[neis] <- ev1.lag[neis] - tmp/sqrt(n.nei[i] *
              n.nei[neis]) + ev1[i]/sqrt(n.nei[i] * n.nei[neis])
          }   
        }
      }
    }
    ortho <- sum(ev1 * nn.nei)
    ev1 <- ev1 - ortho/(n*nn.nei)
    ev1 <- ev1/sqrt(sum(ev1^2))
   

    ev1.lag <- lag.listw(lw, ev1, zero.policy=zero.policy)
    eval.min.new <- sum(ev1 * ev1.lag)/sum(ev1^2)
    if (trace) cat(n.regress, sse.new, eval.min.new, "\n")

    if (n.regress > maxiter) {
      keepgoing6 <- FALSE
      msg <- "iteration limit exceeded"
      break
    } else if (abs(eval.min.new-eval.min.old) < tol) {
      keepgoing6 <- FALSE 
      msg <- "converged"
      break
    }
    eval.min.old <- eval.min.new
  }
  res <- eval.min.new
  attr(res, "n.regress") <- n.regress
  attr(res, "msg") <- msg
  attr(res, "en") <- ev1
  res
}


lextrS <- function(lw, zero.policy=TRUE, control=list()) {
# must be variance-stabilized listw object
# (possibly already transformed by similarity)

  stopifnot(lw$style == "S")
  stopifnot(attr(lw$weights, "mode") == "binary")
  comp <- attr(lw$weights, "comp")
# FIXME
  lwcard <- card(lw$neighbours) #as.integer(attr(lw$weights, "comp")$q) 
  n <- as.integer(length(lwcard))
  stopifnot(can.be.simmed(lw))
  if (lw$style != "S:sim") lw <- similar.listw(lw)

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
  if (attr(resl1, "msg") != "converged") warning("lextrS: l_max not converged")
  lwB <- nb2listw(lw$neighbours, style="B", zero.policy=zero.policy)
  resln_2.1 <- lminC_2.1(lw=lwB, y=attr(resl1, "e1")/c(resl1), crd=lwcard,
    zero.policy=zero.policy, control=control)
  if (attr(resln_2.1, "msg") != "converged") warning("lextrS: 2.1 not converged")
  resln_2.2 <- lminC_2.2(lwB, resln_2.1, crd=lwcard, zero.policy=zero.policy,
    control=control)
  lambda.n <- lminC_2.3(lwB, resln_2.2, attr(resln_2.1, "sse"), crd=lwcard, 
    zero.policy=zero.policy, control=control)
  if (attr(lambda.n, "msg") != "converged") warning("lextrS: 2.3 not converged")

  resln_3 <- lminS_3(lw, ev1=attr(lambda.n, "en"), comp=comp, crd=lwcard,
    zero.policy=zero.policy, control=control)
  if (attr(resln_3, "msg") != "converged") warning("lextrS: 3 not converged")

  res <- c(lambda_n=c(resln_3), lambda_1=c(resl1))
  attr(res, "en1") <- cbind(en=attr(resln_3, "en"),
    e1=attr(resl1, "e1")/c(resl1))
  res
}

lminS_3 <- function(lw, ev1, comp, crd, zero.policy=TRUE,
  control=list(
  trace=TRUE,
  tol=.Machine$double.eps^(1/2),
  maxiter=6*(length(lw$neighbours)-2), useC=FALSE)) {
#  stopifnot(lw$style == "S:sim")
  tol <- control$tol
  if (is.null(tol)) tol <- .Machine$double.eps^(1/2)
  trace <- control$trace
  if (is.null(trace)) trace <- TRUE
# n number of observations
  n <- length(lw$neighbours)
  maxiter <- control$maxiter
  if (is.null(maxiter)) maxiter <- 6*(n-2)

#  n.nei <- card(lw$neighbours)
  q <- comp$q
  Q <- comp$Q
  eff.n <- comp$eff.n

#  n.nei <- ((eff.n)/Q)*n.nei
  n.nei <- q
  n.nei.sq <- sqrt(n.nei)
  nn.nei <- n.nei.sq/sqrt(sum(n.nei.sq^2))
  ortho <- sum(ev1 * nn.nei)
  ev1 <- ev1 - ortho/(n*nn.nei)
  ev1 <- ev1/sqrt(sum(ev1^2))   
  ev1.lag <- lag.listw(lw, ev1, zero.policy=zero.policy)
  eval.min.old <- sum(ev1 * ev1.lag)/sum(ev1^2)
   
  n.regress <- 0L
  keepgoing6 <- TRUE
  RV.lm.fit <- paste(R.version$major, R.version$minor, sep=".") > "3.0.3"
  if (!RV.lm.fit) .lm.fit <- function() {}
  while (keepgoing6) {
    n.regress <- n.regress + 1L
    if (RV.lm.fit) {
      lm.y <- .lm.fit(x=cbind(1,ev1.lag), y=ev1)#lm(y ~ cy)
    } else {
      lm.y <- lm.fit(x=cbind(1,ev1.lag), y=ev1)
    }
#    lm.y <- .lm.fit(x=cbind(1,ev1.lag), y=ev1)
    sse.new <- crossprod(lm.y$residuals)
    beta <- lm.y$coefficients
   
# FIXME need double n.nei
    if (control$useC) {
      uCres3 <- .Call("lmin3S", lw$neighbours, ev1, ev1.lag, n.nei, crd,
        beta, tol, PACKAGE="spdep")
      ev1 <- uCres3[[1]]
    } else {
      for (i in 1:n) {
        neis <- lw$neighbours[[i]]
        if (neis[1] > 0L) {
          yhat <- beta[1] + beta[2] * ev1.lag[i]
          if (abs(ev1[i] - yhat) >= tol) {
            tmp <- ev1[i]
            ev1[i] <- yhat 
            ev1.lag[neis] <- ev1.lag[neis] - tmp/sqrt(n.nei[i] *
              n.nei[neis]) + ev1[i]/sqrt(n.nei[i] * n.nei[neis])
          }
        }   
      }
    }
    ortho <- sum(ev1 * nn.nei)
    ev1 <- ev1 - ortho/(n*nn.nei)
    ev1 <- ev1/sqrt(sum(ev1^2))
   

    ev1.lag <- lag.listw(lw, ev1, zero.policy=zero.policy)
    eval.min.new <- sum(ev1 * ev1.lag)/sum(ev1^2)
    if (trace) cat(n.regress, sse.new, eval.min.new, "\n")

    if (n.regress > maxiter) {
      keepgoing6 <- FALSE
      msg <- "iteration limit exceeded"
      break
    } else if (abs(eval.min.new-eval.min.old) < tol) {
      keepgoing6 <- FALSE 
      msg <- "converged"
      break
    }
    eval.min.old <- eval.min.new
  }
  res <- eval.min.new
  attr(res, "n.regress") <- n.regress
  attr(res, "msg") <- msg
  attr(res, "en") <- ev1
  res
}


