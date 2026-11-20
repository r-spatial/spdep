normalized.moran.test <- function (x, listw, zero.policy = attr(listw, "zero.policy"), alternative = "greater",
                                   na.action = na.fail,
                                   spChk = NULL, adjust.n = TRUE)
{
  alternative <- match.arg(alternative, c("greater", "less",
                                          "two.sided"))
  wname <- deparse(substitute(listw))
  if (!inherits(listw, "listw"))
    stop(wname, "is not a listw object")
  xname <- deparse(substitute(x))
  if (!is.numeric(x))
    stop(xname, " is not a numeric vector")
  if (is.null(attr(listw$weights, "W")))
    message(wname, " is not a listw object of style W, then it has been row-standardized")
  if (is.null(zero.policy))
    zero.policy <- get.ZeroPolicyOption()
  stopifnot(is.logical(zero.policy))
  stopifnot(length(zero.policy) == 1L)
  if (is.null(spChk))
    spChk <- get.spChkOption()
  stopifnot(is.logical(spChk))
  stopifnot(length(spChk) == 1L)
  if (spChk && !chkIDs(x, listw))
    stop("Check of data and weights ID integrity failed")
  stopifnot(length(na.action) == 1L)
  NAOK <- deparse(substitute(na.action)) == "na.pass"
  x <- na.action(x)
  na.act <- attr(x, "na.action")
  if (!is.null(na.act)) {
    subset <- !(1:length(listw$neighbours) %in% na.act)
    listw <- subset(listw, subset, zero.policy = zero.policy)
  }
  n <- length(listw$neighbours)
  if (n != length(x))
    stop("objects of different length")
  W <- listw2mat(listw)
  if (is.null(attr(listw$weights, "W")))
    listw <- spdep::mat2listw(W, style="W", zero.policy = zero.policy)
  I <- normalized.moran(x, listw, zero.policy = zero.policy,
                        NAOK = NAOK)
  one <- rep(1,n)
  In <- diag(rep(1,n))
  P <- In - one %*% t(one)/n
  D <- diag(rowSums(W))
  invD <- solve(D)
  B <- P %*% invD %*% W
  C <- t(W) %*% invD %*% P %*% invD %*% W
  B_tilde <- (B + t(B)) / 2

  tr <- function(M) sum(diag(M))
  trP <- tr(P)
  trC <- tr(C)
  trB <- tr(B)
  trB2_tilde <- tr(B_tilde %*% B_tilde)
  trP2 <- tr(P %*% P)
  trC2 <- tr(C %*% C)
  trBP <- tr(B_tilde %*% P)
  trBC <- tr(B_tilde %*% C)
  trPC <- tr(P %*% C)

  EI <- -1/sqrt((n-1)*trC)
  VI <- (2/(trP*trC)) * (trB2_tilde + trB^2/(4*trP^2)*trP2 + trB^2/(4*trC^2)*trC2 -
                           trB/trP*trBP - trB/trC*trBC +
                           trB^2/(2*trP*trC)*trPC)
  if (VI < 0)
    warning("Negative variance,\ndistribution of variable does not meet test assumptions")
  ZI <- (I - EI)/sqrt(VI)
  statistic <- ZI
  names(statistic) <- "Normalized Moran I statistic standard deviate"
  if (alternative == "two.sided")
    PrI <- 2 * pnorm(abs(ZI), lower.tail = FALSE)
  else if (alternative == "greater")
    PrI <- pnorm(ZI, lower.tail = FALSE)
  else PrI <- pnorm(ZI)
  if (!is.finite(PrI) || PrI < 0 || PrI > 1)
    warning("Out-of-range p-value: reconsider test arguments")
  vec <- c(I, EI, VI)
  names(vec) <- c("Normalized Moran I statistic", "Expectation", "Variance")
  method <- "Normalized Moran I test under normality"
  data.name <- paste(xname, "\nweights:", wname, ifelse(is.null(na.act), "",
                                                        paste("\nomitted:", paste(na.act, collapse = ", "))),
                     ifelse(adjust.n && isTRUE(any(sum(card(listw$neighbours) ==
                                                         0L))), "\nn reduced by no-neighbour observations",
                            ""), "", "\n")
  res <- list(statistic = statistic, p.value = PrI, estimate = vec,
              alternative = alternative, method = method, data.name = data.name)
  if (!is.null(na.act))
    attr(res, "na.action") <- na.act
  class(res) <- "htest"
  res
}
