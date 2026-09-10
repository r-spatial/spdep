normalized.moran.mc <- function(x, listw, nsim, zero.policy = attr(listw, "zero.policy"),
                                alternative = "greater", na.action = na.fail, spChk = NULL,
                                return_boot = FALSE, adjust.n = TRUE)
{
  alternative <- match.arg(alternative, c("greater", "less",
                                          "two.sided"))
  wname <- deparse(substitute(listw))
  if (!inherits(listw, "listw"))
    stop(wname, "is not a listw object")
  xname <- deparse(substitute(x))
  if (!is.numeric(x))
    stop(xname, "is not a numeric vector")
  if (is.null(attr(listw$weights, "W")))
    message(wname, " is not a listw object of style W, then it has been row-standardized")
  if (is.null(zero.policy))
    zero.policy <- get.ZeroPolicyOption()
  stopifnot(is.logical(zero.policy))
  stopifnot(length(zero.policy) == 1L)
  if (missing(nsim))
    stop("nsim must be given")
  if (is.null(spChk))
    spChk <- get.spChkOption()
  stopifnot(is.logical(spChk))
  stopifnot(length(spChk) == 1L)
  if (spChk && !chkIDs(x, listw))
    stop("Check of data and weights ID integrity failed")
  cards <- card(listw$neighbours)
  if (!zero.policy && any(cards == 0))
    stop("regions with no neighbours found")
  stopifnot(length(na.action) == 1L)
  if (deparse(substitute(na.action)) == "na.pass")
    stop("na.pass not permitted")
  x <- na.action(x)
  na.act <- attr(x, "na.action")
  if (!is.null(na.act)) {
    subset <- !(1:length(listw$neighbours) %in% na.act)
    listw <- subset(listw, subset, zero.policy = zero.policy)
    if (return_boot)
      message("NA observations omitted: ", paste(na.act,
                                                 collapse = ", "))
  }
  n <- length(listw$neighbours)
  if (n != length(x))
    stop("objects of different length")
  gamres <- suppressWarnings(nsim > gamma(n + 1))
  if (gamres)
    stop("nsim too large for this number of observations")
  if (nsim < 1)
    stop("nsim too small")
  if (adjust.n)
    n <- n - sum(cards == 0L)
  if (is.null(attr(listw$weights, "W"))) {
    W <- spdep::listw2mat(listw)
    listw <- spdep::mat2listw(W, style="W", zero.policy = zero.policy)
  }
  if (return_boot) {
    moran_boot <- function(var, i, ...) {
      var <- var[i]
      return(normalized.moran(x = var, ...))
    }
    p_setup <- spdep:::parallel_setup(NULL)
    parallel <- p_setup$parallel
    ncpus <- p_setup$ncpus
    cl <- p_setup$cl
    res <- boot::boot(x, statistic = moran_boot, R = nsim, sim = "permutation",
                      listw = listw, zero.policy = zero.policy,
                      parallel = parallel, ncpus = ncpus, cl = cl)
    return(res)
  }
  res <- numeric(length = nsim + 1)
  for (i in 1:nsim) res[i] <- normalized.moran(sample(x), listw,
                                               zero.policy)
  res[nsim + 1] <- normalized.moran(x, listw, zero.policy)
  rankres <- rank(res)
  xrank <- rankres[length(res)]
  diff <- nsim - xrank
  diff <- ifelse(diff > 0, diff, 0)
  if (alternative == "less")
    pval <- punif((diff + 1)/(nsim + 1), lower.tail = FALSE)
  else if (alternative == "greater")
    pval <- punif((diff + 1)/(nsim + 1))
  else pval <- punif(abs(xrank - (nsim + 1)/2)/(nsim + 1),
                     0, 0.5, lower.tail = FALSE)
  if (!is.finite(pval) || pval < 0 || pval > 1)
    warning("Out-of-range p-value: reconsider test arguments")
  statistic <- res[nsim + 1]
  names(statistic) <- "statistic"
  parameter <- xrank
  names(parameter) <- "observed rank"
  method <- "Monte-Carlo simulation of Normalized Moran I"
  data.name <- paste(xname, "\nweights:", wname, ifelse(is.null(na.act),
                                                        "", paste("\nomitted:", paste(na.act, collapse = ", "))),
                     "\nnumber of simulations + 1:", nsim + 1, "\n")
  lres <- list(statistic = statistic, parameter = parameter,
               p.value = pval, alternative = alternative, method = method,
               data.name = data.name, res = res)
  if (!is.null(na.act))
    attr(lres, "na.action") <- na.act
  class(lres) <- c("htest", "mc.sim")
  lres
}
