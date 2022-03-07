LOSH.mc <- function(x, listw, a = 2, nsim = 99, zero.policy = NULL, na.action = na.fail, 
                    spChk = NULL, adjust.n = TRUE, p.adjust.method = "none") {
  
  stopifnot(is.vector(x))
  if (!inherits(listw, "listw")) 
    stop(paste(deparse(substitute(listw)), "is not a listw object"))
  if (is.null(zero.policy)) 
    zero.policy <- get("zeroPolicy", envir = .spdepOptions)
  stopifnot(is.logical(zero.policy))
  if (!is.null(attr(listw$neighbours, "self.included")) && 
      attr(listw$neighbours, "self.included")) 
    stop("Self included among neighbours")
  if (is.null(spChk)) 
    spChk <- get.spChkOption()
  if (spChk && !chkIDs(x, listw)) 
    stop("Check of data and weights ID integrity failed")
  if (!is.numeric(x)) 
    stop(paste(deparse(substitute(x)), "is not a numeric vector"))
  if (missing(nsim)) 
    stop("nsim must be given")
  cards <- card(listw$neighbours)
  if (!zero.policy && any(cards == 0))
    stop("regions with no neighbours found")
  if (deparse(substitute(na.action)) == "na.pass") 
    stop("na.pass not permitted")
  x <- na.action(x)
  na.act <- attr(x, "na.action")
  if (!is.null(na.act)) {
    subset <- !(1:length(listw$neighbours) %in% na.act)
    listw <- subset(listw, subset, zero.policy = zero.policy)
  }
  n <- length(listw$neighbours)
  rn <- attr(listw, "region.id")
  if (n != length(x)) 
    stop("objects of different length")
  gamres <- suppressWarnings(nsim > gamma(n + 1))
  if (gamres) 
    stop("nsim too large for this number of observations")
  if (nsim < 1) 
    stop("nsim too small")
  if (adjust.n) 
    n <- n - sum(cards == 0L)
  
  res <- LOSH(x, listw, a, FALSE, zero.policy, na.action, spChk)
  res <- cbind(res, "Pr()" = 1)
  
  losh_boot <- function(data, indices, curr_i, ...) {
    var <- data[indices]
    var[curr_i] <- data[curr_i]
    return(LOSH(x = var, ...)[curr_i,"Hi"])
  }
  cores <- get.coresOption()
  if (is.null(cores)) {
    parallel <- "no"
  }
  else {
    parallel <- ifelse(get.mcOption(), "multicore", "snow")
  }
  ncpus <- ifelse(is.null(cores), 1L, cores)
  cl <- NULL
  if (parallel == "snow") {
    cl <- parallel::makeCluster(get.coresOption())
    parallel::clusterExport(cl, list("LOSH", "lag.listw"), envir = environment())
    on.exit(parallel::stopCluster(cl))
    if (is.null(cl)) {
      parallel <- "no"
      warning("no cluster in ClusterOption, parallel set to no")
    }
  }

  pvals <- numeric(length = length(x))
  for(curr_i in 1:length(x)) {
    boot_obj <- boot(x, statistic = losh_boot, curr_i = curr_i, R = nsim, sim = "permutation", 
                listw = listw, a = a, var_hi = FALSE, zero.policy = zero.policy, 
                na.action = na.action, spChk = spChk, parallel = parallel, ncpus = ncpus, cl = cl)
    
    resampled <- append(boot_obj$t, res[curr_i, "Hi"])
    rankboot <- rank(unlist(resampled))
    xrank <- rankboot[length(resampled)]
    diff <- nsim - xrank
    diff <- ifelse(diff > 0, diff, 0)
    pval <- punif((diff + 1)/(nsim + 1))
    if (!is.finite(pval) || pval < 0 || pval > 1) 
      warning("Out-of-range p-value: reconsider test arguments")
    pvals[curr_i] <- pval
  }
  
  res[,"Pr()"] <- pvals
  res[,"Pr()"] <- p.adjustSP(res[,"Pr()"], listw$neighbours, method = p.adjust.method)
  if (!is.null(rn)) 
    rownames(res) <- rn
  if (!is.null(na.act)) 
    attr(res, "na.action") <- na.act
  class(res) <- c("LOSH", "mc.sim", class(res))
  res
}
