LOSH.mc <- function(x, listw, a = 2, nsim = 99, zero.policy = NULL, na.action = na.fail, 
                    spChk = NULL, adjust.n = TRUE, p.adjust.method = "none") {
  
  require(foreach) || install.packages("foreach")
  require(doParallel) || install.packages("doParallel")
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
  
  ### parallel backend to use multiple cores
  cores = detectCores()
  print(paste("Randomizing on ", cores[1] - 1, "cores."))
  cl <- makeCluster(cores[1] - 1)
  registerDoParallel(cl)
  clusterCall(cl, function() library(spdep))
  clusterExport(cl, list("LOSH"), envir = environment())
  on.exit(stopCluster(cl))
  pvals <- foreach(i=1:length(x)) %dopar% {
    bootstrap <- numeric(length = nsim + 1)
    for(j in 1:nsim) {
      x_rand <- append(sample(x[-i]), x[i], (i-1))
      bootstrap[j] <- LOSH(x_rand, listw, a, FALSE, zero.policy, na.action, spChk)[i,"Hi"]
    }
    bootstrap[[nsim + 1]] <- res[i, "Hi"]
    rankboot <- rank(unlist(bootstrap))
    xrank <- rankboot[length(bootstrap)]
    diff <- nsim - xrank
    diff <- ifelse(diff > 0, diff, 0)
    pval <- punif((diff + 1)/(nsim + 1))
    if (!is.finite(pval) || pval < 0 || pval > 1) 
      warning("Out-of-range p-value: reconsider test arguments")
    pval
  }
  res[,"Pr()"] <- unlist(pvals)
  res[,"Pr()"] <- p.adjustSP(res[,"Pr()"], listw$neighbours, method = p.adjust.method)
  if (!is.null(rn)) 
    rownames(res) <- rn
  if (!is.null(na.act)) 
    attr(res, "na.action") <- na.act
  class(res) <- c("htest", "mc.sim", "LOSH", class(res))
  res
}
