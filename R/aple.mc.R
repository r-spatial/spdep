aple.mc <- function(x, listw, nsim, override_similarity_check=FALSE,
    useTrace=TRUE) {
    aple.boot <- function(var, i, ...) {
        var <- var[i]
        return(inAple(x=var, ...))
    }
    pre <- preAple(x=x, listw=listw,
        override_similarity_check=override_similarity_check, useTrace=useTrace)
    
    cores <- get.coresOption()
    if (is.null(cores)) {
        parallel <- "no"
    } else {
        parallel <- ifelse (get.mcOption(), "multicore", "snow")
    }
    ncpus <- ifelse(is.null(cores), 1L, cores)
    cl <- NULL
    if (parallel == "snow") {
        cl <- get.ClusterOption()
        if (is.null(cl)) {
            parallel <- "no"
            warning("no cluster in ClusterOption, parallel set to no")
        }
    }
    res <- boot(x, statistic=aple.boot, R=nsim, sim="permutation", pre=pre,
        parallel=parallel, ncpus=ncpus, cl=cl)

    res
}

boot_wrapper_in <- function(cl, nsim) {
      if (requireNamespace("parallel", quietly = TRUE)) {
#        require(rlecuyer)
        rlseed <- get("rlecuyerSeed", envir = .spdepOptions)
        if (storage.mode(rlseed) != "integer") rlseed <- as.integer(rlseed)
        if (length(rlseed) != 6L) rlseed <- rep(12345L, 6)
        parallel::clusterSetRNGStream(cl, iseed=rlseed)
        parallel::clusterEvalQ(cl, library(spdep))
        nnsim <- ceiling(nsim/length(cl))
        nnsim
      } else {
        stop("parallel not available")
      }
}

boot_wrapper_out <- function(lres, mcall) {
        res <- list()
        res$t0 <- lres[[1]]$t0
        res$t <- matrix(c(sapply(lres, function(x) x$t)), ncol=1)
        res$R <- sum(sapply(lres, function(x) x$R))
        res$data <- lres[[1]]$data
        res$seed <- c(sapply(lres, function(x) x$seed))
        res$statistic <- lres[[1]]$statistic
        res$sim <- lres[[1]]$sim
        res$call <- mcall
        res$stype <- lres[[1]]$stype
        res$strata <- lres[[1]]$strata
        class(res) <- "boot"
        res
}
