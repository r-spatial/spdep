nbcosts <- function(nb, data, method=c("euclidean", "maximum", "manhattan",
                                "canberra", "binary", "minkowski",
                                "mahalanobis"), p=2, cov, inverted=FALSE) {
#  if ((!require(parallel)) | (length(nb)<300))
#    clist <- lapply(1:length(nb), function(i)
#                    nbcost(data, i, nb[[i]], method,
#                           p, cov, inverted))
#  else {
#    if (.Platform$OS.type == "windows") {
#      cl <- makeCluster(getOption("cl.cores", 2))
#      clusterEvalQ(cl, library(spdep))
    if (any(card(nb) == 0L)) stop("nbcosts: no-neighbour nodes")
    nc <- n.comp.nb(nb)$nc
    if (nc > 1) stop("nbcosts:", nc, "disjoint connected subgraphs")
    if (missing(cov)) cov <- NULL
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
    if (length(nb)<300) parallel <- "no"
    
    if (parallel == "snow") {
      if (requireNamespace("parallel", quietly = TRUE)) {
#        require(parallel)
        sI <- parallel::splitIndices(length(nb), length(cl))
         env <- new.env()
         assign("nb", nb, envir=env)
         assign("data", data, envir=env)
         assign("method", method, envir=env)
         assign("p", p, envir=env)
         assign("cov", cov, envir=env)
         assign("inverted", inverted, envir=env)
         parallel::clusterExport(cl, varlist=c("nb", "data", "method", "p", "cov",
             "inverted"), envir=env)
         out <- parallel::clusterApply(cl, x = sI, fun=lapply, function(i) {
 	     nbcost(data, i, nb[[i]], method, p, cov, inverted)})
        clist <- do.call("c", out)
        rm(env)
      } else {
        stop("parallel not available")
      }
    } else if (parallel == "multicore") {
      if (requireNamespace("parallel", quietly = TRUE)) {
#        require(parallel)
        sI <- parallel::splitIndices(length(nb), ncpus)
        out <- parallel::mclapply(sI, FUN=lapply, function(i) {nbcost(data, i, nb[[i]],
            method, p, cov, inverted)}, mc.cores=ncpus)
        clist <- do.call("c", out)
      } else {
        stop("parallel not available")
      }
    } else {
      clist <- lapply(1:length(nb),
                   function(i) nbcost(data, i, nb[[i]], method,
                           p, cov, inverted))
    }
    attr(clist, "call") <- match.call()
    attr(clist, "class") <- "nbdist"
    return(clist)
}

nbcost <- function(data, id, id.neigh,
                   method=c("euclidean", "maximum", "manhattan",
                     "canberra", "binary", "minkowski",
                     "mahalanobis"), p=2, cov, inverted=FALSE) {
  if (is.function(method))
    return(method(data, id, id.neigh))
  else {
    method <- match.arg(method)
    data <- as.matrix(data)
    if (method=="mahalanobis")
      return(mahalanobis(data[id.neigh,,drop=FALSE], data[id,,drop=FALSE],
cov, inverted))
    else
      return(dist(rbind(data[id,,drop=FALSE], data[id.neigh,,drop=FALSE]),
method=method,
                p=p)[1:length(id.neigh)])
  }
}


