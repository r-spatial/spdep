prunecost <- function(edges, data,
                      method=c("euclidean", "maximum", "manhattan",
                        "canberra", "binary", "minkowski",
                        "mahalanobis"), p=2, cov, inverted=FALSE) {
  sswt <- ssw(data, unique(as.integer(edges)),
              method, p, cov, inverted)
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
    if (nrow(edges)<300) parallel <- "no"
#    if (parallel == "snow") {
#        parallel <- "no"
#        warning("no parallel calculations available")
#    }
    if (parallel == "snow") {
      if (requireNamespace("parallel", quietly = TRUE)) {
#        require(parallel)
        sI <- parallel::splitIndices(nrow(edges), length(cl))
#    if (.Platform$OS.type == "windows") {
#      cl <- makeCluster(getOption("cl.cores", 2))
#      clusterEvalQ(cl, library(spdep))
        sswp <- do.call("c", parallel::parLapply(cl, sI, sapply, function(i) {
            pruned.ids <- prunemst(rbind(edges[i, ], edges[-i, ]),
                             only.nodes=TRUE)
            sum(sapply(pruned.ids, function(j) 
                 ssw(data, j, method, p, cov, inverted)))
        }))
      } else {
        stop("parallel not available")
      }
    } else if (parallel == "multicore") {
      if (requireNamespace("parallel", quietly = TRUE)) {
#        require(parallel)
        sI <- parallel::splitIndices(nrow(edges), ncpus)
        out <- parallel::mclapply(sI, sapply, function(i) { 
            pruned.ids <- prunemst(rbind(edges[i, ], edges[-i, ]),
                             only.nodes=TRUE)
            sum(sapply(pruned.ids, function(j)
                 ssw(data, j, method, p, cov, inverted)))        
            }, mc.cores=ncpus)
        sswp <- do.call("c", out)
      } else {
        stop("parallel not available")
      }
    } else {
        sswp <- sapply(1:nrow(edges), function(i) {
            pruned.ids <- prunemst(rbind(edges[i, ], edges[-i, ]),
                           only.nodes=TRUE)
            sum(sapply(pruned.ids, function(j)
               ssw(data, j, method, p, cov, inverted)))
        })
    }
    return(sswt - sswp)
}

ssw <- function(data, id, method=c("euclidean", "maximum", "manhattan",
                            "canberra", "binary", "minkowski",
                            "mahalanobis"), p=2, cov, inverted=FALSE) {
  if (is.function(method))
    return(method(data, id))
  else { 
    method <- match.arg(method)
    data <- as.matrix(data)
    if (method=="mahalanobis")
      return(sum(mahalanobis(data[id,,drop=FALSE],
                             colMeans(data[id,,drop=FALSE]),
                             cov, inverted)))
    else
      return(sum(dist(rbind(colMeans(data[id,,drop=FALSE]),
                            data[id,,drop=FALSE]),
                      method, p=p)[1:length(id)]))
  }
}

