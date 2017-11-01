\name{set.mcOption}
\alias{set.ClusterOption}
\alias{get.ClusterOption}
\alias{set.mcOption}
\alias{get.coresOption}
\alias{set.coresOption}
\alias{get.mcOption}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Options for parallel support}
\description{
  Provides support for the use of parallel computation in the parallel package. 
}
\usage{
set.mcOption(value)
get.mcOption()
set.coresOption(value)
get.coresOption()
set.ClusterOption(cl)
get.ClusterOption()
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{value}{valid replacement value}
  \item{cl}{a cluster object created by \code{makeCluster} in \pkg{parallel}}
}
\details{
 Options in the spdep package are held in an environment local to the package namespace and not exported. Option values are set and retrieved with pairs of access functions, get and set. The \code{mc} option is set by default to FALSE on Windows systems, as they cannot fork the R session; by default it is TRUE on other systems, but may be set FALSE. If \code{mc} is FALSE, the \code{Cluster} option is used: if \code{mc} is FALSE and the \code{Cluster} option is NULL no parallel computing is done, or the \code{Cluster} option is passed a \dQuote{cluster} object created by the parallel or snow package for access without being passed as an argument. The \code{cores} option is set to NULL by default, and can be used to store the number of cores to use as an integer. If \code{cores} is NULL, facilities from the parallel package will not be used.
}
\value{
  The option access functions return their current settings, the assignment functions usually return the previous value of the option.
}
\note{An extended example is shown in the documentation of \code{\link{aple.mc}}, including treatment of seeding of RNG for multicore/cluster.}

\author{Roger Bivand \email{Roger.Bivand@nhh.no}}

 
\examples{
ls(envir=spdep:::.spdepOptions)
library(parallel)
nc <- detectCores(logical=FALSE)
nc
# set nc to 1L here
if (nc > 1L) nc <- 1L
#nc <- ifelse(nc > 2L, 2L, nc)
coresOpt <- get.coresOption()
coresOpt
if (!is.na(nc)) {
 invisible(set.coresOption(nc))
 print(exists("aple.mc"))
 if(.Platform$OS.type == "windows") {
# forking not permitted on Windows - start cluster
  print(get.mcOption())
  cl <- makeCluster(get.coresOption())
  print(clusterEvalQ(cl, exists("aple.mc")))
  set.ClusterOption(cl)
  clusterEvalQ(get.ClusterOption(), library(spdep))
  print(clusterEvalQ(cl, exists("aple.mc")))
  clusterEvalQ(get.ClusterOption(), detach(package:spdep))
  set.ClusterOption(NULL)
  print(clusterEvalQ(cl, exists("aple.mc")))
  stopCluster(cl)
 } else {
  mcOpt <- get.mcOption()
  print(mcOpt)
  print(mclapply(1:get.coresOption(), function(i) exists("aple.mc"),
   mc.cores=get.coresOption()))
  invisible(set.mcOption(FALSE))
  cl <- makeCluster(nc)
  print(clusterEvalQ(cl, exists("aple.mc")))
  set.ClusterOption(cl)
  clusterEvalQ(get.ClusterOption(), library(spdep))
  print(clusterEvalQ(cl, exists("aple.mc")))
  clusterEvalQ(get.ClusterOption(), detach(package:spdep))
  set.ClusterOption(NULL)
  print(clusterEvalQ(cl, exists("aple.mc")))
  stopCluster(cl)
  invisible(set.mcOption(mcOpt))
 }
 invisible(set.coresOption(coresOpt))
}
}
\keyword{ spatial }
