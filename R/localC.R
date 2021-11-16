#' @title Compute Local Geary statistic
#'
#' @param x a numeric vector, numeric matrix, or list. See details for more.
#' @param listw a \code{listw} object created for example by \code{nb2listw}
#' @param spChk hould the data vector names be checked against the spatial objects for identity integrity, TRUE, or FALSE, default NULL to use \code{get.spChkOption()}
#'
#' @description
#'
#' The Local Geary is a local adaptation of Geary's C statistic of spatial autocorrelation. The Local Geary uses squared differences to measure dissimilarity unlike the Local Moran. Low values of the Local Geary indicate positive spatial autocorrelation and large, negative.
#'

#'
#' @details
#'
#' The Local Geary can be extended to a multivariate context. When \code{x} is a numeric vector, the univariate Local Geary will be calculated. To calculate the multivariate Local Moran provide either a list or a matrix. When \code{x} is a list, each element must be a numeric vector of the same length and of the same length as the neighbours in \code{listw}. In the case that \code{x} is a matrix the number of rows must be the same as the length of the neighbours in \code{listw}.
#'
#' While not required in the univariate context, the standardized Local Geary is calculated. The multivariate Local Geary is \emph{always} standardized.
#'
#' The univariate Local Geary is calculated as \eqn{c_i = \sum_j w_{ij}(x_i - x_j)^2} and the multivariate Local Geary is calculated as \eqn{c_{k,i} = \sum_{v=1}^{k} c_{v,i}} as described in Anselin (2019).
#'
#' @author Josiah Parry, \email{josiah.parry@gmail.com}
#' @references {Anselin, L. (1995), Local Indicators of Spatial Association—LISA. Geographical Analysis, 27: 93-115. \doi{10.1111/j.1538-4632.1995.tb00338.x}}
#'
#' {Anselin, L. (2019), A Local Indicator of Multivariate Spatial Association: Extending Geary's c. Geogr Anal, 51: 133-150. \doi{10.1111/gean.12164}}
#'
#' @examples
#' columbus <- st_read(system.file("shapes/columbus.shp", package="spData")[1], quiet=TRUE)
#' col.gal.nb <- read.gal(system.file("weights/columbus.gal", package="spData")[1])
#' localC(columbus$CRIME, nb2listw(col.gal.nb))
#' localC_perm(columbus$CRIME, nb2listw(col.gal.nb))
#' @export
localC <-  function(x, listw, spChk=NULL) {

  n <- length(listw$neighbours)

  # checks for list and matrix elements
  if (is.list(x) & !length(unique(lengths(x))) == 1)
    stop(paste("Elements of", deparse(substitute(x)), "must be of equal length"))

  # maxtrix check
  if (is.matrix(x) && !nrow(x) == n)
    stop("Different numbers of observations")

  if(is.vector(x) & is.numeric(x) & n != length(x))
    stop("Different numbers of observations")

  # check listw object
  if (!inherits(listw, "listw"))
    stop(paste(deparse(substitute(listw)), "is not a listw object"))

  # check missing values
  if (any(is.na(x))) stop(paste("NA in ", deparse(substitute(x))))

  if (is.null(spChk)) spChk <- get.spChkOption()
  if (is.numeric(x) & spChk && !chkIDs(x, listw))
    stop("Check of data and weights ID integrity failed")

  if (inherits(x, "list")) {
    x <- scale(Reduce(cbind, x))
  }

  if (inherits(x, "matrix"))
    res <- rowSums(apply(scale(x), 2, localC_calc, listw)) / ncol(x)

  if (inherits(x, "numeric")) res <- localC_calc(scale(x), listw)

  # adding attributes as done in localG
  attr(res, "call") <- match.call()
  class(res) <- c("localC", "numeric")

  res

}




#' @rdname localC
#' @export
localC_perm <- function(x, listw, nsim = 499, spChk = NULL) {
  # checks are inherited from localC no need to implement
  obs <- localC(x, listw)

  if (inherits(x, "list")) {
    x <- Reduce(cbind, x)
    reps <- replicate(nsim, localC(x[sample.int(nrow(x)),], listw))
  }

  if (inherits(x, "matrix")) {
    reps <- replicate(nsim, localC(x[sample.int(nrow(x)),], listw))
  }

  if (is.vector(x) & is.numeric(x)) {
    reps <- replicate(nsim, localC(x[sample.int(length(x))], listw))
  }

  pseudo_p <- (rowSums(reps <= obs)+ 1)/ (nsim + 1)


  attr(res, "call") <- match.call()
  attr(res, "pseudo-p") <- pseudo_p
  class(res) <- c("localC", "numeric")

  res
}

localC_calc <- function(x, listw) {
  xij <- sapply(listw$neighbours, FUN = function(listw) x[listw])
  res <- mapply(function(x, j, wi) sum(wi * (j - x)^2), x, xij, listw$weights)
  res
}
