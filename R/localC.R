#' @title Compute Local Geary statistic
#'
#' @param x a numeric vector the same length as the neighbours list in listw
#' @param listw a \code{listw} object created for example by \code{nb2listw}
#' @param spChk hould the data vector names be checked against the spatial objects for identity integrity, TRUE, or FALSE, default NULL to use \code{get.spChkOption()}
#'
#' @description
#'
#' The Local Geary is a local adaptation of Geary's C statistic of spatial autocorrelation. The Local Geary uses squared differences to measure dissimilarity unlike the Local Moran. Low values of the Local Geary indicate positive spatial autocorrelation and large, negative.
#'
#' @author Josiah Parry, \email{josiah.parry@gmail.com}
#' @references {Anselin, L. (1995), Local Indicators of Spatial Association—LISA. Geographical Analysis, 27: 93-115. \doi{10.1111/j.1538-4632.1995.tb00338.x}}
#'
#' {Anselin, L. (2019), A Local Indicator of Multivariate Spatial Association: Extending Geary's c. Geogr Anal, 51: 133-150. \doi{10.1111/gean.12164}}
#' @export
localC <-  function(x, listw, spChk=NULL) {
  # using checks from localG
  if (!inherits(listw, "listw"))
    stop(paste(deparse(substitute(listw)), "is not a listw object"))
  if (!is.numeric(x))
    stop(paste(deparse(substitute(x)), "is not a numeric vector"))
  stopifnot(is.vector(x))
  if (any(is.na(x))) stop(paste("NA in ", deparse(substitute(x))))
  n <- length(listw$neighbours)
  if (n != length(x))stop("Different numbers of observations")
  if (is.null(spChk)) spChk <- get.spChkOption()
  if (spChk && !chkIDs(x, listw))
    stop("Check of data and weights ID integrity failed")

  x <- scale(x)

  xij <- sapply(listw$neighbours, FUN = function(listw) x[listw])
  res <- mapply(function(x, j, wi) sum(wi * (j - x)^2), x, xij, listw$weights)

  # adding attributes as done in localG
  attr(res, "call") <- match.call()
  class(res) <- c("localC", "numeric")

  res

}


#' @rdname localC
#' @export
localC_perm <- function(x, listw, nsim = 499, spChk = NULL) {
  # checks are inherited from localC no need to implement
  res<- localC(x, listw)

  reps <- replicate(nsim, localC(x[sample.int(length(x))], listw))

  pseudo_p <- (rowSums(reps <= obs)+ 1)/ (nsim + 1)


  attr(res, "call") <- match.call()
  attr(res, "pseudo-p") <- pseudo_p
  class(res) <- c("localC", "numeric")

  res
}

