localGS <- function (x, listw, dmin, dmax, attr, longlat = NULL) {
  if (!is.numeric(x[[attr]]))
    stop(paste(deparse(substitute(x)), "is not a numeric vector"))
  if (any(is.na(x[[attr]]))) 
    stop(paste("NA in ", deparse(substitute(x))))
  if (!is.numeric(dmin) || !is.numeric(dmax) || !(dmax > dmin)) 
    stop("dmax must be larger than dmin and both must be numeric")
  
  if (inherits(x, "Spatial")) {
    sf <- FALSE
    if (!is.numeric(coordinates(x))) stop("Coordinates non-numeric")
    if (!is.matrix(coordinates(x))) stop("Coordinates not in matrix form")
    if (any(is.na(coordinates(x)))) stop("Coordinates include NAs")
    if ((is.null(longlat) || !is.logical(longlat)) 
        && !is.na(is.projected(x)) && !is.projected(x)) {
      longlat <- TRUE
    } else longlat <- FALSE
  } else {
    sf <- TRUE
  if (!is.numeric(st_coordinates(x))) stop("Coordinates non-numeric")
    if (!is.matrix(st_coordinates(x))) stop("Coordinates not in matrix form")
    if (any(is.na(st_coordinates(x)))) stop("Coordinates include NAs")   
    if (inherits(x, "sf")) {
      if (is.null(row.names)) row.names <- row.names(x)
    }
    if (inherits(x, "sfc")) {
      if ((is.null(longlat) || !is.logical(longlat))
          && !is.na(sf::st_is_longlat(x)) && sf::st_is_longlat(x)) {
        longlat <- TRUE
      } else longlat <- FALSE
    }
  }
  if (is.null(longlat) || !is.logical(longlat)) longlat <- FALSE
  if (longlat) {
    if (!.ll_sanity(bbox(x)))
      warning("Coordinates are not geographical: longlat argument wrong")
  }
  
  if(!sf) {
    phi <- suppressWarnings(nb2listw(dnearneigh(coordinates(x), dmin, dmax, bounds=c("GE", "LE"), longlat = longlat), style="B", zero.policy = T))
    } else {
    phi <- suppressWarnings(nb2listw(dnearneigh(st_centroid(x), dmin, dmax, bounds=c("GE", "LE"), longlat = longlat), style = "B", zero.policy = T))
  }
      
  n <- length(x[[attr]])
  phi_times_f <- vector("list", n) 
  res <- rep(NA, n)
  Ai <- 0
  Wi <- 0

	for(i in 1:n) {
	    if(length(phi$weights[[i]]) > 0) {
	        phi_times_f[[i]] <- x[[attr]][i] + c(x[[attr]][phi$neighbours[[i]]])
	    } else {
	        phi_times_f[[i]] <- 0
	    }
	}
	
	gamma <- sum(unlist(phi_times_f)) / 2
	gamma2 <- sum(unlist(phi_times_f)^2) / 2
	big_phi <- sum(unlist(phi$weights)) / 2
	
	if (big_phi < 1) 
	    stop("There are no links corresponding to the selected geometric scale. Please adjust your scale band.")
	
  for(i in 1:n) {
    Ai <- 0
    Wi <- 0
    if(is.null(phi$weights[[i]]))
	    next
    relevant_neighbours <- sort(c(i, phi$neighbours[[i]]))
    for(j in relevant_neighbours) {
      if(j == i)
        next
      for(k in relevant_neighbours) {
        if((k < j) && (k %in% phi$neighbours[[j]])) {
          if(!is.na(match(j, listw$neighbours[[i]])) && !is.na(match(k, listw$neighbours[[i]]))) {
            wij <- listw$weights[[i]][match(j, listw$neighbours[[i]])]
            wik <- listw$weights[[i]][match(k, listw$neighbours[[i]])]
            Ai <- Ai + (wij * wik * phi_times_f[[j]][match(k, phi$neighbours[[j]])])
            Wi <- Wi + (wij * wik)
          }
        }
      }
    }
    numerator <- Ai - ( (Wi / big_phi) * gamma)
    denominator <- sqrt( (Wi / big_phi) * gamma2 + ((Wi * (Wi - 1)) / (big_phi * (big_phi - 1))) * (gamma^2 - gamma2) - ((Wi / big_phi) * gamma)^2)
    if(!is.nan(numerator) && !is.nan(denominator))
      if(denominator != 0)
        res[i] <- numerator / denominator
  }
  res
}
