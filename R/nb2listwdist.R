nb2listwdist <- function(neighbours, x, type="idw", style="raw", alpha = 1, dmax = NULL, longlat = NULL, zero.policy=NULL)
{
  if(!inherits(neighbours, "nb")) stop("Not a neighbours list")
  if (is.null(zero.policy))
    zero.policy <- get("zeroPolicy", envir = .spdepOptions)
  stopifnot(is.logical(zero.policy))
  if (!(type %in% c("idw", "dpd", "exp")))
    stop(paste("type", type, "invalid"))
  if (inherits(x, "Spatial")) {
    sf <- FALSE
    if ((is.null(longlat) || !is.logical(longlat)) 
        && !is.na(is.projected(x)) && !is.projected(x)) {
      longlat <- TRUE
    } else longlat <- FALSE
    if (!is.numeric(coordinates(x))) stop("Coordinates non-numeric")
    if (!is.matrix(coordinates(x))) stop("Coordinates not in matrix form")
    if (any(is.na(coordinates(x)))) stop("Coordinates include NAs")
  } else {
    sf <- TRUE
    if (inherits(x, "sf"))
      if (is.null(row.names)) row.names <- row.names(x)
    if (inherits(x, "sfc")) {
      if ((is.null(longlat) || !is.logical(longlat)) 
          && !is.na(sf::st_is_longlat(x)) && sf::st_is_longlat(x)) {
        longlat <- TRUE
      } else longlat <- FALSE
    }
    if (!is.numeric(st_coordinates(x))) stop("Coordinates non-numeric")
    if (!is.matrix(st_coordinates(x))) stop("Coordinates not in matrix form")
    if (any(is.na(st_coordinates(x)))) stop("Coordinates include NAs")
  } 
  if (is.null(longlat) || !is.logical(longlat)) longlat <- FALSE
  if (longlat) {
    if (!.ll_sanity(bbox(x)))
      warning("Coordinates are not geographical: longlat argument wrong")
  }
  
  n <- length(neighbours)
  if (n < 1) stop("non-positive number of entities")
  cardnb <- card(neighbours)
  if (!zero.policy)
    if (any(cardnb == 0)) stop("Empty neighbour sets found")
  vlist <- vector(mode="list", length=n)
  glist <- vector(mode="list", length=n)
  
  if (!sf) {
    for (i in 1:n)
      if(cardnb[i] > 0) {
        gx <- geometry(x)
        if(longlat)
          glist[[i]] <- as.numeric(spDists(gx[i], gx[neighbours[[i]]], longlat = TRUE)) * 1000
        else
          glist[[i]] <- as.numeric(spDists(geometry(x)[i], geometry(x)[neighbours[[i]]], longlat = FALSE))
        mode(glist[[i]]) <- "numeric"
      }
  } else {
    if(longlat)
      coords_type = "Great Circle"
    else
      coords_type = "Euclidean"
    gx <- st_geometry(x)
    for (i in 1:n)
      if(cardnb[i] > 0) {
        glist[[i]] <- as.numeric(st_distance(gx[i], gx[neighbours[[i]]], which = coords_type))
        mode(glist[[i]]) <- "numeric"
      }
  }
  
  attr(vlist, "mode") <- "distance"
  attr(vlist, as.character(style)) <- TRUE
  
  if (type == "idw") {
    for (i in 1:n) {
      if (cardnb[i] > 0) {
        vlist[[i]] <- glist[[i]]^((-1) * alpha)
        if(!is.null(dmax))
          if(dmax > 0)
            vlist[[i]][which(glist[[i]] > dmax)] <- 0
      }
    }
    max_finite <- max(is.finite(unlist(vlist)))
    for(i in 1:n) {
      vlist[[i]][which(is.infinite(vlist[[i]]))] <- max_finite
    }
  }
  
  if (type == "exp") {
    for (i in 1:n) {
      if (cardnb[i] > 0) {
        vlist[[i]] <- exp(glist[[i]] * ((-1) * alpha))
        if(!is.null(dmax))
          if(dmax > 0)
            vlist[[i]][which(glist[[i]] > dmax)] <- 0
      }
    }
  }
  
  if (type == "dpd") {
    if (is.null(dmax)) stop("DPD weights require a maximum distance threshold")
    if (dmax <= 0) stop("DPD weights require a positive maximum distance threshold")
    for (i in 1:n) {
      if (cardnb[i] > 0) {
        vlist[[i]] <- (1 - (glist[[i]] / dmax)^alpha)^alpha
        vlist[[i]][which(vlist[[i]] < 0)] <- 0
      }
    }
  }
  
  if(style != "raw") {
    
    glist <- vlist
    
    if (zero.policy) {
      eff.n <- n - length(which(cardnb == 0))
      if (eff.n < 1) stop("No valid observations")
    } else eff.n <- n
    
    if (style == "W") {
      d <- unlist(lapply(glist, sum))
      for (i in 1:n) {
        if (cardnb[i] > 0) {
          if (d[i] > 0) vlist[[i]] <- (1/d[i]) * glist[[i]]
          else vlist[[i]] <- 0
        }
      }
      attr(vlist, "comp") <- list(d=d)
    }
    
    if (style == "B") {
      for (i in 1:n) {
        if (cardnb[i] > 0) vlist[[i]] <- as.numeric(I(glist[[i]] > 0))
      }
    }
    
    if (style == "C" || style == "U") {
      D <- sum(unlist(glist))
      if (is.na(D) || !(D > 0))
        stop(paste("Failure in sum of weights:", D))
      for (i in 1:n) {
        if (cardnb[i] > 0) {
          if (style == "C")
            vlist[[i]] <- (eff.n/D) * glist[[i]]
          else
            vlist[[i]] <- (1/D) * glist[[i]]
        }
      }
    }
    
    if (style == "S") {
      glist2 <- lapply(glist, function(x) x^2)
      q <- sqrt(unlist(lapply(glist2, sum)))
      for (i in 1:n) {
        if (cardnb[i] > 0) {
          if (q[i] > 0) glist[[i]] <- (1/q[i]) * glist[[i]]
          else glist[[i]] <- 0 
        }
      }
      Q <- sum(unlist(glist))
      if (is.na(Q) || !(Q > 0))
        stop(paste("Failure in sum of intermediate weights:", Q))
      for (i in 1:n) {
        if (cardnb[i] > 0)
          vlist[[i]] <- (eff.n/Q) * glist[[i]]
      }
      attr(vlist, "comp") <- list(q=q, Q=Q, eff.n=eff.n)
    }
  }
  
  style <- style
  if (!zero.policy)
    if (any(is.na(unlist(vlist))))
      stop ("NAs in coding scheme weights list")
  
  if (style == "minmax") {
    res <- list(style=style, neighbours=neighbours, weights=vlist)
    class(res) <- c("listw", "nb")
    mm <- minmax.listw(res)
    vlist <- lapply(vlist, function(x) (1/c(mm)) * x)
  }
  
  res <- list(style=style, type=type, neighbours=neighbours, weights=vlist)
  class(res) <- c("listw", "nb")
  attr(res, "region.id") <- attr(neighbours, "region.id")
  attr(res, "call") <- match.call()
  if (!is.null(attr(neighbours, "GeoDa")))
    attr(res, "GeoDa") <- attr(neighbours, "GeoDa")
  if (!is.null(attr(res, "GeoDa")$dist)) 
    attr(res, "GeoDa")$dist <- NULL
  res
}
