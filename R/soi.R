# Copyright 2001 by Nicholas Lewin-Koh, modified RSB 2016-05-31
#


#soi.graph <- function(tri.nb, coords){
#  x <- coords
#  if (!is.matrix(x)) stop("Data not in matrix form")
#  if (any(is.na(x))) stop("Data cannot include NAs")
#  np<-length(tri.nb)
#  noedges<-0
#  rad<-nearneigh<-rep(0,np)
#  neigh<-unlist(tri.nb)  
#  noneigh<-unlist(lapply(tri.nb,length))
#  g1<-g2<-rep(0,sum(noneigh))
#  storage.mode(x) <- "double"
# called Computational Geometry in C functions banned by Debian admins
#  answ<-.C("compute_soi", np=as.integer(np), from=as.integer(g1),
#     to=as.integer(g2), nedges=as.integer(noedges),
#     notri.nb=as.integer(noneigh), tri.nb=as.integer(neigh),
#     nn=as.integer(nearneigh), 
#     circles=as.double(rad), x=x[,1], y=x[,2],
#     PACKAGE="spdep")
#  answ$from<-answ$from[1:answ$nedges]
#  answ$to<-answ$to[1:answ$nedges]
#  answ<-list(np=answ$np,nedges=answ$nedges,
#             from=answ$from,to=answ$to,circles=answ$circ)
#  attr(answ, "call") <- match.call()
#  class(answ)<-c("Graph","SOI")
#  answ
#}

soi.graph <- function(tri.nb, coords, quadsegs=10){
  obj <- NULL
  if (!is.matrix(coords)) obj <- coords
  if (inherits(obj, "SpatialPoints")) {
      obj <- sf::st_geometry(sf::st_as_sf(obj))
  } 
  if (inherits(obj, "sfc")) {
      if (!inherits(obj, "sfc_POINT"))
           stop("Point geometries required")
      if (attr(obj, "n_empty") > 0L) 
           stop("Empty geometries found")
      if (!is.na(sf::st_is_longlat(obj)) && sf::st_is_longlat(obj))
           warning("tri2nb: coordinates should be planar")
      coords <- sf::st_coordinates(obj)
  }
  if (!is.matrix(coords)) stop("Data not in matrix form")
  if (any(is.na(coords))) stop("Data cannot include NAs")
  stopifnot(length(tri.nb) == nrow(coords))
  if (requireNamespace("dbscan", quietly = TRUE)) {
    dists_1 <- dbscan::kNN(coords, k=1)$dist[,1]
  } else {
    stop("dbscan required")
  }

  if (is.null(obj)) obj <- st_geometry(st_as_sf(as.data.frame(coords),
    coords=1:2))
  if(inherits(obj, "sfc")) {
    bobj <- st_buffer(obj, dist=dists_1, nQuadSegs=quadsegs)
    gI <- st_intersects(bobj)
  } 
  gI_1 <- lapply(1:length(gI), function(i) {
    x <- gI[[i]][match(tri.nb[[i]], gI[[i]])]; x[!is.na(x)]
  })
  answ <- list(np=length(tri.nb),nedges=length(unlist(gI_1)),
    from=rep(1:length(tri.nb), times=sapply(gI_1, length)), to=unlist(gI_1),
    circles=dists_1)
  attr(answ, "call") <- match.call()
  class(answ) <- c("Graph","SOI")
  answ
}

