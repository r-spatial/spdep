library(spdep)
# https://github.com/r-spatial/spdep/issues/20
if (require("sp", quietly=TRUE)) {
GT <- GridTopology(c(1, 1), c(1, 1), c(10, 50))
SPix <- as(SpatialGrid(GT), "SpatialPixels")
nb_rook_cont <- poly2nb(as(SPix, "SpatialPolygons"), queen=FALSE)
nb_rook_dist <- dnearneigh(coordinates(SPix), 0, 1.01)
expect_true(all.equal(nb_rook_cont, nb_rook_dist, check.attributes=FALSE))
## [1] TRUE
t.nb <- cell2nb(GT, type='rook', legacy=TRUE)
expect_false(isTRUE(all.equal(nb_rook_cont, t.nb, check.attributes=FALSE)))
## [1] FALSE
t.nb <- cell2nb(GT, type='rook')
expect_true(isTRUE(all.equal(nb_rook_cont, t.nb, check.attributes=FALSE)))
## [1] TRUE
# https://github.com/r-spatial/spdep/issues/55
c <- 20:21
r <- 11:12
set.seed(1)
for (i in c) {
  for (j in r) {
    GT <- GridTopology(c(1, 1), c(1, 1), c(i, j))
    print(GT)
    SPix <- as(SpatialGrid(GT), "SpatialPixels")
    x <- rnorm(prod(slot(GT, "cells.dim")))
    SPixDF <- SpatialPixelsDataFrame(SPix, data=data.frame(x=x))
    SGDF <- as(SPixDF, "SpatialGridDataFrame")
    SPolDF <- as(SPixDF, "SpatialPolygonsDataFrame")
    nb_rook_cont <- poly2nb(SPolDF, queen=FALSE)
    M_cont <- moran.test(SPolDF$x, nb2listw(nb_rook_cont))
    nb_rook_dist <- dnearneigh(coordinates(SPix), 0, 1.01)
    M_dist <- moran.test(SPixDF$x, nb2listw(nb_rook_dist))
    print(expect_true(isTRUE(all.equal(nb_rook_cont, nb_rook_dist,
      check.attributes=FALSE))))
    print(expect_equal(M_cont$estimate[1], M_dist$estimate[1]))
    nb_legacy <- cell2nb(GT, type='rook', legacy=TRUE)
    print(expect_false(isTRUE(all.equal(nb_rook_cont, nb_legacy,
      check.attributes=FALSE))))
    ## [1] FALSE
    M_cell_legacy <- moran.test(SGDF$x, nb2listw(nb_legacy))
    print(expect_false(isTRUE(all.equal(M_cont$estimate[1],
      M_cell_legacy$estimate[1]))))
    nb <- cell2nb(GT, type='rook')
    print(expect_true(isTRUE(all.equal(nb_rook_cont, nb,
      check.attributes=FALSE))))
    ## [1] TRUE
    M_cell <- moran.test(SGDF$x, nb2listw(nb))
    print(expect_equal(M_cont$estimate[1], M_cell$estimate[1]))
  }
}
NIN <- data.frame(row = rep(1:11, 22), col = rep(1:22, each = 11))
NIN$id <- paste(NIN$col, NIN$row, sep=":")
xy.rook <- cell2nb(nrow = max(NIN$row), ncol = max(NIN$col), type="rook")
GT <- GridTopology(c(1, 1), c(1, 1), c(max(NIN$col),  max(NIN$row)))
xy.rook1 <- cell2nb(GT, type="rook")
expect_true(isTRUE(all.equal(xy.rook, xy.rook1, check.attributes=FALSE)))
NIN_sf <- st_as_sf(NIN, coords=c("col", "row"))
xy.dist <- dnearneigh(NIN_sf, 0, 1.01, row.names=NIN_sf$id)
expect_false(isTRUE(all.equal(xy.dist, xy.rook, check.attributes=FALSE)))
expect_false(isTRUE(all.equal(attr(xy.rook, "region.id"),
  attr(xy.dist, "region.id"))))
NINa <- data.frame(row = rep(1:11, each = 22), col = rep(1:22, 11))
NINa$id <- paste(NINa$col, NINa$row, sep=":")
xy.rooka <- cell2nb(nrow = max(NINa$row), ncol = max(NINa$col), type="rook")
expect_true(isTRUE(all.equal(xy.rook, xy.rooka, check.attributes=FALSE)))
NIN_sfa <- st_as_sf(NINa, coords=c("col", "row"))
xy.dista <- dnearneigh(NIN_sfa, 0, 1.01, row.names=NIN_sfa$id)
expect_true(isTRUE(all.equal(xy.dista, xy.rooka, check.attributes=FALSE)))
expect_true(isTRUE(all.equal(attr(xy.rooka, "region.id"),
  attr(xy.dista, "region.id"))))
}

