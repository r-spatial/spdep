# Copyright 2001-2019 by Roger Bivand and Elias Krainski
#


plot.nb <- function(x, coords, col="black", points=TRUE, add=FALSE, 
	arrows=FALSE, length=0.1, xlim=NULL, ylim=NULL, ...) {
	nb <- x

        if (inherits(coords, "sfc")) {
            if (!inherits(coords, "sfc_POINT")) {
#                    if (inherits(coords, "sfc_POLYGON") || #144
#                        inherits(coords, "sfc_MULTIPOLYGON")) 
#
                 if (all(st_dimension(coords) == 2))
                        coords <- st_point_on_surface(coords)
                    else stop("Point-conforming geometries required")
                }
            if (attr(coords, "n_empty") > 0L) 
                stop("Empty geometries found")
            coords <- sf::st_coordinates(coords)
        } else if (inherits(coords, "Spatial")) {
            coords <- coordinates(coords)
        }
        
        stopifnot(length(nb) == nrow(coords))
	sym <- is.symmetric.nb(nb, verbose = FALSE, force = FALSE)
	x <- coords[,1, drop=TRUE]
	y <- coords[,2, drop=TRUE]
	n <- length(nb)
	if (n < 1) stop("non-positive number of entities")
	if (!add) {
		plot.new()
		if (is.null(xlim)) xlim <- range(x)
		if (is.null(ylim)) ylim <- range(y)
        	plot.window(xlim = xlim, ylim = ylim, log="", asp=1)
	}
	cardnb <- card(nb)
	if (length(col) < n) col <- rep(col[1], n)
#	for (i in 1:n) {
#		if (cardnb[i] > 0) {
#       		inb <- nb[[i]]
#        		for (j in inb) {
#				if (sym) {
#					lines(c(x[i], x[j]), c(y[i], y[j]),
#						col=col[i], ...)
#				} else {
#					if (arrows) 
#						arrows(x[i], y[i], x[j], y[j], 
#						col=col[i], length=length, ...)
#					else lines(c(x[i], x[j]), c(y[i], y[j]),
#						col=col[i], ...)
#				}
#
#			}
#		}
#	}

# Elias Krainski Tue, 21 May 2013

   i <- rep(1:n, cardnb)
   j <- unlist(nb)
   if (arrows)
     arrows(x[i], y[i], x[j], y[j], col=col[i], length = length, ...)
   else segments(x[i], y[i], x[j], y[j], col=col[i], ...)


	if (points) points(x, y, ...)
}

plot.listw <- function(x, coords, col="black", points=TRUE, add=FALSE, 
	arrows=FALSE, length=0.1, xlim=NULL, ylim=NULL, ...) {
	plot.nb(x$neighbours, coords=coords, col=col, points=points, add=add, 
	arrows=arrows, length=length, xlim=xlim, ylim=ylim, ...)
}
