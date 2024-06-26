# Copyright 2001-2010 by Roger Bivand
#


tri2nb <- function(coords, row.names = NULL) {
#	require("tripack")
#	require("deldir")
        if (inherits(coords, "SpatialPoints")) {
            if (!is.na(is.projected(coords)) && !is.projected(coords)) {
                warning("tri2nb: coordinates should be planar")
            }
            coords <- coordinates(coords)
        } else if (inherits(coords, "sfc")) {
            if (!inherits(coords, "sfc_POINT"))
                stop("Point geometries required")
            if (attr(coords, "n_empty") > 0L) 
                stop("Empty geometries found")
            if (!is.na(sf::st_is_longlat(coords)) && sf::st_is_longlat(coords))
                warning("tri2nb: coordinates should be planar")
            coords <- sf::st_coordinates(coords)
        }

	n <- nrow(coords)
	if (n < 3) stop("too few coordinates")
#	left <- function(x) {
#		res <- (x[3,1]-x[2,1])*(x[1,2]-x[2,2]) >= 
#			(x[1,1]-x[2,1])*(x[3,2]-x[2,2])
#		res
#	}
#	if (left(coords[1:3,])) stop("first three coordinates collinear")
    	if (!is.null(row.names)) {
		if(length(row.names) != n)
            		stop("row.names wrong length")
		if (length(unique(row.names)) != length(row.names))
	    		stop("non-unique row.names given")
    	}
    	if (is.null(row.names)) row.names <- as.character(1:n)
        stopifnot(!anyDuplicated(coords))
#	tri <- tri.mesh(x=coords[,1], y=coords[,2])
        tri <- deldir::deldir(x=coords[, 1, drop=TRUE], y=coords[,2, drop=TRUE])
        from <- c(tri$delsgs[,5], tri$delsgs[,6])
        to <- c(tri$delsgs[,6], tri$delsgs[,5])
        df <- data.frame(from=as.integer(from), to=as.integer(to), weight=1)
        attr(df, "n") <- tri$n.data
        class(df) <- c(class(df), "spatial.neighbour")
        df1 <- df[order(df$from),]
        nb <- sn2listw(df1, style="B")$neighbours
#	nb <- neighbours(tri)
 	attr(nb, "region.id") <- row.names
	class(nb) <- "nb"
	attr(nb, "tri") <- TRUE
	attr(nb, "call") <- match.call()
	nb <- sym.attr.nb(nb)
        NE <- n + sum(card(nb))
        if (get.SubgraphOption() && get.SubgraphCeiling() > NE) {
            ncomp <- n.comp.nb(nb)
            attr(nb, "ncomp") <- ncomp
            if (ncomp$nc > 1) warning("neighbour object has ", ncomp$nc, " sub-graphs")
        }
	nb
}

