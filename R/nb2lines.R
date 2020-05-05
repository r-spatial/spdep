# Copyright 2005-7 by Roger Bivand
#

nb2lines <- function(nb, wts, coords, proj4string=NULL, as_sf=FALSE) {

        if (inherits(coords, "sfc")) {
            if (!inherits(coords, "sfc_POINT")) {
                if (inherits(coords, "sfc_POLYGON") || 
                    inherits(coords, "sfc_MULTIPOLYGON")) 
                    coords <- sf::st_point_on_surface(coords)
                else stop("Point-conforming geometries required")
            }
            if (attr(coords, "n_empty") > 0L) 
                stop("Empty geometries found")
            as_sf <- TRUE
            proj4string <- sf::st_crs(coords)
            coords <- sf::st_coordinates(coords)
        } else if (inherits(coords, "Spatial")) {
            proj4string <- slot(coords, "proj4string")
            coords <- coordinates(coords)
        }
	x <- coords[,1]
	y <- coords[,2]
	n <- length(nb)
	if (n < 1) stop("zero length neighbour list")
	ID <- as.character(attr(nb, "region.id"))
	cardnb <- card(nb)
	totlinks <- sum(cardnb)
	ll <- vector(mode="list", length=totlinks)
	df <- data.frame(i=integer(totlinks), j=integer(totlinks),
		i_ID=character(totlinks), j_ID=character(totlinks),
		wt=numeric(totlinks), stringsAsFactors=FALSE)
	line = 1
	for (i in 1:n) {
		if (cardnb[i] > 0) {
        		inb <- nb[[i]]
			if (!missing(wts)) iwts <- wts[[i]]
        		for (j in 1:cardnb[i]) {
				jj <- inb[j]
				xx <- c(x[i], x[jj])
				yy <- c(y[i], y[jj])
				xy <- cbind(xx, yy)
#				ll[[line]] <- cbind(xx, yy)
                                if (as_sf) {
                                  ll[[line]] <- sf::st_linestring(xy, dim="XY")
                                } else {
				  Ll <- list(Line(xy))
				  ll[[line]] <- Lines(Ll, ID=as.character(line))
                                }
				df[line, "i"] <- i
				df[line, "i_ID"] <- ID[i]
				df[line, "j"] <- jj
				df[line, "j_ID"] <- ID[jj]
				if (missing(wts))
				    df[line, "wt"] <- 1
				else
				    df[line, "wt"] <- iwts[j]
				line <- line + 1
			}
		}
	}
	row.names(df) <- as.character(1:(line-1))
        if (as_sf) {
            res <- df
            if (!inherits(proj4string, "crs"))
                proj4string <- sf::st_crs(proj4string)
            sf::st_geometry(res) <- sf::st_as_sfc(ll, crs=proj4string)
        } else {
            if (!inherits(proj4string, "CRS")) proj4string <- CRS(proj4string)
	    res <- SpatialLinesDataFrame(SpatialLines(ll,
                proj4string=proj4string), data=df)
        }
        res
}

listw2lines <- function(listw, coords, proj4string=NULL, as_sf=FALSE) {
	nb2lines(listw$neighbours, listw$weights, coords, proj4string, as_sf)
}

df2sn <- function(df, i="i", i_ID="i_ID", j="j", wt="wt") {
	IDs <- unique(df[c(i, i_ID)])
	res <- df[c(i, j, wt)]
	names(res) <- c("from", "to", "weights")
	attr(res, "n") <- nrow(IDs)
	attr(res, "region.id") <- as.character(IDs$i_ID)
	class(res) <- c("spatial.neighbour", "data.frame")
	res
}

