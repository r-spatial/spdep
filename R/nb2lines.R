# Copyright 2005-7 by Roger Bivand
#

nb2lines <- function(nb, wts, coords, proj4string=CRS(as.character(NA))) {

	x <- coords[,1]
	y <- coords[,2]
	n <- length(nb)
	if (n < 1) stop("zero length neighbour list")
	ID <- as.character(attr(nb, "region.id"))
	cardnb <- card(nb)
	totlinks <- sum(cardnb)
	ll <- vector(mode="list", length=totlinks)
	df <- data.frame(i=integer(totlinks), j=integer(totlinks),
		i_ID=I(character(totlinks)), j_ID=I(character(totlinks)),
		wt=numeric(totlinks))
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
				Ll <- list(Line(xy))
				ll[[line]] <- Lines(Ll, ID=as.character(line))
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
	SpatialLinesDataFrame(SpatialLines(ll, proj4string=proj4string),
		data=df)
#	list(ll=ll, df=df)
}

listw2lines <- function(listw, coords, proj4string=CRS(as.character(NA))) {
	nb2lines(listw$neighbours, listw$weights, coords, proj4string)
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

