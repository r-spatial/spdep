# Copyright 2003-6 by Luc Anselin and Roger Bivand
#

# LA 6/28/03 read.gwt
# LA 7/12/03 revised, sorted ids
# LA 9/30/03 use match to correct orders

read.gwt2nb <- function(file, region.id=NULL) {
	con <- file(file, open="r")   #opens the file
	firstline <- unlist(strsplit(readLines(con,1)," "))
	if (length(firstline) == 4L) {
		n <- as.integer(firstline[2])
		shpfile <- firstline[3]
		ind <- firstline[4]
		if (ind != deparse(substitute(region.id)))
			warning(paste("region.id not named", ind))
	} else if (length(firstline) == 1L) {
		n <- as.integer(firstline[1])
		shpfile <- as.character(NA)
		ind <- as.character(NA)
		warning("Old-style GWT file")
	} else stop("Invalid header line format for GWT file")
	close(con)
	if (n < 1) stop("non-positive number of entities")
	nseq <- 1:n
	if (is.null(region.id)) region.id <- nseq
	if (n != length(region.id))
		stop("Mismatch in dimensions of GWT file and region.id")
	if (length(unique(region.id)) != length(region.id))
	    stop("non-unique region.id given")
	odij <- read.table(file, skip=1)
	    # convert region.id to order
	regodij <- match(odij[,1], region.id)
	regddij <- match(odij[,2], region.id)
	odij <- cbind(regodij, regddij, odij[,3])
	qorder <- order(odij[,1],odij[,2])
	odij <- odij[qorder,]
	origvec <- unique(odij[,1])
	if (!all(nseq %in% origvec))
		warning(paste(paste(region.id[which(!(nseq %in% origvec))], 
			collapse=", "), "are not origins"))
	destvec <- unique(odij[,2])
	if (!all(nseq %in% destvec))
		warning(paste(paste(region.id[which(!(nseq %in% destvec))], 
			collapse=", "), "are not destinations"))

	res <- vector(mode="list", length=n)
	vlist <- vector(mode="list", length=n)
	rle.sn <- rle(odij[,1])
	cs1.sn <- cumsum(rle.sn$lengths)
	cs0.sn <- c(1, cs1.sn[1:(n-1)]+1)
	ii <- 1
	for (i in 1:n) {
# Bug hit by Thomas Halvorsen 10/2006, was already fixed in sn2listw()
		if (!is.na(rle.sn$value[ii]) && rle.sn$value[ii] == i) {
			res[[i]] <- as.integer(odij[cs0.sn[ii]:cs1.sn[ii],2])
			vlist[[i]] <- as.double(odij[cs0.sn[ii]:cs1.sn[ii],3])
			ii <- ii+1
		} else {
			res[[i]] <- 0L
		}
	}

	class(res) <- c("nb", "GWT")
	attr(res, "region.id") <- region.id
	attr(res, "neighbours.attrs") <- as.character(NA)
	attr(res, "weights.attrs") <- as.character(NA)
	attr(res, "GeoDa") <- list(dist=vlist, shpfile=shpfile, ind=ind)
	attr(res, "call") <- match.call()
	attr(res, "n") <- n
	res <- sym.attr.nb(res)
	res
}

write.sn2gwt <- function(sn, file, shpfile=NULL, ind=NULL, useInd=FALSE, legacy=FALSE) {
	if(!inherits(sn, "spatial.neighbour")) 
	    stop("not a spatial.neighbour object")
	n <- attr(sn, "n")
	if (n < 1) stop("non-positive number of entities")
	if (is.null(shpfile)) {
		tmp <- attr(sn, "GeoDa")$shpfile
		if (is.null(tmp)) shpfile <- "unknown"
		else shpfile <- tmp
	} else {
            stopifnot(is.character(shpfile))
            stopifnot(length(shpfile) == 1L)
        }
	if (is.null(ind)) {
		tmp <- attr(sn, "GeoDa")$ind
		if (is.null(tmp)) ind <- "unknown"
		else ind <- tmp
	} else {
            stopifnot(is.character(ind))
            stopifnot(length(ind) == 1L)
        }
        if (useInd) {
            rid <- attr(sn, "region.id")
            sn$from <- rid[sn$from]
            sn$to <- rid[sn$to]
        }
	con <- file(file, open="w")
	if (legacy) writeLines(format(n), con)
        else writeLines(paste("0", n, shpfile, ind, sep=" "), con)
	write.table(as.data.frame(sn), file=con, append=TRUE,
		row.names=FALSE, col.names=FALSE, quote=FALSE)
	close(con)
}

write.sn2dat <- function(sn, file) {
	if(!inherits(sn, "spatial.neighbour")) 
	    stop("not a spatial.neighbour object")
	write.table(data.frame(sn[order(sn[,2]), ]), 
	    file=file, col.names=FALSE, row.names=FALSE)
}

read.dat2listw <- function(file) {
	wmat <- read.table(file)
        stopifnot(ncol(wmat) == 3)
        stopifnot(is.numeric(wmat[,3]))
        if (storage.mode(wmat[,1]) != "integer")
            storage.mode(wmat[,1])<- "integer"
        if (storage.mode(wmat[,2]) != "integer")
            storage.mode(wmat[,2]) <- "integer"
	sn <- wmat[order(wmat[,1]),]
	IDS <- unique(sn[,1])
	class(sn) <- c("spatial.neighbour", "data.frame")
	attr(sn, "n") <- length(IDS)
	attr(sn, "region.id") <- as.character(IDS)
	listw <- sn2listw(sn)
	listw
}


write.sn2Arc <- function(sn, file, field=NULL) {
	if(!inherits(sn, "spatial.neighbour")) 
	    stop("not a spatial.neighbour object")
	if (is.null(field)) stop("field must be given")
	n <- attr(sn, "n")
	if (n < 1) stop("non-positive number of entities")
	nms <- as.character(attr(sn, "region.id"))
	sn[,1] <- nms[sn[,1]]
	sn[,2] <- nms[sn[,2]]
	con <- file(file, open="w")
	writeLines(field, con)
	write.table(as.data.frame(sn), file=con, append=TRUE,
		row.names=FALSE, col.names=FALSE, quote=FALSE)
	close(con)
}

# Copyright 2011 Virgilio Gomez-Rubio
# a function to export from nb object to a particular file format
# which is used by INLA when fitting spatial models for lattice data
# revision Marcos Prates 2011-10-07

nb2INLA <-function(file, nb)
{
        
        n<-length(nb)

        if(!file.create(file))
        {
                stop("Cannot open file")
        }

        txt<-paste(n, "\n", sep="")
        cat(txt, file=file, append = TRUE)
        crd <- card(nb)

        for(i in 1:length(nb))
        {
                if (crd[i] == 0) txt <- paste(c(i, 0), collapse=" ")
# Marcos Prates 2011-10-07 no-neighbour case
                else txt<-paste(c(i, length(nb[[i]]), nb[[i]]), collapse=" ")
                txt<-paste(txt, "\n",sep="")
                cat(txt, file=file, append = TRUE)
        }
}

