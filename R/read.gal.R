# Copyright 2001-2003 by Roger Bivand and Luc Anselin
#

read.gal <- function(file, region.id=NULL, override.id=FALSE) 
{
	con <- file(file, open="r")
	line <- unlist(strsplit(readLines(con, 1), " "))
	x <- subset(line, nchar(line) > 0)
	if (length(x) == 1L) {
		n <- as.integer(x[1])
		shpfile <- as.character(NA)
		ind <- as.character(NA)
	} else if (length(x) == 4L) {
		n <- as.integer(x[2])
		shpfile <- as.character(x[3])
		ind <- as.character(x[4])
	} else stop ("Invalid header line in GAL file")
	if (n < 1) stop("Non-positive number of regions")
	if (!is.null(region.id))
		if (length(unique(region.id)) != length(region.id))
	    		stop("non-unique region.id given")
	if (is.null(region.id)) region.id <- as.character(1:n)
	if (n != length(region.id))
		stop("Mismatch in dimensions of GAL file and region.id")
    	rn <- character(n)
	res <- vector(mode="list", length=n)
	for (i in 1:n) {
		line <- unlist(strsplit(readLines(con, 1), " "))
		x <- subset(line, nchar(line) > 0)
		rn[i] <- x[1]
		line <- unlist(strsplit(readLines(con, 1), " "))
		y <- subset(line, nchar(line) > 0)
		if(length(y) != as.integer(x[2])) {
			close(con)
			stop(paste("GAL file corrupted at region", i))
		}
		res[[i]] <- y
	}
	close(con)
	if (!override.id) {
		if (!all(rn %in% region.id)) {
			stop("GAL file IDs and region.id differ")
		}
	} else region.id <- rn
	mrn <- match(rn, region.id)
	res1 <- vector(mode="list", length=n)
	for (i in 1:n) {
            if (length(res[[i]]) > 0) {    
		x <- match(res[[i]], region.id)
		if (any(is.na(x)) | (length(x) != length(res[[i]]))) {
			stop(paste("GAL file corrupted at region", i))
		}
		if(any(x < 0) || any(x > n))
			stop("GAL file corrupted")

		res1[[mrn[i]]] <- sort(x)
            } else {
                res1[[mrn[i]]] <- 0L
            }
	}
	class(res1) <- "nb"
    	attr(res1, "region.id") <- region.id
	attr(res1, "GeoDa") <- list(shpfile=shpfile, ind=ind)
	attr(res1, "gal") <- TRUE
	attr(res1, "call") <- TRUE
	res1 <- sym.attr.nb(res1)
	res1
}

write.nb.gal <- function(nb, file, oldstyle=TRUE, shpfile=NULL, ind=NULL) {
# class to inherits Jari Oksanen 080603
  	if (!inherits(nb, "nb")) stop("not a neighbours list")
	n <- length(nb)
	if (n < 1) stop("non-positive number of entities")
	cn <- card(nb)
	rn <- attr(nb, "region.id")
	if (is.null(shpfile)) {
		tmp <- attr(nb, "GeoDa")$shpfile
		if (is.null(tmp)) shpfile <- "unknown"
		else shpfile <- tmp
	}
	if (is.null(ind)) {
		tmp <- attr(nb, "GeoDa")$ind
		if (is.null(tmp)) ind <- "unknown"
		else ind <- tmp
	}
	con <- file(file, open="w")
	if (oldstyle) writeLines(paste(n), con)
	else writeLines(paste("0", n, shpfile, ind, sep=" "), con)
	for (i in 1:n) {
		if (oldstyle) 
			writeLines(paste(i, cn[i], 
				collapse=" "), con)
		else writeLines(paste(rn[i], cn[i],
			collapse=" "), con)
		if (oldstyle) writeLines(ifelse(cn[i] > 0, 
			paste(nb[[i]], collapse=" "), ""), con)
		else writeLines(ifelse(cn[i] > 0,
			paste(rn[nb[[i]]], collapse=" "), ""), con)
	}
	close(con)
}


# read.geoda
# helper function to read GeoDa export files
# LA 7/10/03
# specify input file = file
# default is no row names, specify row names as second parameter
# example: balt <- read.geoda("baltim.txt")
#          balt <- read.geoda("baltim.txt","STATION")

read.geoda <- function(file, row.names=NULL, skip=0)
{
	res <- read.csv(file=file, header=TRUE, skip=skip, row.names=row.names)
	if (ncol(res) < 2) warning("data frame possibly malformed") 
	res
}

