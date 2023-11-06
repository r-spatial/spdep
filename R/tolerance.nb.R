`tolerance.nb` <-
function (coords, unit.angle = "degrees", max.dist, tolerance, rot.angle, plot.sites=FALSE) {
	coords <- as.matrix(coords)
	if (missing(rot.angle)) {
		rot.angle <- 0
	}else{
		if (rot.angle == 0){
			if (unit.angle == "degrees" | unit.angle == "radians") {
				rot.angle<-0
			}else{
				stop("Unit for angles be either 'degrees' or 'radians'")
			}
		}

		if (rot.angle != 0){
			if (unit.angle == "degrees") {
				rot.angle <- pi/180 * rot.angle
			}else{
				if (unit.angle == "radians") {
					rot.angle <- rot.angle
				}else{
					stop("Unit for angles be either 'degrees' or 'radians'")
				}
			}
		}
	}
	coords <- Rotation(coords, rot.angle)
	if (plot.sites) {
		plot(coords, pch = 19, asp = 1)
	}

	dist.coords <- dist(coords)
	angles <- find.angles(coords)
	if (unit.angle == "degrees") {
		angles <- (angles * 180)/pi
	}

	no.good <- which((angles - tolerance) > 0, arr.ind = TRUE)
	if(nrow(no.good) > 0){
		for (i in 1:nrow(no.good)) {
			angles[no.good[i, 1], no.good[i, 2]] <- NA
		}
	}

	if (missing(max.dist)) {
		max.dist <- max(dist.coords)
	}else{
		too.far <- which(as.matrix(dist.coords) > max.dist, arr.ind = TRUE)
		if(nrow(too.far) > 0){
			for (i in 1:nrow(too.far)) {
				angles[too.far[i, 1], too.far[i, 2]] <- NA
			}
		}
	}

	no.na <- which(!is.na(angles), arr.ind = TRUE)
	if(nrow(no.na) > 0){
		for (i in 1:nrow(no.na)) {
			angles[no.na[i, 1], no.na[i, 2]] <- 1
			angles[no.na[i, 2], no.na[i, 1]] <- 1
		}
	}

	na.all <- which(is.na(angles), arr.ind = TRUE)
	if(nrow(na.all) > 0){
		for (i in 1:nrow(na.all)) {
			angles[na.all[i, 1], na.all[i, 2]] <- 0
		}
	}

	nb.obj <- mat2nb(angles)

	return(nb.obj)
}

mat2nb <- function(x, row.names=NULL) {
    	n <- nrow(x)
	if (n < 1) stop("non-positive number of entities")
	m <- ncol(x)
	if (n != m) stop("x must be a square matrix")
	if (any(x < 0)) stop("values in x cannot be negative")
	if (any(is.na(x))) stop("NA values in x not allowed")
    	if (!is.null(row.names)) {
	    if(length(row.names) != n)
            	stop("row.names wrong length")
	    if (length(unique(row.names)) != length(row.names))
	    	stop("non-unique row.names given")
    	}
    	if (is.null(row.names)) {
	    if (!is.null(row.names(x))) {
		row.names <- row.names(x)
	    } else {
		row.names <- as.character(1:n)
	    }
	}
	neighbours <- vector(mode="list", length=n)
	for (i in 1:n) {
	    nbs  <- which(x[i,] > 0.0)
	    if (length(nbs) > 0) {
		neighbours[[i]] <- nbs
	    } else {
		neighbours[[i]] <- 0L
	    }
        }
	class(neighbours) <- "nb"
	attr(neighbours, "region.id") <- row.names
 	attr(neighbours, "call") <- NA
        attr(neighbours, "sym") <- is.symmetric.nb(neighbours, 
		verbose=FALSE, force=TRUE)
        neighbours
}

`find.angles` <-
function (coords)
{
    n.angles <- (nrow(coords) * (nrow(coords) - 1))/2
    angles <- vector(length = n.angles)
    opp <- matrix(nrow = nrow(coords), ncol = nrow(coords))
    adj <- matrix(nrow = nrow(coords), ncol = nrow(coords))
    for (i in 1:nrow(coords)) {
        for (j in 1:nrow(coords)) {
            opp[i, j] <- coords[j, 2] - coords[i, 2]
            if (opp[i, j] != abs(opp[i, j])) {
                opp[i, j] <- NA
            }
            adj[i, j] <- abs(coords[j, 1] - coords[i, 1])
        }
    }
    opp.adj <- opp/adj
    angles <- pi/2 - atan(opp.adj)
    return(angles)
}
