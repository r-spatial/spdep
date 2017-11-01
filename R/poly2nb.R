# Copyright 2001-2010 by Roger Bivand 
#
#
# Modified by Micah Altman 2010
	


poly2nb <- function(pl, row.names=NULL, snap=sqrt(.Machine$double.eps),
	queen=TRUE, useC=TRUE, foundInBox=NULL) {
        verbose <- get("verbose", envir = .spdepOptions)
        .ptime_start <- proc.time()
        stopifnot(extends(class(pl), "SpatialPolygons"))

	n <- length(slot(pl, "polygons"))
	if (n < 1) stop("non-positive number of entities")
	if (is.null(row.names)) regid <- row.names(pl)
	else regid <- NULL
	if (is.null(regid)) {
		if(is.null(row.names)) regid <- as.character(1:n)
		else {
			if(length(row.names) != n)
				stop("row.names wrong length")
			else if (length(unique(row.names)) != length(row.names))
	    			stop("non-unique row.names given")
			else regid <- row.names
		}
	}
        if (snap < 0) snap <- abs(snap)
        if (snap < .Machine$double.eps) {
            bsnap <- .Machine$double.eps
        } else { 
            bsnap <- snap
        }
        vbsnap <- c(-bsnap, snap)
        if (verbose)
            cat("handle IDs:", (proc.time() - .ptime_start)[3], "\n")
        .ptime_start <- proc.time()

        xpl <- slot(pl, "polygons")
        xxpl <- vector(mode="list", length=length(xpl))
        for (i in 1:length(xpl)) {
            xpli <- slot(xpl[[i]], "Polygons")
            zz <- lapply(xpli, function(j) slot(j, "coords")[-1,])
            xxpl[[i]] <- do.call("rbind", zz)
        }
        nrs <- sapply(xxpl, nrow)
        bb <- t(sapply(xxpl, function(x) {
            rx <- range(x[,1]) + vbsnap
            ry <- range(x[,2]) + vbsnap
            c(rbind(rx, ry))
        }))
        if (verbose)
            cat("massage polygons:", (proc.time() - .ptime_start)[3], "\n")
        .ptime_start <- proc.time()
#poly2bbs <- function(pl) {
#    t(sapply(slot(pl, "polygons"), bbox))
#}
        genBBIndex<-function(bb) { 
            n <- nrow(bb)
            bxv <- as.vector(bb[,c(1,3)])
            byv <- as.vector(bb[,c(2,4)])
            obxv <- order(bxv)
            rbxv <- c(1:(n*2))[obxv]
            mbxv <- match(1:(n*2),obxv)
            obyv <- order(byv)
            rbyv <- c(1:(n*2))[obyv]
            mbyv <- match(1:(n*2),obyv)

            return(list(bb=bb, bxv=bxv, byv=byv, obxv=obxv, obyv=obyv, 
                mbxv=mbxv, mbyv=mbyv, rbyv=rbyv, rbxv=rbxv))
        }


	dbsnap <- as.double(bsnap)
        dsnap <- as.double(snap)
        if (is.null(foundInBox)) {
	    BBindex <- genBBIndex(bb)
        if (verbose)
            cat("size of BBindex:", object.size(BBindex), "\n")

        } else {
            stopifnot(is.list(foundInBox))
            stopifnot(length(foundInBox) == (n-1L))
            stopifnot(all(unlist(sapply(foundInBox,
                function(x) {if(!is.null(x)) is.integer(x)}))))
            nfIBB <- sum(sapply(foundInBox, length))
        }
        if (verbose)
            cat("generate BBs:", (proc.time() - .ptime_start)[3], "\n")
        .ptime_start <- proc.time()

#	nrs <- integer(n)
#	for (i in 1:n) {
#		pl[[i]] <- na.omit(pl[[i]][-1,])
#		nrs[i] <- as.integer(nrow(pl[[i]]))
#		pl[[i]] <- as.double(pl[[i]])
#	}
	
#        findInBox <- function(i, sp, bigger=TRUE) {
#	    n <- dim(sp$bb)[1]
#	
#	# ! i1 > j3 --> i1 <= j3
#	    tmp1 <- sp$rbxv[sp$mbxv[i]:(n*2)] 
#	    tmp1 <- tmp1[which(tmp1>n)] - n
#	# ! i2 > j4 --> i2 <= bj4
#	    tmp2 <- sp$rbyv[sp$mbyv[i]:(n*2)] 
#	    tmp2 <- tmp2[which(tmp2>n)] - n
#	# ! i3 < j1 -> i3 >= j1
#	    tmp3 <- sp$rbxv[1:sp$mbxv[i+n]] 
#	    tmp3 <- tmp3[which(tmp3<=n)]
#	# ! i4 < j2 -> i4 >= j2
#	    tmp4 <- sp$rbyv[1:sp$mbyv[i+n]] 
#	    tmp4 <- tmp4[which(tmp4<=n)]
#	    result <- intersect(intersect(tmp1,tmp2), intersect(tmp3,tmp4))
#	    if (bigger) {
#		result <- result[which(result>i)]
#	    }
#	    return(sort(result))
#        }

	polypoly2 <- function(poly1, nrs1, poly2, nrs2, snap) {
		if (any(nrs1 == 0 || nrs2 == 0)) return(0L)
		res <- .Call("polypoly", poly1, nrs1, poly2, 
			nrs2, snap, PACKAGE="spdep")
		res
	}

        if (is.null(foundInBox)) {
            foundInBox <- lapply(1:(n-1), function(i) findInBox(i, BBindex))
            if (verbose) {
                cat("findInBox:", (proc.time() - .ptime_start)[3])
            }
          nfIBB <- sum(sapply(foundInBox, length))
          if (verbose) cat(" list size", nfIBB, "\n")

          .ptime_start <- proc.time()
        }
	criterion <- ifelse(queen, 0, 1)
        if (useC) {
#            if (justC) {
              ans <- .Call("poly_loop2", as.integer(n), foundInBox, bb, xxpl,
                as.integer(nrs), as.double(dsnap), as.integer(criterion),
                as.integer(nfIBB), PACKAGE="spdep")
#            } else {
#              ans <- .Call("poly_loop", as.integer(n), i_findInBox, bb, pl, 
#                nrs, as.double(dsnap), as.integer(criterion), as.integer(10),
#                PACKAGE="spdep")
#            }
        } else {
	    ans <- vector(mode="list", length=n)
	    for (i in 1:n) ans[[i]] <- integer(0)
	    for (i in 1:(n-1)) {
		#for (j in (i+1):n) {
#		for (j in findInBox(i,BBindex)) {
		for (j in foundInBox[[i]]) {
			jhit <- .Call("spOverlap", bb[i,], 
				bb[j,], PACKAGE="spdep")
			if (jhit > 0) {
			    khit <- 0
			    khit <- polypoly2(xxpl[[i]], nrs[i], xxpl[[j]], 
				nrs[j], dsnap)

			    if (khit > criterion) {
				ans[[i]] <- c(ans[[i]], j)
				ans[[j]] <- c(ans[[j]], i)
			    }
			}
		}
	    }
	    for (i in 1:n) {
                if (length(ans[[i]]) == 0L) ans[[i]] <- 0L
                if (length(ans[[i]]) > 1L) ans[[i]] <- sort(ans[[i]])
            }
        }
        if (verbose)
            cat("work loop:", (proc.time() - .ptime_start)[3], "\n")
        .ptime_start <- proc.time()
	class(ans) <- "nb"
	attr(ans, "region.id") <- regid
	attr(ans, "call") <- match.call()
	if (queen) attr(ans, "type") <- "queen"
	else attr(ans, "type") <- "rook"
	ans <- sym.attr.nb(ans)
        if (verbose)
            cat("done:", (proc.time() - .ptime_start)[3], "\n")
        .ptime_start <- proc.time()
	ans
}	


# faster findInBox

qintersect<-function(x,y) {
	    # streamlined intersect function for unique vectors
    as.integer(y[match(x, y, 0L)])
}

findInBox<-function(i, sp, bigger=TRUE) {
    n <- dim(sp$bb)[1]

# use index structure to identify which other BB's fall in i's BB
# by getting id's of polygons with BBmin_j < BBmax_i, BBmax_j > BBmin_i for x and y 
# then taking the intersection of these four lists of id's

    tmp<-vector(mode="list", length=4)
        # ! i1 > j3 --> i1 <= j3
    tmp[[1]] <- sp$rbxv[sp$mbxv[i]:(n*2)]
    tmp[[1]]<- tmp[[1]][which(tmp[[1]]>n)] - n
        # ! i2 > j4 --> i2 <= bj4
    tmp[[2]] <- sp$rbyv[sp$mbyv[i]:(n*2)]
    tmp[[2]]<- tmp[[2]][which(tmp[[2]]>n)] - n
        # ! i3 < j1 -> i3 >= j1
    tmp[[3]] <- sp$rbxv[1:sp$mbxv[i+n]]
    tmp[[3]] <- tmp[[3]][which(tmp[[3]]<=n)]
        # ! i4 < j2 -> i4 >= j2
    tmp[[4]] <- sp$rbyv[1:sp$mbyv[i+n]]
    tmp[[4]]<- tmp[[4]][which(tmp[[4]]<=n)]

	# for performance, order the comparison of the lists

    lentmp <- order(sapply(tmp,length))

	# use qintersect, since these are already vectors and unique 
    result <- qintersect(tmp[[lentmp[2]]],tmp[[lentmp[1]]])
    result <- qintersect(tmp[[lentmp[3]]],result)
    result <- qintersect(tmp[[lentmp[4]]],result)

    if (bigger) {
        result<-result[which(result>i)]
    }
    return(sort(result))
}


  
