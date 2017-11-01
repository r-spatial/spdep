# Copyright 2001-10 by Roger Bivand 
#


is.symmetric.nb <- function(nb, verbose=NULL, force=FALSE)
{
	if(!inherits(nb, "nb")) stop("Not neighbours list")
        if (is.null(verbose)) verbose <- get("verbose", envir = .spdepOptions)
        stopifnot(is.logical(verbose))
	nbsym <- attr(nb, "sym")
	if(!is.null(nbsym)) res <- nbsym
	if(force || is.null(nbsym)) {
		res <- .Call("symtest", nb=nb, card=as.integer(card(nb)),
			verbose=as.logical(verbose), PACKAGE="spdep")
	}
	if(!res && verbose) cat("Non-symmetric neighbours list\n")
	res
}

is.symmetric.glist <- function(nb, glist)
{
	if(!inherits(nb, "nb")) stop("Not neighbours list")
	nbsym <- attr(nb, "sym")
	if(is.null(nbsym)) nbsym <- is.symmetric.nb(nb)
	if (!nbsym) {
		res0 <- vector(mode="list", length=2)
		res0[[1]] <- FALSE
		res0[[2]] <- Inf
	} else {
		if (length(nb) != length(glist)) stop("list lengths differ")
		cnb <- as.integer(card(nb))
		gnb <- as.integer(sapply(glist, length))
		if (!(identical(cnb, gnb))) 
			stop("different vector lengths in lists")
		res0 <- .Call("gsymtest", nb=nb, glist=glist, card=cnb, 
			PACKAGE="spdep")
	}
	res <- res0[[1]]
	attr(res, "d") <- res0[[2]]
	res
}


sym.attr.nb <- function(nb) {
	if(!inherits(nb, "nb")) stop("Not neighbours list")
	nbsym <- attr(nb, "sym")
	if(is.null(nbsym))
		attr(nb, "sym") <- is.symmetric.nb(nb, verbose=FALSE,
			force=TRUE)
	nb
}

include.self <- function(nb) {
	if (!is.null(attributes(nb)$self.included) &&
		(as.logical(attributes(nb)$self.included)))
		stop("Self already included")
	n <- length(nb)
	nc <- card(nb)
	for (i in 1:n) {
		if (nc[i] > 0) {
			nb[[i]] <- sort(c(i, nb[[i]]))
		} else {
			nb[[i]] <- i
		}
	}
		
	attr(nb, "self.included") <- TRUE
	nb
}

is.selfneighbour <- function(nb) {
    res <- sapply(seq(along=nb), function(i) i %in% nb[[i]]) 
}

# Copyright 2001-7 by Nicholas Lewin-Koh and Roger Bivand

old.make.sym.nb <- function(nb){
	if(!inherits(nb, "nb")) stop("Not neighbours list")
	if (is.symmetric.nb(nb, FALSE, TRUE)) {
		res <- nb
	} else {
#        	k <- unlist(lapply(nb,length))
# problems handling no-neighbour entities
		k <- card(nb)
        	to <- unlist(nb)
		to <- to[to > 0]
        	from <- NULL
        	res <- vector(mode="list", length=length(nb))
        	for(i in 1:length(nb)){
        		from <- c(from,rep(i,k[i]))
        	}
        	for(i in 1:length(nb)){
        		res[[i]] <- sort(unique(c(to[from==i],from[to==i])))
        		if(length(res[[i]]) == 0L) res[[i]] <- 0L
        	}
        	attr(res, "region.id") <- attr(nb,"region.id")
        	attr(res, "call") <- attr(nb, "call")
        	attr(res, "type") <- attr(nb, "type")
        	attr(res, "sym") <- TRUE
        	class(res) <- "nb"
	}
	res
}

# Copyright 2009 by Bjarke Christensen and Roger Bivand

make.sym.nb <- function (nb)
{
    if (!inherits(nb, "nb"))
        stop("Not neighbours list")
    if (any(card(nb) == 0)) return(old.make.sym.nb(nb))
    res <- nb
    if (!is.symmetric.nb(nb, FALSE, TRUE)) {
      for (i in 1:length(res)) {
#Which of observation i's neighbors have i amongst _its_ neighbors?
        refersback <- sapply(res[res[[i]]], function(x) i %in% x)
#Add i to the neighborhood of those of i's neighbors who don't refer back
        res[ res[[i]][!refersback] ] <- lapply(res[ res[[i]][!refersback]],
          function(x) sort(c(i, x)))
      }
      attributes(res) <- attributes(res)[!(names(attributes(res))=='knn-k')]
      attr(res, "sym") <- TRUE
    }
    res
}

# Idea due to Roberto Patuelli

aggregate.nb <- function(x, IDs, remove.self=TRUE, ...) {
    stopifnot(length(x) == length(IDs))
    in_reg.ids <- attr(x, "region.id")
    mtch <- tapply(in_reg.ids, IDs, function(i) c(i))
    out_reg.ids <- names(mtch)
    nb_short <- vector(mode="list", length=length(mtch))
    for (i in seq(along=mtch)) {
        nb_short[[i]] <- 0L
        imtch <- match(mtch[[i]], in_reg.ids)
        res <- unlist(x[imtch])
        nb_short[[i]] <- as.integer(sort(unique(match(IDs[res], out_reg.ids))))
        if (remove.self && i %in% nb_short[[i]]) {
            nb_short[[i]] <- nb_short[[i]][-(match(i, nb_short[[i]]))]
            if (length(nb_short[[i]]) < 1L) nb_short[[i]] <- 0L
        }
    }
    attr(nb_short, "region.id") <- out_reg.ids
    class(nb_short) <- "nb"
    nb_short <- sym.attr.nb(nb_short)
    nb_short
}


