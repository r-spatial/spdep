# Copyright 2001-6 by Nicholas Lewin-Koh and Roger S. Bivand.
#


graph2nb <- function(gob, row.names=NULL,sym=FALSE) {
	if (!inherits(gob, "Graph")) stop("Not a Graph object")
	res <- vector(mode="list", length=gob$np)
    	if (!is.null(row.names)) {
		if(length(row.names) != gob$np)
            		stop("row.names wrong length")
		if (length(unique(row.names)) != length(row.names))
	    		stop("non-unique row.names given")
    	}
	if (gob$np < 1) stop("non-positive gob$np")
    	if (is.null(row.names)) row.names <- as.character(1:gob$np)
        if(sym){
          for (i in 1:gob$np) {
		res[[i]] <- sort(unique(c(gob$to[gob$from==i],
                                       gob$from[gob$to==i])))
	  	if(length(res[[i]]) == 0L) res[[i]] <- 0L
	  }
        }
        else{
	  for (i in 1:gob$np) {
		res[[i]] <- sort(gob$to[gob$from==i])
	  	if(length(res[[i]]) == 0L) res[[i]] <- 0L
	  }
        }
        attr(res, "region.id") <- row.names
 	attr(res, "call") <- attr(gob, "call")
 	attr(res, "type") <- attr(gob, "type")
	class(res) <- "nb"
	res <- sym.attr.nb(res)
	res
}
