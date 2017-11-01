# Copyright 2001 by Roger Bivand
#


knn2nb <- function(knn, row.names=NULL, sym=FALSE) {
	if (class(knn) != "knn") stop("Not a knn object")
	res <- vector(mode="list", length=knn$np)
    	if (!is.null(row.names)) {
		if(length(row.names) != knn$np)
            		stop("row.names wrong length")
		if (length(unique(row.names)) != length(row.names))
	    		stop("non-unique row.names given")
    	}
	if (knn$np < 1) stop("non-positive number of spatial units")
    	if (is.null(row.names)) row.names <- as.character(1:knn$np)
        if(sym){
          to<-as.vector(knn$nn)
          from<-rep(1:knn$np,knn$k)
          for (i in 1:knn$np)res[[i]] <- sort(unique(c(to[from==i],
                                                       from[to==i]))) 
        } else {
          for (i in 1:knn$np) res[[i]] <- sort(knn$nn[i,])
        }
 	attr(res, "region.id") <- row.names
 	attr(res, "call") <- attr(knn, "call")
        attr(res, "sym") <- sym
	attr(res, "type") <- "knn"
 	attr(res, "knn-k") <- knn$k
	class(res) <- "nb"
	res
}
