# Copyright 2001-2006 by Nicholas Lewin-Koh and Roger Bivand
#


union.nb<-function(nb.obj1, nb.obj2){
  if(!inherits(nb.obj1,"nb") | !inherits(nb.obj2,"nb")){
    stop("Both arguments must be of class nb")
  }
  if(any(attr(nb.obj1,"region.id")!= attr(nb.obj2,"region.id"))){
   stop("Both neighbor objects must be \n generated from the same coordinates")
  }
  n <- length(nb.obj1)
  if (n != length(nb.obj2)) stop("Both arguments must be of same length")
  if (n < 1) stop("non-positive number of entities")
  card1 <- card(nb.obj1)
  card2 <- card(nb.obj2)
  new.nb<-vector(mode="list", length=n)
  for(i in 1:n) {
    if (card1[i] == 0) {
      if (card2[i] == 0) new.nb[[i]] <- 0L
      else new.nb[[i]] <- nb.obj2[[i]]
    } else {
      if (card2[i] == 0) new.nb[[i]] <- nb.obj1[[i]]
      else new.nb[[i]]<-sort(union(nb.obj1[[i]], nb.obj2[[i]]))
    }
  }
  attr(new.nb,"region.id")<-attr(nb.obj1,"region.id")
  attr(new.nb,"type")<-paste("union(",attr(nb.obj1,"type"),
                             ",",attr(nb.obj2,"type"),")")
  class(new.nb)<-"nb"
  new.nb
 }

intersect.nb<-function(nb.obj1, nb.obj2){
  if(!inherits(nb.obj1,"nb") | !inherits(nb.obj2,"nb")){
    stop("Both arguments must be of class nb")
  }
  if(any(attr(nb.obj1,"region.id")!= attr(nb.obj2,"region.id"))){
   stop("Both neighbor objects must be \n generated from the same coordinates")
  }
  n <- length(nb.obj1)
  if (n != length(nb.obj2)) stop("Both arguments must be of same length")
  if (n < 1) stop("non-positive number of entities")
  card1 <- card(nb.obj1)
  card2 <- card(nb.obj2)
  new.nb<-vector(mode="list", length=n)
  for(i in 1:n) {
    if (card1[i] > 0 && card2[i] > 0) {
      res <- sort(intersect(nb.obj1[[i]], nb.obj2[[i]]))
      if(length(res) == 0L) new.nb[[i]] <- 0L
      else new.nb[[i]] <- res
    } else new.nb[[i]] <- 0L
  }
  attr(new.nb,"region.id")<-attr(nb.obj1,"region.id")
  attr(new.nb,"type")<-paste("intersect(",attr(nb.obj1,"type"),
                             ",",attr(nb.obj2,"type"),")")
  class(new.nb)<-"nb"
  new.nb
}
setdiff.nb<-function(nb.obj1, nb.obj2){
  	if(!inherits(nb.obj1,"nb") | !inherits(nb.obj2,"nb")){
    		stop("Both arguments must be of class nb")
  	}
  	if(any(attr(nb.obj1,"region.id")!= attr(nb.obj2,"region.id"))){
   		stop("Both neighbor objects must be \n generated from the same coordinates")
  	}
  	n <- length(nb.obj1)
  	if (n != length(nb.obj2)) stop("Both arguments must be of same length")
	if (n < 1) stop("non-positive number of entities")
  	card1 <- card(nb.obj1)
  	card2 <- card(nb.obj2)
  	new.nb<-vector(mode="list", length=n)
  	for(i in 1:n) {
    		if (card1[i] == 0) {
      			if (card2[i] == 0) new.nb[[i]] <- 0L
      			else new.nb[[i]] <- nb.obj2[[i]]
    		} else {
      			if (card2[i] == 0) new.nb[[i]] <- nb.obj1[[i]]
      			else {
            			if (card2[i] == 0)
                			new.nb[[i]] <- nb.obj1[[i]]
            			else {
                			if (card1[i] >= card2[i]) {
                  				a <- nb.obj1[[i]]
                  				b <- nb.obj2[[i]]
                			} else {
                  				b <- nb.obj1[[i]]
                  				a <- nb.obj2[[i]]
                			}
                			res <- sort(setdiff(a, b))
			                if(length(res) == 0L) 
					    new.nb[[i]] <- 0L
                			else new.nb[[i]] <- res
	    			}
        		}
    		}
  	}
  	attr(new.nb,"region.id")<-attr(nb.obj1,"region.id")
  	attr(new.nb,"type")<-paste("setdiff(",attr(nb.obj1,"type"),
                             ",",attr(nb.obj2,"type"),")")
  	class(new.nb)<-"nb"
  	new.nb
}

complement.nb<-function(nb.obj){
   if(!inherits(nb.obj,"nb")){
    stop("Argument must be of class nb")
   }
  n <- length(nb.obj)
  if (n < 1) stop("non-positive number of entities")
  card1 <- card(nb.obj)
  new.nb<-vector(mode="list", length=n)
  cmp<-1:n
  attributes(new.nb)<-attributes(nb.obj)
  for(i in 1:n) {
    if (card1[i] == 0) new.nb[[i]] <- cmp
    else {
      res <- sort(cmp[-nb.obj[[i]]])
      if(length(res) == 0L) new.nb[[i]] <- 0L
      else new.nb[[i]] <- res
    }
  }
  attr(new.nb,"type")<-paste("complement(",attr(nb.obj,"type"),")")
  class(new.nb)<-"nb"
  new.nb
 }
