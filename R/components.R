# Copyright 2001 by Nicholas Lewin-Koh 
#


n.comp.nb <- function(nb.obj){
  if(!inherits(nb.obj,"nb"))stop("not a neighbours list")
  nb.obj <- make.sym.nb(nb.obj)
  comp <- rep(0,length(nb.obj))
  comp <- .Call("g_components", nb.obj, as.integer(comp), PACKAGE="spdep")
  answ <- list(nc=length(unique(comp)), comp.id=comp)
  answ
}

