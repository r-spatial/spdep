# Copyright 2001 by Nicholas Lewin-Koh, igraph added RSB 2024
#


n.comp.nb <- function(nb.obj, igraph=FALSE){
  if(!inherits(nb.obj,"nb"))stop("not a neighbours list")
  stopifnot(is.logical(igraph))
  stopifnot(length(igraph) == 1L)
  nb.sym <- is.symmetric.nb(nb.obj)
  if (igraph) {
    if (!requireNamespace("igraph", quietly=TRUE)) {
      igraph <- !igraph
      warning("igraph not available, set FALSE")
    }
  }
  if (!igraph) {
    if (!nb.sym) nb.obj <- make.sym.nb(nb.obj)
    comp <- rep(0,length(nb.obj))
    comp <- .Call("g_components", nb.obj, as.integer(comp), PACKAGE="spdep")
    answ <- list(nc=length(unique(comp)), comp.id=comp)
  } else {
    stopifnot(requireNamespace("igraph", quietly=TRUE))
    stopifnot(requireNamespace("spatialreg", quietly=TRUE))
    B <- as(nb2listw(nb.obj, style="B", zero.policy=TRUE), "CsparseMatrix")
    
    g1 <- igraph::graph_from_adjacency_matrix(B,
      mode=ifelse(nb.sym, "undirected", "directed"))
    c1 <- igraph::components(g1, mode="weak")
    answ <- list(nc=c1$no, comp.id=c1$membership)
  }
  answ
}

