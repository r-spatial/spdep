# Copyright 2001 by Nicholas Lewin-Koh, igraph added RSB 2024
#


n.comp.nb <- function(nb.obj){
  if(!inherits(nb.obj,"nb")) stop("not a neighbours list")
  if (sum(card(nb.obj)) == 0L) {
    return(list(nc=length(nb.obj), comp.id=1:length(nb.obj)))
  }
  nb.sym <- is.symmetric.nb(nb.obj)
  igraph <- FALSE
  if (requireNamespace("igraph", quietly=TRUE) &&
      requireNamespace("spatialreg", quietly=TRUE)) {
      igraph <- TRUE
  }
  if (!igraph) {
    if (!nb.sym) nb.obj <- make.sym.nb(nb.obj)
    comp0 <- rep(0,length(nb.obj))
    comp <- .Call("g_components", nb.obj, as.integer(comp0), PACKAGE="spdep")
    answ <- list(nc=length(unique(comp)), comp.id=comp)
  } else {
    stopifnot(requireNamespace("igraph", quietly=TRUE))
    stopifnot(requireNamespace("spatialreg", quietly=TRUE))
    B <- as(nb2listw(nb.obj, style="B", zero.policy=TRUE), "CsparseMatrix")
    g1 <- igraph::graph_from_adjacency_matrix(B,
      mode=ifelse(nb.sym, "undirected", "directed"))
    c1 <- igraph::components(g1, mode="weak")
    answ <- list(nc=c1$no, comp.id=unname(c1$membership))
  }
  answ
}

