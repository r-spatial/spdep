`mstree` <-
function(nbw, ini=NULL) {
  n <- length(nbw[[2]])
  nodes <- cbind(FALSE, 0, rep(Inf,n))
  if (is.null(ini))
    ini <- sample(1:n, 1)
  nodes[ini, 1] <- TRUE
  nodes[nbw$neighbours[[ini]], 2] <- ini
  nodes[nbw$neighbours[[ini]], 3] <- nbw$weights[[ini]]
  
  mst <- matrix(0, n-1, 3)
  for (i in 1:(n-1)){
    id.min <- which.min(nodes[,3])
    if (!is.finite(nodes[id.min,3]))
      stop("Graph is not connected!")
    nodes[id.min, 1] <- TRUE
    mst[i, ] <- c(nodes[id.min, 2], id.min, nodes[id.min, 3])
    id.out <- !nodes[nbw$neighbours[[id.min]], 1]
    node.can <- nbw$neighbours[[id.min]][id.out]
    node.cost <- nbw$weights[[id.min]][id.out]
    id.best <- node.cost<nodes[node.can,3]
    nodes[node.can[id.best], 2] <- id.min
    nodes[node.can[id.best], 3] <- node.cost[id.best]
    nodes[id.min, 3] <- Inf
  }
  attr(mst, "class") <- c("mst", "matrix")
  mst
}

