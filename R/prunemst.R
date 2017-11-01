prunemst <- function(edges, only.nodes=TRUE) {
  id1 <- .C("prunemst", as.integer(edges[,1]),
            as.integer(edges[,2]),
            as.integer(nrow(edges)),
            integer(nrow(edges)), PACKAGE="spdep")[[4]]
  no1 <- unique(c(edges[1,1], as.integer(edges[id1==1,])))
  no2 <- setdiff(edges, no1)
  if(only.nodes)
    return(list(node1=no1, node2=no2))
  else {
    edges <- edges[-1, , drop=FALSE]
    return(list(list(node=no1,
                     edge=edges[id1[-1]==1,, drop=FALSE]), 
                list(node=no2,
                     edge=edges[id1[-1]==0,, drop=FALSE])))
  }
}    

