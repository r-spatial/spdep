`plot.mst` <-
function(x, coords, label.areas=NULL,
                     cex.circles=1, cex.labels=1, ...){
###  funcao para plotar o grafo da arvore geradora minima
###  vec.argem e' o vetor com os indices do vizinho de conexao
###       de cada area na arvore geradora minima
###  
##
   n <- nrow(coords)
   if (is.null(label.areas))
     label.areas <- as.character(1:n) 
   symbols(coords[,1], coords[,2], circles=rep(cex.circles,n),
           inches=FALSE, xlab=" ", ylab=" ", xaxt="n", yaxt="n", ...)
   text(coords[,1], coords[,2], label.areas, cex=cex.labels)
   segments(coords[x[,1],1], coords[x[,1],2],
            coords[x[,2],1], coords[x[,2],2], ...)
   invisible()
}

