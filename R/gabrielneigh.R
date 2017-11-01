# Copyright 2001 by Nicholas Lewin-Koh 
#

gabrielneigh <- function(coords, nnmult=3) {
    x <- coords
    if (!is.matrix(x)) stop("Data not in matrix form")
    if (any(is.na(x))) stop("Data cannot include NAs")
    np <- nrow(x)
    if(ncol(x)!=2) stop("Planar graphs only work in 2d")
    ngaballoc <- np*nnmult
    g1<-g2<-rep(0,ngaballoc)
    nogab <- 0
    storage.mode(x) <- "double"
    z <- .C("compute_gabriel", np=as.integer(np), from=as.integer(g1),
             to=as.integer(g2), nedges=as.integer(nogab), 
             ngaballoc=as.integer(ngaballoc), x=x[,1], 
             y=x[,2], PACKAGE="spdep")
    z$from<-z$from[1:z$nedges]
    z$to<-z$to[1:z$nedges]
    attr(z, "call") <- match.call()
    class(z)<-c("Graph","Gabriel")
    z
}

plot.Gabriel<-function(x, show.points=FALSE, add=FALSE,
                       linecol=par(col), ...)
{
  if(!add) plot(x$x,x$y,type='n',...)
  segments(x$x[x$from], x$y[x$from],
           x$x[x$to], x$y[x$to], col=linecol)
  if(show.points) points(x$x,x$y)
}

