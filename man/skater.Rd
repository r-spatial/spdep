\name{skater}
\alias{skater}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Spatial 'K'luster Analysis by Tree Edge Removal}
\description{
  This function implements a SKATER procedure for spatial clustering
  analysis. This procedure essentialy begins with an edges set, a data set and
  a number of cuts. The output is an object of 'skater' class and is
  valid for input again.
}
\usage{
skater(edges, data, ncuts, crit, vec.crit, method = c("euclidean", 
    "maximum", "manhattan", "canberra", "binary", "minkowski", 
    "mahalanobis"), p = 2, cov, inverted = FALSE) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{edges}{A matrix with 2 colums with each row is an edge}
  \item{data}{A data.frame with data observed over nodes.}
  \item{ncuts}{The number of cuts}
  \item{crit}{A scalar or two dimensional vector with
    criteria for groups. Examples: limits of group size or
    limits of population size. If scalar, is the minimum
    criteria for groups.}
  \item{vec.crit}{A vector for evaluating criteria.}
  \item{method}{Character or function to declare distance method.
    If \code{method} is character, method must be "mahalanobis" or
    "euclidean", "maximum", "manhattan", "canberra", "binary"
    or "minkowisk".
    If \code{method} is one of "euclidean", "maximum",
    "manhattan", "canberra", "binary" or "minkowski", see
    \code{\link{dist}} for details,
     because this function as used to compute the distance.
     If \code{method="mahalanobis"}, the mahalanobis distance
     is computed between neighbour areas.
     If \code{method} is a \code{function}, this function is
     used to compute the distance.}
  \item{p}{The power of the Minkowski distance.}
  \item{cov}{The covariance matrix used to compute the mahalanobis
   distance.}
  \item{inverted}{logical.  If 'TRUE', 'cov' is supposed to contain the
          inverse of the covariance matrix.}
}

\value{
  A object of \code{skater} class with:
  
  \item{groups}{A vector with length equal the number of nodes. Each
    position identifies the group of node} 
  \item{edges.groups}{A list of length equal the number of groups with
    each element is a set of edges} 
  \item{not.prune}{A vector identifying the groups with are not
    candidates to partition.}
  \item{candidates}{A vector identifying the groups with are candidates
    to partition.} 
  \item{ssto}{The total dissimilarity in each step of edge removal.}

}
\references{
  Assuncao, R.M., Lage J.P., and Reis, E.A. (2002). Analise de
  conglomerados espaciais via arvore geradora minima. Revista Brasileira
  de Estatistica, 62, 1-23. 

  Assuncao, R. M, Neves, M. C., Camara, G. and Freitas, C. da
  C. (2006). Efficient regionalization techniques for socio-economic
  geographical units using minimum spanning trees. International Journal
  of Geographical Information Science Vol. 20, No. 7, August 2006,
  797-811 
}
\author{Renato M. Assuncao and Elias T. Krainski}
%%\note{}
\seealso{See Also as \code{\link{mstree}}}
\examples{
### loading data
GDAL37 <- numeric_version(unname(sf::sf_extSoftVersion()["GDAL"]), strict=FALSE)
(GDAL37 <- ifelse(is.na(GDAL37), FALSE, GDAL37 >= "3.7.0"))
file <- "etc/shapes/bhicv.gpkg.zip"
zipfile <- system.file(file, package="spdep")
if (GDAL37) {
    bh <- st_read(zipfile)
} else {
    td <- tempdir()
    bn <- sub(".zip", "", basename(file), fixed=TRUE)
    target <- unzip(zipfile, files=bn, exdir=td)
    bh <- st_read(target)
}
### data standardized 
dim(bh)
dpad <- data.frame(scale(as.data.frame(bh)[,5:8]))

### neighboorhod list
bh.nb <- poly2nb(bh)
bh.nb

### calculating costs
lcosts <- nbcosts(bh.nb, dpad)
head(lcosts)

### making listw
nb.w <- nb2listw(bh.nb, lcosts, style="B")
nb.w

### find a minimum spanning tree
mst.bh <- mstree(nb.w,5)
str(mst.bh)

### the mstree plot
par(mar=c(0,0,0,0))
plot(st_geometry(bh), border=gray(.5))
pts <- st_coordinates(st_centroid(bh))
plot(mst.bh, pts, col=2, 
     cex.lab=.6, cex.circles=0.035, fg="blue", add=TRUE)

### three groups with no restriction
res1 <- spdep::skater(edges=mst.bh[,1:2], data=dpad, ncuts=2)

### groups size
table(res1$groups)

### the skater plot
opar <- par(mar=c(0,0,0,0))
plot(res1, pts, cex.circles=0.035, cex.lab=.7)

### the skater plot, using other colors
plot(res1, pts, cex.circles=0.035, cex.lab=.7,
     groups.colors=heat.colors(length(res1$ed)))

### the Spatial Polygons plot
plot(st_geometry(bh), col=heat.colors(length(res1$edg))[res1$groups])

#par(opar)
### EXPERT OPTIONS

### more one partition
res1b <- spdep::skater(res1, dpad, 1)

### length groups frequency
table(res1$groups)

table(res1b$groups)

### thee groups with minimum population 
res2 <- spdep::skater(mst.bh[,1:2], dpad, 2, 200000, bh$Pop)
table(res2$groups)

### thee groups with minimun number of areas
res3 <- spdep::skater(mst.bh[,1:2], dpad, 2, 3, rep(1,nrow(bh)))
table(res3$groups)

### thee groups with minimun and maximun number of areas
res4 <- spdep::skater(mst.bh[,1:2], dpad, 2, c(20,50), rep(1,nrow(bh)))
table(res4$groups)

### if I want to get groups with 20 to 40 elements
res5 <- spdep::skater(mst.bh[,1:2], dpad, 2,
   c(20,40), rep(1,nrow(bh))) ## DON'T MAKE DIVISIONS 
table(res5$groups)

### In this MST don't have groups with this restrictions
### In this case, first I do one division
### with the minimun criteria
res5a <- spdep::skater(mst.bh[,1:2], dpad, 1, 20, rep(1,nrow(bh))) 
table(res5a$groups)

### and do more one division with the full criteria
res5b <- spdep::skater(res5a, dpad, 1, c(20, 40), rep(1,nrow(bh)))
table(res5b$groups)

### and do more one division with the full criteria
res5c <- spdep::skater(res5b, dpad, 1, c(20, 40), rep(1,nrow(bh)))
table(res5c$groups)

### It don't have another divison with this criteria
res5d <- spdep::skater(res5c, dpad, 1, c(20, 40), rep(1,nrow(bh)))
table(res5d$groups)

\dontrun{
data(boston, package="spData")
bh.nb <- boston.soi
dpad <- data.frame(scale(boston.c[,c(7:10)]))
### calculating costs
system.time(lcosts <- nbcosts(bh.nb, dpad))
### making listw
nb.w <- nb2listw(bh.nb, lcosts, style="B")
### find a minimum spanning tree
mst.bh <- mstree(nb.w,5)
### three groups with no restriction
system.time(res1 <- spdep::skater(mst.bh[,1:2], dpad, 2))
library(parallel)
nc <- max(2L, detectCores(logical=FALSE), na.rm = TRUE)-1L
# set nc to 1L here
if (nc > 1L) nc <- 1L
coresOpt <- get.coresOption()
invisible(set.coresOption(nc))
if(!get.mcOption()) {
# no-op, "snow" parallel calculation not available
  cl <- makeCluster(get.coresOption())
  set.ClusterOption(cl)
}
### calculating costs
system.time(plcosts <- nbcosts(bh.nb, dpad))
all.equal(lcosts, plcosts, check.attributes=FALSE)
### making listw
pnb.w <- nb2listw(bh.nb, plcosts, style="B")
### find a minimum spanning tree
pmst.bh <- mstree(pnb.w,5)
### three groups with no restriction
system.time(pres1 <- spdep::skater(pmst.bh[,1:2], dpad, 2))
if(!get.mcOption()) {
  set.ClusterOption(NULL)
  stopCluster(cl)
}
all.equal(res1, pres1, check.attributes=FALSE)
invisible(set.coresOption(coresOpt))
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{cluster}
\keyword{tree}% __ONLY ONE__ keyword per line
