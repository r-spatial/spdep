grid2nb <- function(grid, d = grid@cells.dim,
                    queen = TRUE, nb = TRUE, self = FALSE) {
### supply either a grid (GridTopology)
###    or d (number of cells in each dimention) 
### choose queen and nb (see the poly2nb() function)
    n <- as.integer(prod(d))
    if (n==1) {
        if (nb) {
            nn <- list(self*1L)
            class(nn) <- 'nb'
            return(nn)
        } else return(matrix(0L, 1, 1))
    }
    nd <- as.integer(length(d))
    if (nd==1) {
        nn <- cbind(c(0L, 1L:(n-1)),
                    c(2L:n, 0L))
    } else {
        if (queen) {
            exd <- expand.grid(
                lapply(1L:nd, function(j)
                    c(1,d)[j]^(j-1)*(-1L:1L)))
            ll <- as.integer(rowSums(exd))
            gc(reset=TRUE)
            nn <- matrix(1L:n, n, length(ll))
            for (j in 1:ncol(nn))
                nn[, j] <- nn[, j] + ll[j]
            for (j in 1:nd) {
                ii <- rep(rep(1:d[j],
                              each=prod(c(1,d)[1:j])
                              ), prod(d[-(1:j)]))
                id <- which(ii==1)
                for (jj in which(exd[,j]<0))
                    nn[id, jj] <- 0L
                id <- which(ii==d[j])
                for (jj in which(exd[,j]>0))
                    nn[id, jj] <- 0L
            }
            if (!self) 
                nn[, which(ll==0)] <- 0L
        } else {
            nn <- matrix(0L, n, self + 2*nd)
            id <- 1L:n
            if (self) nn[, 1] <- id
            for (j in 1:nd) {
                ii <- rep(rep(1L:d[j],
                              each=prod(c(1,d)[1:j])
                              ), prod(d[-(1:j)]))
                nn[, self + 2*(j-1)+1] <- as.integer((
                    id-c(1,d)[j]^(j-1))*(ii!=1)) 
                nn[, self + 2*(j-1)+2] <- as.integer((
                    id+c(1,d)[j]^(j-1))*(ii!=d[j])) 
            }
        }
    }
    if (nb) {
        nn <- lapply(1L:n, function(i)
            nn[i, nn[i,]>0L])
        class(nn) <- 'nb'
        return(nn)
    } else {
        if (self) return(nn)
        else return(nn[, -which(ll==0)])
    }
}
