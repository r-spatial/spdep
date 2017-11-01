preAple <- function(x, listw, override_similarity_check=FALSE, useTrace=TRUE) {
    stopifnot(isTRUE(all.equal(mean(x), 0.0)))
    stopifnot(is.vector(x))
    W <- as(listw, "CsparseMatrix")
    n <- dim(W)[1]
    if (useTrace) {
        trWW <- sum(diag(W %*% W))
    } else {
        if (listw$style %in% c("W", "S") && !override_similarity_check) {
            can.sim <- can.be.simmed(listw)
            eig <- eigenw(similar.listw(listw))
        } else {
            can.sim <- FALSE
            eig <- eigenw(listw)
        }
        if (is.complex(eig)) trWW <- Re(crossprod(eig))
        else trWW <- crossprod(eig)
# modified 110414 RSB
#            eig <- Re(eig)
#        trWW <- crossprod(eig)
    }
    corterm <- (trWW/n) * Diagonal(n)
    corterm <- as(corterm, "CsparseMatrix")
    WU <- ((W + t(W))/2)
    W2 <- crossprod(W) + corterm
    res <- list(W=W, corterm=corterm, W2=W2, WU=WU, n=n)
    res
}

inAple <- function(x, pre) {
    xwx <- crossprod(x, (pre$WU %*% x))
    xwwx <- crossprod(x, (pre$W2 %*% x))
    res <- c(as.matrix(xwx/xwwx))
    res
}

aple <- function(x, listw, override_similarity_check=FALSE, useTrace=TRUE) {
    pre <- preAple(x=x, listw=listw,
        override_similarity_check=override_similarity_check, useTrace=useTrace)
    res <- inAple(x=x, pre=pre)
    res
}


