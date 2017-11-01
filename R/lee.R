#Lee (2001)'s bivariate association statistic
#Based on code by Roger Bivand for moran's I

lee <- function(x, y, listw, n, S2=NULL, zero.policy=NULL, NAOK=FALSE) {
        if (is.null(zero.policy))
            zero.policy <- get("zeroPolicy", envir = .spdepOptions)
        stopifnot(is.logical(zero.policy))
        n1 <- length(listw$neighbours)
        x <- c(x)
	y <- c(y)
        if (n1 != length(x) | n1 != length(y)) stop("objects of different length")
        xx <- mean(x, na.rm=NAOK)
        yy <- mean(y, na.rm=NAOK)

        z  <- x - xx
        zz <- sum(z^2, na.rm=NAOK)
	zy <- y - yy
	zzy<- sum(zy^2, na.rm=NAOK)

#        K <- NA#(length(x)*sum(z^4, na.rm=NAOK))/(zz^2)
        lz <- lag.listw(listw, z, zero.policy=zero.policy, NAOK=NAOK)
        lzy <- lag.listw(listw, zy, zero.policy=zero.policy, NAOK=NAOK)
##       I <- (n / S0) * ((t(z) %*% lz) / zz)
#        I <- (n / S0) * ((sum(z*lz, na.rm=NAOK)) / zz)
#        res <- list(I=I, K=K)
#        res

	if(is.null(S2))
		S2<-sum ( (unlist(lapply(listw$weights, sum)))^2 )

	L<- (n/S2)* (sum(lz*lzy))/(sqrt(zz)*sqrt(zzy))

	localL<-n*lz*lzy/(sqrt(zz)*sqrt(zzy))

#	res<-list(L=L, K=K, localL=localL)
	res<-list(L=L, localL=localL)
	return(res)
}
