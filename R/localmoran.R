# Copyright 2001-18 by Roger Bivand, 2021 Jeff Sauer and Levi Wolf (conditional code)
#

localmoran <- function(x, listw, zero.policy=NULL, na.action=na.fail,
        conditional=FALSE, alternative = "greater", p.adjust.method="none",
        mlvar=TRUE, spChk=NULL, adjust.x=FALSE) {
        stopifnot(is.vector(x))
	if (!inherits(listw, "listw"))
		stop(paste(deparse(substitute(listw)), "is not a listw object"))
        if (is.null(zero.policy))
            zero.policy <- get("zeroPolicy", envir = .spdepOptions)
        stopifnot(is.logical(zero.policy))
	if (!is.null(attr(listw$neighbours, "self.included")) &&
		attr(listw$neighbours, "self.included"))
		stop("Self included among neighbours")
	if (is.null(spChk)) spChk <- get.spChkOption()
	if (spChk && !chkIDs(x, listw))
		stop("Check of data and weights ID integrity failed")
	if (!is.numeric(x))
		stop(paste(deparse(substitute(x)), "is not a numeric vector"))
	NAOK <- deparse(substitute(na.action)) == "na.pass"
	x <- na.action(x)
	na.act <- attr(x, "na.action")
        rn <- attr(listw, "region.id")
	if (!is.null(na.act)) {
	    subset <- !(1:length(listw$neighbours) %in% na.act)
	    listw <- subset(listw, subset, zero.policy=zero.policy)
	    excl <- class(na.act) == "exclude"
	}
	n <- length(listw$neighbours)
	if (n != length(x))stop("Different numbers of observations")
	res <- matrix(nrow=n, ncol=5)
        if (alternative == "two.sided") Prname <- "Pr(z != 0)"
        else if (alternative == "greater") Prname <- "Pr(z > 0)"
        else Prname <- "Pr(z < 0)"
	colnames(res) <- c("Ii", "E.Ii", "Var.Ii", "Z.Ii", Prname)
	if (adjust.x) {
          nc <- card(listw$neighbours) > 0L
	  xx <- mean(x[nc], na.rm=NAOK)
	} else {
	  xx <- mean(x, na.rm=NAOK)
	}
	z <- x - xx 
	lz <- lag.listw(listw, z, zero.policy=zero.policy, NAOK=NAOK)

	if (mlvar) {
          if (adjust.x) {
            s2 <- sum(z[nc]^2, na.rm=NAOK)/sum(nc)
          } else {
            s2 <- sum(z^2, na.rm=NAOK)/n
          }
	} else {
          if (adjust.x) {
            s2 <- sum(z[nc]^2, na.rm=NAOK)/(sum(nc)-1) 
          } else {
            s2 <- sum(z^2, na.rm=NAOK)/(n-1) 
          }
        }
	res[,1] <- (z/s2) * lz
	Wi <- sapply(listw$weights, sum) 
	if (conditional){	
	  m2 <- sum(z * z) / n
	  res[, 2] <- -(z ** 2 * Wi) / ((n - 1) * m2)
	} else {
	  res[, 2] <- -Wi / (n-1) 
	}  
	if (mlvar)  {
          if (adjust.x) {
            b2 <- (sum(z[nc]^4, na.rm=NAOK)/sum(nc))/(s2^2)
          } else {
            b2 <- (sum(z^4, na.rm=NAOK)/n)/(s2^2)
          }
	} else {
          if (adjust.x) {
            b2 <- (sum(z[nc]^4, na.rm=NAOK)/(sum(nc)-1))/(s2^2) 
          } else {
            b2 <- (sum(z^4, na.rm=NAOK)/(n-1))/(s2^2)
          }
        }
        
	Wi2 <- sapply(listw$weights, function(x) sum(x^2)) 
	A <- (n-b2) / (n-1)
	B <- (2*b2 - n) / ((n-1)*(n-2))
	if (conditional){
	  res[, 3] <- ((z / m2) ** 2 *
	                 (n / (n - 2)) *
	                 (Wi2 - (Wi ** 2 / (n - 1))) *
	                 (m2 - (z ** 2 / (n - 1))))
	} else {
	  res[,3] <- A*Wi2 + B*(Wi^2 - Wi2) - res[,2]^2
	}
# Changed to Sokal (1998) VIi
#	 C <- Wi^2 / ((n-1)^2) # == res[,2]^2
#	 Wikh2 <- sapply(listw$weights, function(x) {
#	   if(is.null(x)) 0 else 1 - (crossprod(x,x))
#	 })
#        res[,3] <- A*Wi2 + B*Wikh2 - C
	res[,4] <- (res[,1] - res[,2]) / sqrt(res[,3])
	if (alternative == "two.sided") {
	  pv <- 2 * pnorm(abs(res[,4]), lower.tail=FALSE)
	} else if (alternative == "greater") {
	  pv <- pnorm(res[,4], lower.tail=FALSE)
	} else {
	  pv <- pnorm(res[,4])
	}
	res[,5] <- p.adjustSP(pv, listw$neighbours, method=p.adjust.method)
	if (!is.null(na.act) && excl) {
		res <- naresid(na.act, res)
	}
        if (!is.null(rn)) rownames(res) <- rn
	attr(res, "call") <- match.call()
	if (!is.null(na.act)) attr(res, "na.action") <- na.act
	class(res) <- c("localmoran", class(res))
	res
}



