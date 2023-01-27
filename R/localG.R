# Copyright 2001-18 by Roger Bivand 
#

localG <- function(x, listw, zero.policy=NULL, spChk=NULL, return_internals=TRUE, GeoDa=FALSE, alternative = "two.sided") {
	if (!inherits(listw, "listw"))
		stop(paste(deparse(substitute(listw)), "is not a listw object"))
	if (!is.numeric(x))
		stop(paste(deparse(substitute(x)), "is not a numeric vector"))
        if (is.null(zero.policy))
            zero.policy <- get("zeroPolicy", envir = .spdepOptions)
        stopifnot(is.logical(zero.policy))
        stopifnot(is.vector(x))
	if (any(is.na(x))) stop(paste("NA in ", deparse(substitute(x))))
	n <- length(listw$neighbours)
	if (n != length(x))stop("Different numbers of observations")
	if (is.null(spChk)) spChk <- get.spChkOption()
	if (spChk && !chkIDs(x, listw))
		stop("Check of data and weights ID integrity failed")
	gstari <- FALSE
	if (!is.null(attr(listw$neighbours, "self.included")) &&
	    attr(listw$neighbours, "self.included")) gstari <- TRUE
        Gi_str <- ifelse(gstari, "G*i", "Gi")
	lx <- lag.listw(listw, x, zero.policy=zero.policy)
        if (GeoDa) {
            if (gstari) {
                x_star <- sum(x[card(listw$neighbours) > 1L])
            } else {
                x_star <- sum(x[card(listw$neighbours) > 0L])
            }
        } else {
            x_star <- sum(x)
        }
	if (gstari) {
		xibar <- rep(mean(x), n)
		si2 <- rep(sum(scale(x, scale=FALSE)^2)/n, n)
	} else {
		xibar <- (rep(sum(x),n) - x) / (n - 1)
		si2 <- ((rep(sum(x^2),n) - x^2) / (n - 1)) - xibar^2
	}
	Wi <- sapply(listw$weights, sum)
	S1i <- sapply(listw$weights, function(x) sum(x^2))
        EG <- Wi*xibar
	res <- (lx - EG)
	if (gstari) {
                VG <- si2*((n*S1i - Wi^2)/(n-1))
	} else {
                VG <- si2*(((n-1)*S1i - Wi^2)/(n-2))
	}
        res <- res / sqrt(VG)
        if (return_internals) {
            alternative <- match.arg(alternative, c("two.sided", "greater",
                "less"))
            if (alternative == "two.sided") Prname <- paste0("Pr(z != E(",
                Gi_str, "))")
            else if (alternative == "greater") Prname <- paste0("Pr(z > E(",
                Gi_str, "))")
            else Prname <- paste0("Pr(z < E(", Gi_str, "))")
            if (alternative == "two.sided") {
	        pv <- 2 * pnorm(abs(res), lower.tail=FALSE)
	    } else if (alternative == "greater") {
	        pv <- pnorm(res, lower.tail=FALSE)
	    } else {
	        pv <- pnorm(res)
	    }
            if (gstari) {
                ints <- cbind(G=lx/x_star, EG=EG/x_star, VG=VG/x_star^2,
                    ZG=res, pv=pv)
            } else {
                ints <- cbind(G=lx/(x_star-c(x)), EG=EG/(x_star-c(x)),
                    VG=VG/(x_star-c(x))^2, ZG=res, pv=pv)
            }
            colnames(ints) <- c(paste(c("", "E(", "V(", "Z("), Gi_str,
                c("", ")", ")", ")"), sep=""), Prname)
            attr(res, "internals") <- ints
	}
        attr(res, "cluster") <- cut(x, c(-Inf, mean(x), Inf),
            labels = c("Low", "High"))
        attr(res, "gstari") <- gstari
	attr(res, "call") <- match.call()
	class(res) <- "localG"
	res
}

# "localG" cluster c("Low", "High")
