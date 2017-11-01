# Copyright 2001-2010 by Roger Bivand 
#

lm.morantest <- function(model, listw, zero.policy=NULL, 
	    alternative = "greater", spChk=NULL, resfun=weighted.residuals,		    naSubset=TRUE) {
        if (is.null(zero.policy))
            zero.policy <- get("zeroPolicy", envir = .spdepOptions)
        stopifnot(is.logical(zero.policy))
	alternative <- match.arg(alternative, c("greater", "less", "two.sided"))
	listw_name <- deparse(substitute(listw))
	if (!inherits(listw, "listw")) stop(paste(listw_name,
		"is not a listw object"))
	if(!inherits(model, "lm")) stop(paste(deparse(substitute(model)),
		"not an lm object"))
	na.act <- model$na.action
	if (!is.null(na.act) && naSubset) {
	    subset <- !(1:length(listw$neighbours) %in% na.act)
	    listw <- subset(listw, subset, zero.policy=zero.policy)
	}
# 101124 Aleksandr Andreev
 	N <- as.double(length(listw$neighbours))
	u <- resfun(model)
	if (N != length(u)) 
            stop("objects of different length")
	if (is.null(spChk)) spChk <- get.spChkOption()
	if (spChk && !chkIDs(u, listw))
		stop("Check of data and weights ID integrity failed")
	
	u <- as.vector(u)
	listw.U <- listw2U(listw)

	S0 <- sum(unlist(listw.U$weights))
	S1 <- 0.5 * sum((2*unlist(listw.U$weights))^2)
	lu <- lag.listw(listw.U, u, zero.policy=zero.policy)
# 101125 Aleksandr Andreev
	if (zero.policy)
            N <- as.double(length(which(card(listw$neighbours) > 0L)))
	I <- (N/S0) * ((t(u) %*% lu) / (t(u) %*% u))
	p <- model$rank
	p1 <- 1:p
	nacoefs <- which(is.na(coefficients(model)))
	XtXinv <- chol2inv(model$qr$qr[p1, p1, drop = FALSE])
	X <- model.matrix(terms(model), model.frame(model))
# fixed after looking at TOWN dummy in Boston data
	if (length(nacoefs) > 0L) X <- X[,-nacoefs]
	if (!is.null(wts <- weights(model))) {
#		X <- sqrt(diag(wts)) %*% X
		X <- drop(t(sapply(1:length(wts), 
			function(i) sqrt(wts[i])*X[i,])))
	}
# Cliff/Ord 1981, p. 203
	Z <- lag.listw(listw.U, X, zero.policy=zero.policy)
	C1 <- t(X) %*% Z
	trA <- (sum(diag(XtXinv %*% C1)))
	EI <- -((N * trA) / ((N-p) * S0))
# minus changed from trA to EI (Luis Galvis, Dec 2, 2003)
	C2 <- t(Z) %*% Z
	C3 <- XtXinv %*% C1
	trA2 <- sum(diag(C3 %*% C3))
	trB <- sum(diag(4*(XtXinv %*% C2)))
	VI <- (((N*N)/((S0*S0)*(N-p)*(N-p+2))) *
		(S1 + 2*trA2 - trB - ((2*(trA^2))/(N-p))))
	ZI <- (I - EI) / sqrt(VI)
    	if (alternative == "two.sided") pv <- 2 * pnorm(abs(ZI), 
		lower.tail=FALSE)
    	else if (alternative == "greater")
	        pv <- pnorm(ZI, lower.tail=FALSE)
    	else pv <- pnorm(ZI)
	if (!is.finite(pv) || pv < 0 || pv > 1) 
		warning("Out-of-range p-value: reconsider test arguments")
    	statistic <- ZI
    	attr(statistic, "names") <- "Moran I statistic standard deviate"
    	p.value <- pv
    	estimate <- c(I, EI, VI)
    	attr(estimate, "names") <- c("Observed Moran I", "Expectation",
	    "Variance")
    	method <- "Global Moran I for regression residuals"
    	data.name <- paste("\n", paste(strwrap(paste("model: ",
	    gsub("[[:space:]]+", " ", 
	    paste(deparse(model$call), sep="", collapse="")))), collapse="\n"),
    	    "\nweights: ", listw_name, "\n", sep="")
    	res <- list(statistic = statistic, p.value = p.value,
	       estimate = estimate, method = method,
		alternative = alternative, data.name = data.name)
	class(res) <- "htest"
    	res
}
