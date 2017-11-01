# Copyright 2002-3 by Hisaji ONO and Roger Bivand 
#
# General G Statistics
#
#
globalG.test <- function(x, listw, zero.policy=NULL,
	alternative="greater", spChk=NULL, adjust.n=TRUE, B1correct=TRUE) {
        if (is.null(zero.policy))
            zero.policy <- get("zeroPolicy", envir = .spdepOptions)
        stopifnot(is.logical(zero.policy))
        stopifnot(is.vector(x))
	alternative <- match.arg(alternative, c("greater", "less", "two.sided"))
	if (!inherits(listw, "listw"))
	stop(paste(deparse(substitute(listw)), "is not a listw object"))
	if (is.na(match(listw$style, c("B", "C", "U")))) 
          warning("Binary weights recommended (sepecially for distance bands)")
	if (!is.numeric(x))
	stop(paste(deparse(substitute(x)), "is not a numeric vector"))
	if (any(is.na(x))) stop(paste("NA in ", deparse(substitute(x))))
	if (any(x < 0.0)) 
		stop(paste("Negative value in ", deparse(substitute(x))))
	n <- length(listw$neighbours)
	if (n != length(x))stop("Different numbers of observations")
	if (is.null(spChk)) spChk <- get.spChkOption()
	if (spChk && !chkIDs(x, listw))
		stop("Check of data and weights ID integrity failed")

	wc <- spweights.constants(listw, zero.policy=zero.policy, 
		adjust.n=adjust.n)
	n1 <- n - 1
	n2 <- n - 2
	n3 <- n - 3
	nn <- n * n
	S0 <- wc$S0
	S1 <- wc$S1
	S2 <- wc$S2
	S02 <- S0*S0
	G <- (t(x) %*% lag.listw(listw, x, zero.policy=zero.policy)) /
		(sum(x %x% x) - (t(x) %*% x))

	E.G <- S0 / (n * n1)

	B0 <- ((nn - 3*n + 3)*S1) - (n*S2) + (3*S02)
# added 141222 to permit comparison with CrimeStat IV
	B1 <- -(((nn - n)*S1) - (2*n*S2) + (ifelse(B1correct, 6, 3)*S02))
	B2 <- -((2*n*S1) - ((n+3)*S2) + (6*S02))
	B3 <- (4*n1*S1) - (2*(n+1)*S2) + (8*S02)
	B4 <- S1 - S2 + S02
	sx <- sum(x)
	sx2 <- sum(x^2)
	sx3 <- sum(x^3)
	sx4 <- sum(x^4)

	var.G <- ((B0*(sx2^2) + B1*sx4 + B2*(sx^2)*sx2 + B3*sx*sx3 +
		 B4*(sx^4)) / ((((sx^2) - sx2)^2)*n*n1*n2*n3)) - (E.G^2)

	statistic <- (G - E.G) / sqrt(var.G)
	names(statistic) <- "standard deviate"
	method <- "Getis-Ord global G statistic"
	if (alternative == "two.sided") PrG <- 2 * pnorm(abs(statistic), 
# swirched -abs() to abs() 141121 RSB comment Tomasz Kossowski
		lower.tail=FALSE)
        else if (alternative == "greater")
            PrG <- pnorm(statistic, lower.tail=FALSE)
        else PrG <- pnorm(statistic)
	if (!is.finite(PrG) || PrG < 0 || PrG > 1) 
		warning("Out-of-range p-value: reconsider test arguments")
	vec <- c(G, E.G, var.G)
	names(vec) <- c("Global G statistic", "Expectation", "Variance")
	data.name <- paste(deparse(substitute(x)), "\nweights:",
	    deparse(substitute(listw)), "\n")
	res <- list(statistic=statistic, p.value=PrG, estimate=vec, 
	    alternative=alternative, data.name=data.name, method=method)
	class(res) <- "htest"
	res
}

