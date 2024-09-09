# Copyright 2002-2024 by Hisaji ONO and Roger Bivand 
#
# General G Statistics
#
#
globalG.test <- function(x, listw, zero.policy=attr(listw, "zero.policy"),
	alternative="greater", spChk=NULL, adjust.n=TRUE, B1correct=TRUE,
        adjust.x=TRUE, Arc_all_x=FALSE, na.action=na.fail) {
        if (is.null(zero.policy))
            zero.policy <- get.ZeroPolicyOption()
        stopifnot(is.logical(zero.policy))
        stopifnot(is.vector(x))
	alternative <- match.arg(alternative, c("greater", "less", "two.sided"))
	wname <- deparse(substitute(listw))
	if (!inherits(listw, "listw")) stop(wname, "is not a listw object")
	if (is.na(match(listw$style, c("B", "C", "U")))) 
          warning("Binary weights recommended (especially for distance bands)")
        xname <- deparse(substitute(x))
	if(!is.numeric(x)) stop(xname, " is not a numeric vector")
	if (deparse(substitute(na.action)) == "na.pass")
	    stop("na.pass not permitted")
	x <- na.action(x)
	na.act <- attr(x, "na.action")
	if (!is.null(na.act)) {
	    subset <- !(1:length(listw$neighbours) %in% na.act)
	    listw <- subset(listw, subset, zero.policy=zero.policy)
	}
	if (any(x < 0.0)) stop("Negative value in ", xname)
        
	n <- length(listw$neighbours)
	if (n != length(x)) stop("Different numbers of observations")
	if (is.null(spChk)) spChk <- get.spChkOption()
	if (spChk && !chkIDs(x, listw))
		stop("Check of data and weights ID integrity failed")

	wc <- spweights.constants(listw, zero.policy=zero.policy, 
		adjust.n=adjust.n)
        n <- wc$n
	n1 <- wc$n1
	n2 <- wc$n2
	n3 <- wc$n3
	nn <- wc$nn
	S0 <- wc$S0
	S1 <- wc$S1
	S2 <- wc$S2
	S02 <- S0*S0

        if (adjust.x) ax <- x[card(listw$neighbours) > 0L]
        else ax <- x

        numer <- (t(x) %*% lag.listw(listw, x, zero.policy=zero.policy))
        if (Arc_all_x) {
          outer_prod <- sapply(ax, function(axi) sum(x*axi))
#(sum(x %x% ax)
        } else {
          outer_prod <- sapply(ax, function(axi) sum(ax*axi))
#(sum(ax %x% ax)
        }
        denom <- sum(outer_prod) - crossprod(ax)
	G <- numer / denom
		
	E.G <- S0 / (n * n1)

	B0 <- ((nn - 3*n + 3)*S1) - (n*S2) + (3*S02)
# added 141222 to permit comparison with CrimeStat IV
	B1 <- -(((nn - n)*S1) - (2*n*S2) + (ifelse(B1correct, 6, 3)*S02))
	B2 <- -((2*n*S1) - ((n+3)*S2) + (6*S02))
	B3 <- (4*n1*S1) - (2*(n+1)*S2) + (8*S02)
	B4 <- S1 - S2 + S02
	sx <- sum(ax)
	sx2 <- sum(ax^2)
	sx3 <- sum(ax^3)
	sx4 <- sum(ax^4)

	var.G <- ((B0*(sx2^2) + B1*sx4 + B2*(sx^2)*sx2 + B3*sx*sx3 +
		 B4*(sx^4)) / ((((sx^2) - sx2)^2)*n*n1*n2*n3)) - (E.G^2)

	statistic <- (G - E.G) / sqrt(var.G)
	names(statistic) <- "standard deviate"
	method <- "Getis-Ord global G statistic"
	if (alternative == "two.sided") PrG <- 2 * pnorm(abs(statistic), 
# switched -abs() to abs() 141121 RSB comment Tomasz Kossowski
		lower.tail=FALSE)
        else if (alternative == "greater")
            PrG <- pnorm(statistic, lower.tail=FALSE)
        else PrG <- pnorm(statistic)
	if (!is.finite(PrG) || PrG < 0 || PrG > 1) 
		warning("Out-of-range p-value: reconsider test arguments")
	vec <- c(G, E.G, var.G)
	names(vec) <- c("Global G statistic", "Expectation", "Variance")
	data.name <- paste(xname, "\nweights:", wname, ifelse(is.null(na.act),
            "", paste("\nomitted:", paste(na.act, collapse=", "))),
            ifelse(adjust.n && isTRUE(any(sum(card(listw$neighbours) == 0L))),
            "\nn reduced by no-neighbour observations", ""), "\n")
	res <- list(statistic=statistic, p.value=PrG, estimate=vec, 
	    alternative=alternative, data.name=data.name, method=method)
	class(res) <- "htest"
	res
}

