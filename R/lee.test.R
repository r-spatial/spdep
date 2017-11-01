# Copyright 2014 by Roger Bivand , Virgilio Gómez-Rubio
#

lee.test <- function(x, y, listw, #randomisation=TRUE, 
	zero.policy=NULL,
	alternative="greater", 
	#rank = FALSE, 
	na.action=na.fail, spChk=NULL#, 
	#adjust.n=TRUE
	) {
	alternative <- match.arg(alternative, c("greater", "less", "two.sided"))
	if (!inherits(listw, "listw")) stop(paste(deparse(substitute(listw)),
		"is not a listw object"))
	if (!is.numeric(x)) stop(paste(deparse(substitute(x)),
		"is not a numeric vector"))
	if (!is.numeric(y)) stop(paste(deparse(substitute(y)),
		"is not a numeric vector"))
        if (is.null(zero.policy))
            zero.policy <- get("zeroPolicy", envir = .spdepOptions)
        stopifnot(is.logical(zero.policy))
	if (is.null(spChk)) spChk <- get.spChkOption()
	if (spChk && !chkIDs(x, listw) && !chkIDs(y, listw))
		stop("Check of data and weights ID integrity failed")
#	if (any(is.na(x))) stop("NA in X")
	xname <- deparse(substitute(x))
	yname <- deparse(substitute(y))
	wname <- deparse(substitute(listw))
	NAOK <- deparse(substitute(na.action)) == "na.pass"

	#Check NA's in both vectors
        na.act <- attr(na.action(cbind(x,y)), "na.action")

        x[na.act]<-NA
        y[na.act]<-NA

        x<-na.action(x)
        y<-na.action(y)

	if (!is.null(na.act)) {
	    subset <- !(1:length(listw$neighbours) %in% na.act)
	    listw <- subset(listw, subset, zero.policy=zero.policy)
	}
	n <- length(listw$neighbours)
	if (n != length(x)) stop("objects of different length")
	
#	wc <- spweights.constants(listw, zero.policy=zero.policy, 
#		adjust.n=adjust.n)
#	S02 <- wc$S0*wc$S0

        S2<-sum ( (unlist(lapply(listw$weights, sum)))^2 )

	res <- lee(x, y, listw, n, S2, zero.policy=zero.policy, 
		NAOK=NAOK)
#	I <- res$I
#	K <- res$K
	L<-res$L


#	if (rank) K <- (3*(3*wc$n^2 -7))/(5*(wc$n^2 - 1))
#	EI <- (-1) / wc$n1

	#Compute asymptotic mean EI and variance VI
	W <- as(listw, "CsparseMatrix")

	#See Lee (2004)
	dEG<-EGamma(x, y, W)
	dVarG<-VarGamma(x, y, W)

	EL<-dEG$EGammaon+dEG$EGammaoff
	VL<-dVarG$varGammaon+dVarG$varGammaoff+2*dVarG$varGammaonoff

#	if(randomisation) {
#		VI <- wc$n*(wc$S1*(wc$nn - 3*wc$n + 3) - wc$n*wc$S2 + 3*S02)
#		tmp <- K*(wc$S1*(wc$nn - wc$n) - 2*wc$n*wc$S2 + 6*S02)
#                if (tmp > VI) warning("Kurtosis overflow,\ndistribution of variable does not meet test assumptions")
#		VI <- (VI - tmp) / (wc$n1*wc$n2*wc$n3*S02)
#                tmp <- (VI - EI^2)
#                if (tmp < 0) warning("Negative variance,\ndistribution of variable does not meet test assumptions")
#		VI <- tmp
#	} else {
#		VI <- (wc$nn*wc$S1 - wc$n*wc$S2 + 3*S02) / (S02*(wc$nn - 1))
#                tmp <- (VI - EI^2)
#                if (tmp < 0) warning("Negative variance,\ndistribution of variable does not meet test assumptions")
#		VI <- tmp
#	}

	ZL <- (L - EL) / sqrt(VL)
	statistic <- ZL
	names(statistic) <- "Lee's L statistic standard deviate"
        if (alternative == "two.sided") 
		PrL <- 2 * pnorm(abs(ZL), lower.tail=FALSE)
        else if (alternative == "greater")
            PrL <- pnorm(ZL, lower.tail=FALSE)
        else PrL <- pnorm(ZL)
	if (!is.finite(PrL) || PrL < 0 || PrL > 1) 
		warning("Out-of-range p-value: reconsider test arguments")
	vec <- c(L, EL, VL)
	names(vec) <- c("Lee's L statistic", "Expectation", "Variance")
#	method <- paste("Lee's L test under", ifelse(randomisation,
#	    "randomisation", "normality"))
	method <- "Lee's L statistic randomisation"
#	data.name <- paste(xname, ifelse(rank,
#		"using rank correction",""), "\nweights:",
#		wname, ifelse(is.null(na.act), "", paste("\nomitted:", 
#	    paste(na.act, collapse=", "))), "\n")
	data.name <- paste(xname, ", ", yname,
		"\nweights:",
		wname, ifelse(is.null(na.act), "", paste("\nomitted:", 
	    paste(na.act, collapse=", "))), "\n")
	res <- list(statistic=statistic, p.value=PrL, estimate=vec, 
	    alternative=alternative, method=method, data.name=data.name)
	if (!is.null(na.act)) attr(res, "na.action") <- na.act
	class(res) <- "htest"
	res
}

