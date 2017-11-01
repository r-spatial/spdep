# Copyright 2001-6 by Roger Bivand, bugfix large n Ronnie Babigumira
#

joincount <- function(dums, listw) {
	nc <- which(colSums(dums) > 1)
#	n <- length(listw$neighbours)
	cardnb <- card(listw$neighbours)
	res <- as.numeric(rep(0, ncol(dums)))
	for (lev in nc) {
		res[lev] <- .Call("jcintern", listw$neighbours,
			listw$weights, as.integer(dums[,lev]),
			as.integer(cardnb), PACKAGE="spdep")
	}
	res
}

joincount.test <- function(fx, listw, zero.policy=NULL,
	alternative="greater", #adjust.n=TRUE, 
	spChk=NULL, adjust.n=TRUE) {
        if (is.null(zero.policy))
            zero.policy <- get("zeroPolicy", envir = .spdepOptions)
        stopifnot(is.logical(zero.policy))
	alternative <- match.arg(alternative, c("greater", "less", "two.sided"))
	if (!inherits(listw, "listw")) stop(paste(deparse(substitute(listw)),
		"is not a listw object"))
	if (!is.factor(fx)) stop(paste(deparse(substitute(x)),
		"is not a factor"))
	if (any(is.na(fx))) stop("NA in factor")
	n <- length(listw$neighbours)
	if (n != length(fx)) stop("objects of different length")
	cards <- card(listw$neighbours)
	if (!zero.policy && any(cards == 0))
		stop("regions with no neighbours found")
	if (is.null(spChk)) spChk <- get.spChkOption()
	if (spChk && !chkIDs(fx, listw))
		stop("Check of data and weights ID integrity failed")

	wc <- spweights.constants(listw, zero.policy=zero.policy, 
		adjust.n=adjust.n)
	S02 <- wc$S0*wc$S0

	ff <- ~ fx - 1
	dums <- model.matrix(ff, model.frame(ff))
	BB <- joincount(dums, listw)
	nBB <- length(BB)
	if (nBB < 1) stop("non-positive BB length")
	res <- vector(mode="list", length=nBB)
	tab <- table(fx)
	BB5 <- 0.5 * BB
	ntab <- as.numeric(as.vector(tab))
# comment and bug report by Tomoki NAKAYA about no-neighbour observations
#	if (adjust.n) {
		N <- wc$n
#	} else {
#		N <- n
#		wc$n1 <- N-1
#		wc$n2 <- N-2
#		wc$n3 <- N-3
#	}
	Ejc <- (wc$S0*(ntab*(ntab-1))) / (2*N*wc$n1)
	Vjc <- (wc$S1*(ntab*(ntab-1))) / (N*wc$n1)
	Vjc <- Vjc + (((wc$S2 - 2*wc$S1)*ntab*(ntab-1)*(ntab-2)) /
		(N*wc$n1*wc$n2))
	Vjc <- Vjc + (((S02 + wc$S1 - wc$S2)*ntab*(ntab-1)*(ntab-2)*
		(ntab-3)) / (N*wc$n1*wc$n2*wc$n3))
	Vjc <- (0.25 * Vjc) - Ejc^2
	for (i in 1:nBB) {
		estimate <- c(BB5[i], Ejc[i], Vjc[i])
		names(estimate) <- c("Same colour statistic",
			"Expectation", "Variance")
		statistic <- (BB5[i] - Ejc[i]) / sqrt(Vjc[i])
		names(statistic) <- paste("Std. deviate for", names(tab)[i])
		p.value <- NA
		if (is.finite(statistic)) {
		    if (alternative == "two.sided") 
			p.value <- 2 * pnorm(abs(statistic), lower.tail=FALSE)
		    else if (alternative == "greater")
			p.value <- pnorm(statistic, lower.tail=FALSE)
		    else p.value <- pnorm(statistic)
		    if (!is.finite(p.value) || p.value < 0 || p.value > 1) 
		      warning("Out-of-range p-value: reconsider test arguments")
		}
		method <- "Join count test under nonfree sampling"
		data.name <- paste(deparse(substitute(fx)), "\nweights:",
			deparse(substitute(listw)), "\n")
		res[[i]] <- list(statistic=statistic, p.value=p.value,
			estimate=estimate, method=method,
			alternative=alternative, data.name=data.name)
		class(res[[i]]) <- "htest"
	}
	class(res) <- "jclist"
	res
}

print.jclist <- function(x, ...) {
	for (i in seq(along=x)) print(x[[i]], ...)
	invisible(x)
}

joincount.mc <- function(fx, listw, nsim, zero.policy=FALSE,
	alternative="greater", spChk=NULL) {
	alternative <- match.arg(alternative, c("greater", "less"))
	if(!inherits(listw, "listw")) stop(paste(deparse(substitute(listw)),
		"is not a listw object"))
	if(!is.factor(fx)) stop(paste(deparse(substitute(fx)),
		"is not a factor"))
	if(missing(nsim)) stop("nsim must be given")
	if (any(is.na(fx))) stop("NA in factor")
	n <- length(listw$neighbours)
	if (n != length(fx)) stop("objects of different length")
	cards <- card(listw$neighbours)
	if (!zero.policy && any(cards == 0))
		stop("regions with no neighbours found")
	if (is.null(spChk)) spChk <- get.spChkOption()
	if (spChk && !chkIDs(fx, listw))
		stop("Check of data and weights ID integrity failed")
        gamres <- suppressWarnings(nsim > gamma(n + 1))
        if (gamres) stop("nsim too large for this number of observations")

	ff <- ~ fx - 1
	dums <- model.matrix(ff, model.frame(ff))
	nc <- ncol(dums)
	if (nc < 1) stop("non-positive nc")
	if (nsim < 1) stop("non-positive nsim")
	res <- matrix(0, nrow=nsim+1, ncol=nc)
	res[nsim+1,] <- 0.5 * joincount(dums, listw)
	tab <- table(fx)
	for (i in 1:nsim) {
		fxi <- sample(fx)
		ff <- ~ fxi - 1
		dums <- model.matrix(ff, model.frame(ff))
		res[i,] <- 0.5 * joincount(dums, listw)
	}
	rankres <- apply(res, 2, rank)
	xrank <- rankres[nrow(rankres),]
	lres <- vector(mode="list", length=nc)
	for (i in 1:nc) {
		statistic <- res[nrow(res), i]
		names(statistic) <- paste("Join-count statistic for",
			names(tab)[i])
		parameter <- xrank[i]
		names(parameter) <- "rank of observed statistic"
		diff <- nsim - xrank[i]
		diff <- ifelse(diff > 0, diff, 0)
        	if (alternative == "less") 
        		pval <- punif((diff + 1)/(nsim + 1), lower.tail=FALSE)
    		else if (alternative == "greater") 
        		pval <- punif((diff + 1)/(nsim + 1))
		if (!is.finite(pval) || pval < 0 || pval > 1) 
		    warning("Out-of-range p-value: reconsider test arguments")

		method <- "Monte-Carlo simulation of join-count statistic"
		data.name <- paste(deparse(substitute(fx)), "\nweights:",
			deparse(substitute(listw)),
			"\nnumber of simulations + 1:", nsim+1, "\n")
		estimate <- c(mean(res[-(nrow(res)), i]),
			var(res[-(nrow(res)), i]))
		names(estimate) <- c("mean of simulation",
			"variance of simulation")
		lres[[i]] <- list(statistic=statistic, parameter=parameter,
			method=method, data.name=data.name, p.value=pval, 
			alternative=alternative, estimate=estimate, res=res[,i])
		class(lres[[i]]) <- c("htest", "mc.sim")
		
	}
	class(lres) <- "jclist"
	lres
}



joincount.multi <- function(fx, listw, zero.policy=FALSE, #adjust.n=TRUE,
	spChk=NULL, adjust.n=TRUE) {
	if(!inherits(listw, "listw")) stop(paste(deparse(substitute(listw)),
		"is not a listw object"))
	if(!is.factor(fx)) stop(paste(deparse(substitute(fx)),
		"is not a factor"))
	if (any(is.na(fx))) stop("NA in factor")
	n <- length(listw$neighbours)
	if (n != length(fx)) stop("objects of different length")
	cards <- card(listw$neighbours)
	if (!zero.policy && any(cards == 0))
		stop("regions with no neighbours found")
	if (is.null(spChk)) spChk <- get.spChkOption()
	if (spChk && !chkIDs(fx, listw))
		stop("Check of data and weights ID integrity failed")
	ifx <- as.integer(fx)
	k <- length(levels(fx))
	if (k < 2) stop("must be at least two levels in factor")

	sn <- listw2sn(listw)
	y <- factor(paste(ifx[sn[,1]], ":", ifx[sn[,2]], sep=""), 
		levels=as.vector(outer(1:k, 1:k, 
			FUN=function(X,Y) paste(X,Y,sep=":"))))
	res <- matrix(tapply(sn[,3], y, sum), ncol=k)/2
		
	res[is.na(res)] <- 0
	rownames(res) <- colnames(res) <- levels(fx)

	tab <- table(fx)
	ntab <- as.numeric(as.vector(tab))
	wc <- spweights.constants(listw, zero.policy=zero.policy, 
		adjust.n=adjust.n)
# comment and bug report by Tomoki NAKAYA about no-neighbour observations
#	if (adjust.n) {
		N <- wc$n
#	} else {
#		N <- n
#		wc$n1 <- N-1
#		wc$n2 <- N-2
#		wc$n3 <- N-3
#	}
	S02 <- wc$S0*wc$S0

	Ejc <- (wc$S0*(ntab*(ntab-1))) / (2*N*wc$n1)

	Vjc <- (wc$S1*(ntab*(ntab-1))) / (N*wc$n1)
	Vjc <- Vjc + (((wc$S2 - 2*wc$S1)*ntab*(ntab-1)*(ntab-2)) /
		(N*wc$n1*wc$n2))
	Vjc <- Vjc + (((S02 + wc$S1 - wc$S2)*ntab*(ntab-1)*(ntab-2)*
		(ntab-3)) / (N*wc$n1*wc$n2*wc$n3))
	Vjc <- (0.25 * Vjc) - Ejc^2

	nrns <- function(x, op="*") {
		k <- length(x)
		res <- numeric(((k^2) - k)/2)
		ii <- 1
		for (i in 2:k) {
			for (j in 1:(i-1)) {
				if (is.character(op) && op == "*") {
					res[ii] <- x[i]*x[j]
				} else if (is.character(op) && op == "+") {
					res[ii] <- x[i]+x[j]
				}
				ii <- ii+1
			}
		}
		res
	}

	ldiag <- numeric(((k^2) - k)/2)
	diffcolnames <- character(((k^2) - k)/2)
	ii <- 1
	for (i in 2:k) {
		for (j in 1:(i-1)) {
			ldiag[ii] <- res[i,j] + res[j,i]
			diffcolnames[ii] <- paste(levels(fx)[i], 
				levels(fx)[j], sep=":")
			ii <- ii+1
		}
	}

	Exp <- (wc$S0*(nrns(ntab, op="*"))) / (N*wc$n1)
	Var <- (2*wc$S1*nrns(ntab, op="*"))/(N*wc$n1)
	Var <- Var + (((wc$S2 - 2*wc$S1)*nrns(ntab, op="*")*
		(nrns(ntab, op="+")-2))/(N*wc$n1*wc$n2))
	Var <- Var + ((4*(S02 + wc$S1 - wc$S2)*nrns((ntab*(ntab-1)), op="*")) /
		(N*wc$n1*wc$n2*wc$n3))
	Var <- (0.25 * Var) - Exp^2
	Jtot <- sum(ldiag)
	JtotExp <- sum(Exp)
	Jvar <- ((wc$S2/(N*wc$n1))-((4*(S02 + wc$S1 - wc$S2)*wc$n1) /
		(N*wc$n1*wc$n2*wc$n3)))*sum(nrns(ntab, op="*"))
	Jvar <- Jvar + 4*(((wc$S1 - wc$S2)/(N*wc$n1*wc$n2*wc$n3)) + 
		((2*S02*(2*n-3))/((N*wc$n1)*(N*wc$n1*wc$n2*wc$n3))))*
		sum(nrns(ntab^2, op="*"))
	if(k>2) {
		ntnsnr <- as.numeric(0)
		for (r in 1:(k-2)) {
			for (s in (r+1):(k-1)) {
				for (t in (s+1):(k)) {
					ntnsnr <- ntnsnr +
						ntab[r]*ntab[s]*ntab[t]
				}
			}
		}
		Jvar <- Jvar + (((2*wc$S1 - 5*wc$S2)/(N*wc$n1*wc$n2))+
		((12*(S02 + wc$S1 - wc$S2))/(N*wc$n1*wc$n2*wc$n3))+
		((8*S02)/((N*wc$n1*wc$n2)*wc$n1)))*ntnsnr
	}
	if(k>3) {
		nuntnsnr <- as.numeric(0)
		for (r in 1:(k-3)) {
			for (s in (r+1):(k-2)) {
				for (t in (s+1):(k-1)) {
					for (u in (t+1):(k)) {
						nuntnsnr <- nuntnsnr +
						ntab[r]*ntab[s]*ntab[t]*ntab[u]
					}
				}
			}
		}
		Jvar <- Jvar - 8*(((wc$S1 - wc$S2)/(N*wc$n1*wc$n2*wc$n3))+
		((2*S02*(2*N-3))/((N*wc$n1)*(N*wc$n1*wc$n2*wc$n3))))*nuntnsnr
	}
	Jvar <- (0.25 * Jvar)
	statistic <- (c(diag(res), ldiag, Jtot) - c(Ejc, Exp, JtotExp)) / 
		sqrt(c(Vjc, Var, Jvar))
	lres <- cbind(c(diag(res), ldiag, Jtot), c(Ejc, Exp, JtotExp), 
		c(Vjc, Var, Jvar), statistic)
	colnames(lres) <- c("Joincount", "Expected", "Variance", 
		"z-value")
	rownames(lres) <- c(paste(levels(fx), ":", levels(fx), sep=""), 
		diffcolnames, "Jtot")
	class(lres) <- c("jcmulti", "matrix")
	lres
}

print.jcmulti <- function(x, ...) {
	printCoefmat(x, ...)
}




