# Copyright 2001-7 by Roger Bivand 
#

lm.LMtests <- function(model, listw, zero.policy=NULL, test="LMerr",
	spChk=NULL, naSubset=TRUE) {

	if (inherits(model, "lm")) na.act <- model$na.action
	else na.act <- attr(model, "na.action")

	listw_name <- deparse(substitute(listw))

	if (!inherits(listw, "listw")) stop(paste(listw_name,
		"is not a listw object"))
        if (is.null(zero.policy))
            zero.policy <- get("zeroPolicy", envir = .spdepOptions)
        stopifnot(is.logical(zero.policy))
	if (!is.null(na.act) && naSubset) {
	    subset <- !(1:length(listw$neighbours) %in% na.act)
	    listw <- subset(listw, subset, zero.policy=zero.policy)
	}

	if (length(test) == 1L && test[1] == "LMerr") {
		res <- lm.LMErr(model=model, listw=listw, 
			zero.policy=zero.policy, spChk=spChk) 
		if (inherits(model, "lm")) res$data.name <- paste("\n", 
	    		paste(strwrap(paste("model: ",
	    		gsub("[ ]+", " ", paste(deparse(model$call), 
	    		sep="", collapse="")))), collapse="\n"),
    	     		"\nweights: ", listw_name, "\n", sep="")
		else res$data.name <- paste("\nresiduals: ", 
			deparse(substitute(model)), "\nweights: ", 
			listw_name, "\n", sep="")
		tres <- vector(mode="list", length=1)
		names(tres) <- test
		tres[[1]] <- res
		class(tres) <- "LMtestlist"
		return(tres)
	}
	if(!inherits(model, "lm")) stop(paste(deparse(substitute(model)),
		"not an lm object"))
	N <- length(listw$neighbours)
	u <- resid(model)
	if (N != length(u)) stop("objects of different length")
	if (is.null(spChk)) spChk <- get.spChkOption()
	if (spChk && !chkIDs(u, listw))
		stop("Check of data and weights ID integrity failed")
	u <- as.vector(u)

	if (is.null(attr(listw$weights, "W")) || !attr(listw$weights, "W"))
		warning("Spatial weights matrix not row standardized")
	all.tests <- c("LMerr", "LMlag", "RLMerr", "RLMlag", "SARMA")
	if (test[1] == "all") test <- all.tests
	if (!all(test %in% all.tests))
		stop("Invalid test selected - must be either \"all\" or a vector of tests")		

	y <- model.response(model.frame(model))
	X <- model.matrix(terms(model), model.frame(model))
	yhat <- as.vector(fitted(model))
	p <- model$rank
	p1 <- 1:p
	nacoefs <- which(is.na(coefficients(model)))
# fixed after looking at TOWN dummy in Boston data
	if (length(nacoefs) > 0L) X <- X[,-nacoefs]
	XtXinv <- chol2inv(model$qr$qr[p1, p1, drop = FALSE])
	sigma2 <- (t(u) %*% u) / N
	TrW <- tracew(listw)
	Wu <- lag.listw(listw, u, zero.policy)
	Wy <- lag.listw(listw, y, zero.policy)
	Wyhat <- lag.listw(listw, yhat, zero.policy)
	XtWyhat <- t(X) %*% Wyhat
	dutWu <- (t(u) %*% Wu) / sigma2
	resa <- (dutWu ^ 2) / TrW
	J <- (1/(N*sigma2)) *
		((t(Wyhat) %*% Wyhat) -
		(t(XtWyhat) %*% XtXinv %*% XtWyhat) +
		(TrW * sigma2))
	dutWy <- (t(u) %*% Wy) / sigma2
	nt <- length(test)
	if (nt < 1) stop("non-positive number of tests")
	tres <- vector(mode="list", length=nt)
	names(tres) <- test
	for (i in 1:nt) {
		testi <- test[i]
		zz <- switch(testi,
		LMerr = vec <- c(resa, 1),
		LMlag = vec <- c((dutWy ^ 2) / (N * J), 1),
		RLMerr = vec <- c(((dutWu - (TrW*((N*J)^-1))*dutWy)^2) /
			(TrW * (1 - TrW*((N*J)^-1))), 1),
		RLMlag = vec <- c(((dutWy - dutWu)^2)/ ((N*J) - TrW), 1),
		SARMA = vec <- c(((dutWy - dutWu)^2)/ ((N*J) - TrW) + resa, 2)
		)
		if (is.null(zz)) stop(paste(testi, ": no such test", sep=""))
		statistic <- vec[1]
		names(statistic) <- testi
		parameter <- vec[2]
		names(parameter) <- "df"
		p.value <- 1 - pchisq(statistic, parameter)
		if (!is.finite(p.value) || p.value < 0 || p.value > 1) 
		    warning("Out-of-range p-value: reconsider test arguments")
		names(p.value) <- ""
		method <- "Lagrange multiplier diagnostics for spatial dependence"
		data.name <- paste("\n", paste(strwrap(paste("model: ",
		    gsub("[ ]+", " ", paste(deparse(model$call), 
		    sep="", collapse="")))), collapse="\n"),
    	            "\nweights: ", listw_name, "\n", sep="")
		tres[[i]] <- list(statistic=statistic, parameter=parameter,
			p.value=p.value, method=method, data.name=data.name)
		class(tres[[i]]) <- "htest"
	}
	class(tres) <- "LMtestlist"
	tres
}

print.LMtestlist <- function(x, ...) {
	for (i in seq(along=x)) print(x[[i]])
	invisible(x)
}

summary.LMtestlist <- function(object, p.adjust.method="none", ...) {
    res <- as.data.frame(t(sapply(object, "[", 1:3)))
    res[,3] <- p.adjust(res[,3], method=p.adjust.method)
    object$results <- res
    class(object) <- "LMtestlist.summary"
    object
}

print.LMtestlist.summary <- function(x, digits=max(3, getOption("digits") - 2), ...) {
    cat(strwrap(x[[1]]$method, prefix = "\t"), sep = "\n")
    cat("data: ", x[[1]]$data.name, "\n")
    printCoefmat(x$results, has.Pvalue=TRUE, digits=digits, ...)
    invisible(x)
}

tracew <- function (listw) {
	dlmtr <- 0
	n <- length(listw$neighbours)
	if (n < 1) stop("non-positive n")
	ndij <- card(listw$neighbours)
	dlmtr <- 0
	for (i in 1:n) {
		dij <- listw$neighbours[[i]]
		wdij <- listw$weights[[i]]
		for (j in seq(length=ndij[i])) {
			k <- dij[j]
# Luc Anselin 2006-11-11 problem with asymmetric listw
			    dk <- which(listw$neighbours[[k]] == i)
			    if (length(dk) > 0L && dk > 0L &&
				dk <= length(listw$neighbours[[k]]))
				wdk <- listw$weights[[k]][dk]
			    else wdk <- 0
			    dlmtr <- dlmtr + (wdij[j]*wdij[j]) + (wdij[j]*wdk)
		}
	}
	dlmtr
}

lm.LMErr <- function(model, listw, zero.policy=FALSE, spChk=NULL) {
	if (!inherits(listw, "listw")) stop("listw is not a listw object")
	N <- length(listw$neighbours)
	if (inherits(model, "lm")) u <- resid(model)
	else if (is.numeric(model) && length(model) == N) {
		u <- model
		if (!isTRUE(all.equal(mean(u), 0.0)))
		    warning("mean of externally provided residuals not zero")
	} else stop(paste(deparse(substitute(model)),
		"not an lm object or a numeric vector of correct length"))

	if (N != length(u)) stop("objects of different length")
	if (is.null(spChk)) spChk <- get.spChkOption()
	if (spChk && !chkIDs(u, listw))
		stop("Check of data and weights ID integrity failed")
	u <- as.vector(u)

	if (is.null(attr(listw$weights, "W")) || !attr(listw$weights, "W"))
		warning("Spatial weights matrix not row standardized")
	TrW <- tracew(listw)
	Wu <- lag.listw(listw, u, zero.policy)
	sigma2 <- (t(u) %*% u) / N
	dutWu <- (t(u) %*% Wu) / sigma2
	resa <- (dutWu ^ 2) / TrW
	statistic <- resa
	names(statistic) <- "LMErr"
	parameter <- 1
	names(parameter) <- "df"
	p.value <- 1 - pchisq(statistic, parameter)
	if (!is.finite(p.value) || p.value < 0 || p.value > 1) 
	    warning("Out-of-range p-value: reconsider test arguments")
	names(p.value) <- ""
	method <- "Lagrange multiplier diagnostics for spatial dependence"
	if (inherits(model, "lm")) data.name <- paste("\n", 
	    paste(strwrap(paste("model: ",
	    gsub("[ ]+", " ", paste(deparse(model$call), 
	    sep="", collapse="")))), collapse="\n"),
    	     "\nweights: ", deparse(substitute(listw)), "\n", sep="")
        else data.name <- paste("\nresiduals: ", deparse(substitute(model)),
    	     "\nweights: ", deparse(substitute(listw)), "\n", sep="")
	res <- list(statistic=statistic, parameter=parameter,
		p.value=p.value, method=method, data.name=data.name)
	class(res) <- "htest"
	res
}
