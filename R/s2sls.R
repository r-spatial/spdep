# Copyright 2006 by Luc Anselin and Roger Bivand
# modified by Gianfranco Piras on December 11, 2009 (added the argument legacy)
# and on March 12, 2010 (added the argument W2X)
stsls <- function(formula, data = list(), listw, zero.policy=NULL,
	na.action=na.fail, robust=FALSE, HC=NULL, legacy=FALSE, W2X=TRUE) {


    	if (!inherits(listw, "listw")) 
        	stop("No neighbourhood list")

        if (is.null(zero.policy))
            zero.policy <- get("zeroPolicy", envir = .spdepOptions)
        stopifnot(is.logical(zero.policy))
        if (class(formula) != "formula") formula <- as.formula(formula)
    	mt <- terms(formula, data = data)
    	mf <- lm(formula, data, na.action=na.action, method="model.frame")
    	na.act <- attr(mf, "na.action")
    	if (!is.null(na.act)) {
        	subset <- !(1:length(listw$neighbours) %in% na.act)
        	listw <- subset(listw, subset, zero.policy=zero.policy)
    	}

    	y <- model.extract(mf, "response")
    	if (any(is.na(y))) stop("NAs in dependent variable")
    	X <- model.matrix(mt, mf)
    	if (any(is.na(X))) stop("NAs in independent variable")
        if (robust) {
            if (is.null(HC)) HC <- "HC0"
            if (!any(HC %in% c("HC0", "HC1")))
                stop("HC must be one of HC0, HC1")
        }
# modified to pass zero.policy Juan Tomas Sayago 100913
	Wy <- lag.listw(listw, y, zero.policy=zero.policy)
	dim(Wy) <- c(nrow(X),1)
	colnames(Wy) <- c("Rho")

#	WX <- lag.listw(W,X[,2:ncol(X)])
	n <- NROW(X)
	m <- NCOL(X)
	xcolnames <- colnames(X)
	K <- ifelse(xcolnames[1] == "(Intercept)", 2, 1)
	if (m > 1) {
	    WX <- matrix(nrow=n, ncol=(m-(K-1)))
	    if(W2X) WWX <- matrix(nrow = n, ncol = ncol(WX) ) 
	    for (k in K:m) {
		wx <- lag.listw(listw, X[,k], zero.policy=zero.policy)
                if(W2X) wwx <- lag.listw(listw, wx, zero.policy = zero.policy)
		if (any(is.na(wx)))
		    stop("NAs in lagged independent variable")
		WX[,(k-(K-1))] <- wx
		if(W2X) WWX[, (k - (K - 1))] <- wwx		
	    }
            if(W2X) inst <- cbind(WX, WWX)
            else inst <- WX
	}
	if (K == 2 && listw$style != "W") {
# modified to meet other styles, email from Rein Halbersma
		wx1 <- as.double(rep(1, n))
		wx <- lag.listw(listw, wx1, zero.policy=zero.policy)
		if(W2X) wwx <- lag.listw(listw, wx, zero.policy=zero.policy)
                if (m > 1) {
                    inst <- cbind(wx, inst)
                    if(W2X) inst <- cbind(wwx, inst)
		} else {
                    inst <- matrix(wx, nrow=n, ncol=1)
                    if(W2X) inst <- cbind(inst, wwx)
                }
#		colnames(inst) <- xcolnames

	}
#	if (listw$style == "W") colnames(WX) <- xcolnames[-1]
        result <- tsls(y=y, yend=Wy, X=X, Zinst=inst, robust=robust, HC=HC,
            legacy=legacy)
	result$zero.policy <- zero.policy
	result$robust <- robust
        if (robust) result$HC <- HC
	result$legacy <- legacy
        result$listw_style <- listw$style
	result$call <- match.call()
	class(result) <- "stsls"
	result
}
#	    result <- list(coefficients=biv,var=varb,s2=s2,
#	          residuals=e)

print.stsls <- function(x, ...) {
	cat("\nCall:\n")
	print(x$call)
	cat("\nCoefficients:\n")
	print(coef(x))
	cat("\n")
	invisible(x)
}

summary.stsls <- function(object, correlation = FALSE, ...) {
	rest.se <- sqrt(diag(object$var))
	object$Coef <- cbind(object$coefficients, rest.se, 
		object$coefficients/rest.se,
		2*(1-pnorm(abs(object$coefficients/rest.se))))
	if (object$robust) colnames(object$Coef) <- c("Estimate", 
		paste(object$HC, "std. Error"), "z value", "Pr(>|z|)")
	else colnames(object$Coef) <- c("Estimate", "Std. Error", 
		"t value", "Pr(>|t|)")

	rownames(object$Coef) <- names(object$coefficients)
	if (correlation) {
		object$correlation <- diag((diag(object$var))
			^(-1/2)) %*% object$var %*% 
			diag((diag(object$var))^(-1/2))
		dimnames(object$correlation) <- dimnames(object$var)
	}
	structure(object, class=c("summary.stsls", class(object)))
}

print.summary.stsls <- function(x, digits = max(5, .Options$digits - 3),
	signif.stars = FALSE, ...) {
	cat("\nCall:", deparse(x$call),	sep = "", fill=TRUE)
	cat("\nResiduals:\n")
	resid <- residuals(x)
	nam <- c("Min", "1Q", "Median", "3Q", "Max")
	rq <- if (length(dim(resid)) == 2L) 
		structure(apply(t(resid), 1, quantile), dimnames = list(nam, 
			dimnames(resid)[[2]]))
	else structure(quantile(resid), names = nam)
	print(rq, digits = digits, ...)
	if (x$zero.policy) {
		zero.regs <- attr(x, "zero.regs")
		if (!is.null(zero.regs))
			cat("Regions with no neighbours included:\n",
			zero.regs, "\n")
	}
	cat("\nCoefficients:", x$coeftitle, "\n")
	coefs <- x$Coef
	printCoefmat(coefs, signif.stars=signif.stars, digits=digits,
		na.print="NA")
    	correl <- x$correlation
    	cat("\n")
        if (x$robust && x$legacy) cat("Asymptotic robust residual variance: ")
#	if (x$legacy) cat("Asymptotic robust residual variance: ")
	else cat("Residual variance (sigma squared): ")
	cat(format(signif(x$sse/x$df, digits)), ", (sigma: ", 
		format(signif(sqrt(x$sse/x$df), digits)), ")\n", sep="")
	
    	if (!is.null(correl)) {
        	p <- NCOL(correl)
        	if (p > 1) {
            		cat("\nCorrelation of Coefficients:\n")
                	correl <- format(round(correl, 2), nsmall = 2, 
                  	digits = digits)
                	correl[!lower.tri(correl)] <- ""
                	print(correl[-1, -p, drop = FALSE], quote = FALSE)
            	}
    	}
    	cat("\n")
        invisible(x)

}

residuals.stsls <- function(object, ...) {
	if (is.null(object$na.action))
		object$residuals
	else napredict(object$na.action, object$residuals)
}

coef.stsls <- function(object, ...) object$coefficients

deviance.stsls <- function(object, ...) object$sse

# Copyright 2004 by Luc Anselin
# spatial two stage least squares
# Usage:
#    stsls(listw,y,X,robust)
# Arguments:
#    listw: spatial weights file as listw object
#    y: dependent variable as vector
#    X: explanatory variables as matrix using cbind(1,var1,...)
#    robust: flag for heteroskedastic robust estimator
# Details:
#    calls tsls with y as dependent variable, spatial lag of y
#    as endogenous, X as exogenous variables, spatial lags of
#    X as instruments and robust as specified
# Value:
# a list as returned by tsls
#   coefficients: coefficient estimates
#   se: (asymptotic) standard error of estimates
#   t:  value of asymptotic t-test statistic
#   p:  probability of t-test (tail, two-sided)
#   var: coefficient variance matrix
#   s2: residual variance (using degrees of freedom N-K)
#   residuals: observed y - predicted y, to be used in diagnostics

stsls_old <- function(W,y,X,robust=FALSE) {
	Wy <- lag.listw(W,y)
	dim(Wy) <- c(nrow(X),1)
	colnames(Wy) <- c("Rho")
	WX <- lag.listw(W,X[,2:ncol(X)])
    	result <- tsls(y,Wy,X,WX,robust)
	result
}


# Copyright 2004 by Luc Anselin
# heteroskedastic two stage least squares
# helper function, called from tsls
# Usage:
#    htsls(y,Z,Q,e)
# Arguments:
#    y: dependent variable as vector
#    Z: matrix of endogenous and exogenous variables
#    Q: matrix of instruments
#    e: vector of 2SLS residuals
# Details:
#    uses White consistent estimator for XOmegaX 
# Value:
# a list with results
#   coefficients: coefficient estimates
#   se: (asymptotic) standard error of coefficients
#   t: value of asymptotic t-test statistic
#   p: probability of t-test (tail, two-sided)
#   var: coefficient variance matrix
#   s2: residual variance (using N)
#   residuals: observed y - predicted y

htsls <- function(y,Z,Q,e) {
	e2 <- e^2
	oQ <- e2[,1] * Q
	QoQ <- crossprod(Q,oQ)
	QoQi <- solve(QoQ)
	QZ <- crossprod(Q,Z)
	ZQoQ <- crossprod(QZ,QoQi)
	v <- ZQoQ %*% QZ
	vi <- solve(v)
	Qy <- crossprod(Q,y)
	ZQy <- ZQoQ %*% Qy
	biv <- vi %*% ZQy
        yp <- Z %*% biv
    	e <- y - yp
	biv <- biv[,1,drop=TRUE]
    	sse <- c(crossprod(e,e)) # / nrow(Z)
	df <- nrow(Z)
#    	sebiv <- sqrt(diag(vi))
#    	tbiv <- biv / sebiv
#    	pbiv <- pnorm(abs(tbiv),lower.tail=FALSE) * 2
    	result <- list(coefficients=biv,
#            se=sebiv,t=tbiv,p=pbiv,
	    var=vi,sse=sse,residuals=c(e),df=df)
    	result
}


# Copyright 2004 by Luc Anselin
# two stage least squares
# Usage:
#    tsls(y,yend,X,Zinst,robust=FALSE)
# Arguments:
#    y: dependent variable as vector
#    yend: endogenous variables as vector or matrix (using cbind)
#    X: matrix of exogenous variables, including constant
#    Zinst: matrix of instruments (using cbind)
#    robust: flag for heteroskedastic robust estimator
# Details:
#    standard two stage least squares, using explicit two stages
#    uses degrees of freedom in computation of residual variance (N-K not N)
#    calls htsls when robust is TRUE
# Value:
# a list with results:
#   coefficients: coefficient estimates
#   se: (asymptotic) standard error of estimates
#   t:  value of asymptotic t-test statistic
#   p:  probability of t-test (tail, two-sided)
#   var: coefficient variance matrix
#   s2: residual variance (using degrees of freedom N-K)
#   residuals: observed y - predicted y, to be used in diagnostics

tsls <- function(y,yend,X,Zinst,robust=FALSE, HC="HC0", legacy=FALSE) {
#	colnames(X) <- c("CONSTANT",colnames(X)[2:ncol(X)])
	Q <- cbind(X,Zinst)
	Z <- cbind(yend,X)
	df <- nrow(Z) - ncol(Z)
#	QQ <- crossprod(Q,Q)
	Qye <- crossprod(Q,yend)
        Qr <- qr(Q)
        bz <- chol2inv(Qr$qr)%*% Qye
#	bz <- solve(QQ,Qye)
	yendp <- Q %*% bz
	Zp <- cbind(yendp,X)
        Qr <- qr(Zp)
#	ZpZp <- crossprod(Zp,Zp)
#	ZpZpi <- solve(ZpZp)
	ZpZpi <- chol2inv(Qr$qr)
	Zpy <- crossprod(Zp,y)
        biv <- ZpZpi %*% Zpy
#	biv <- crossprod(ZpZpi,Zpy)
	yp <- Z %*% biv
	biv <- biv[,1,drop=TRUE]
        names(biv) <- colnames(Zp)
	e <- y - yp
	if (robust) {
		if (legacy) {		
		result <- htsls(y,Z,Q,e)
		} else {
	        	sse <- c(crossprod(e,e))
                        if (HC == "HC0") omega <- as.numeric(e^2)
                        else if (HC == "HC1")
                            omega <- (nrow(X)/df) * as.numeric(e^2)
                        else stop("invalid HC choice")
			ZoZ<-crossprod(Zp,(Zp*omega))
			varb<-ZpZpi%*%ZoZ%*%ZpZpi
	   
	   		result <- list(coefficients=biv,
				var=varb,
				sse=sse,
	        		residuals=c(e),
				df=df)

		}
	} else {	
	    sse <- c(crossprod(e,e))
    	    s2 <- sse / df
	    varb <- ZpZpi * s2
#	    sebiv <- sqrt(diag(varb))
#	    tbiv <- biv / sebiv
#	    pbiv <- pnorm(abs(tbiv),lower.tail=FALSE) * 2
	    result <- list(coefficients=biv,
#		  se=sebiv,t=tbiv,p=pbiv,
		var=varb,
		sse=sse,
	        residuals=c(e),
		df=df)
	}
	result
}
