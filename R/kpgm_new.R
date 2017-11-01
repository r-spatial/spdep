# Copyright 2005-8 by Luc Anselin and Roger Bivand
# Kelejian-Prucha generalized moments equations
# for spatial SAR error model
# main function
# Usage:
#    GMerrorsar(formula, data = list(), listw, na.action=na.fail, zero.policy=FALSE, control=list())
# Arguments:
#    formula: standard model formula
#    data: which data frame to search for model variables
#    listw: spatial weights file as list object
#    na.action: standard value
#    zero.policy: allow no-neighbour observations if TRUE
#    control: list of control arguments to optim (such as list(trace=1))
# Details:
#    initializes with ols, calls helper function kpwuwu to build
#    the G and g matrices, calls optim unconstrained optimizer with
#    kpgm as function and plausible starting values to get estimate
#    for lambda, then finds results with spatially weighted least squares
# Value:
# an S3 "gmsar" object

GMerrorsar <- function(#W, y, X, 
	formula, data = list(), listw, na.action=na.fail, 
	zero.policy=NULL, method="nlminb", arnoldWied=FALSE, 
        control=list(), pars, scaleU=FALSE, verbose=NULL, legacy=FALSE,
        se.lambda=TRUE, returnHcov=FALSE, pWOrder=250, tol.Hcov=1.0e-10) {
#	ols <- lm(I(y) ~ I(X) - 1)
        if (is.null(verbose)) verbose <- get("verbose", envir = .spdepOptions)
        stopifnot(is.logical(verbose))
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

	if (!inherits(listw, "listw")) stop("No neighbourhood list")

	y <- model.extract(mf, "response")
	if (any(is.na(y))) stop("NAs in dependent variable")
	x <- model.matrix(mt, mf)
	if (any(is.na(x))) stop("NAs in independent variable")
	if (NROW(x) != length(listw$neighbours))
	    stop("Input data and neighbourhood list have different dimensions")

	# added aliased after trying boston with TOWN dummy
	lm.base <- lm(y ~ x - 1)
	aliased <- is.na(coefficients(lm.base))
	cn <- names(aliased)
	names(aliased) <- substr(cn, 2, nchar(cn))
	if (any(aliased)) {
		nacoef <- which(aliased)
		x <- x[,-nacoef]
	}
	ols <- lm(y ~ x - 1)
        ukp <- residuals(ols)
	vvo <- .kpwuwu(listw, ukp, zero.policy=zero.policy,
            arnoldWied=arnoldWied, X=x)
	if (missing(pars)) {
	    scorr <- c(crossprod(lag.listw(listw, ukp,
                zero.policy=zero.policy), ukp) / crossprod(ukp, ukp))
            scorr <- scorr / (sum(unlist(listw$weights)) / length(ukp))
            if (scaleU) ukp <- scale(ukp)
            pars <- c(scorr, var(ukp))
        }
        if (length(pars) !=2L || !is.numeric(pars))
            stop("invalid starting parameter values")
	vv <- .kpwuwu(listw, ukp, zero.policy=zero.policy,
            arnoldWied=arnoldWied, X=x)
#	nlsres <- nlm(.kpgm, pars, print.level=print.level, gradtol=gradtol, steptol=steptol, iterlim=iterlim, v=vv, verbose=verbose)
#	lambda <- nlsres$estimate[1]
        if (method == "nlminb")
            optres <- nlminb(pars, .kpgm, v=vv, verbose=verbose,
               control=control)
        else 
	    optres <- optim(pars, .kpgm, v=vv, verbose=verbose,
                method=method, control=control)
        if (optres$convergence != 0)
            warning(paste("convergence failure:", optres$message))
	lambda <- optres$par[1]
	names(lambda) <- "lambda"
        GMs2 <- optres$par[2]

#        Gn <- vv$bigG
#        Gn2 <- vv$litg
#        pars <- c(lambda, lambda^2, GMs2)
#        Hfun <- function(pars, Gn, Gn2) {
#            val <- Gn2 - Gn %*% pars
#            sum(val^2)
#        }
#        e1 <- Gn2 - Gn %*% pars
#        vare1 <- sd(e1)^2

#        Hess <- fdHess(pars=pars, fun=Hfun, Gn=Gn, Gn2=Gn2)$Hessian
#        res <- solve(Hess)
#        lambda.se <- sqrt(vare1*diag(res))[1]

       lambda.se <- NULL

	wy <- lag.listw(listw, y, zero.policy=zero.policy)
	if (any(is.na(wy)))
	    stop("NAs in lagged dependent variable")
	n <- NROW(x)
	m <- NCOL(x)
	xcolnames <- colnames(x)
	K <- ifelse(xcolnames[1] == "(Intercept)", 2, 1)
	if (any(is.na(wy)))
	    stop("NAs in lagged dependent variable")
	if (m > 1) {
	    WX <- matrix(nrow=n,ncol=(m-(K-1)))
	    for (k in K:m) {
		wx <- lag.listw(listw, x[,k], zero.policy=zero.policy)
		if (any(is.na(wx)))
		    stop("NAs in lagged independent variable")
		WX[,(k-(K-1))] <- wx
	    }
	}
	if (K == 2) {
# modified to meet other styles, email from Rein Halbersma
		wx1 <- as.double(rep(1, n))
		wx <- lag.listw(listw, wx1, zero.policy=zero.policy)
		if (m > 1) WX <- cbind(wx, WX)
		else WX <- matrix(wx, nrow=n, ncol=1)
	}
	colnames(WX) <- xcolnames
	rm(wx)
	lm.target <- lm(I(y - lambda*wy) ~ I(x - lambda*WX) - 1)
	coef.lambda <- coefficients(lm.target)
	names(coef.lambda) <- xcolnames
        if (legacy) {
	    SSE <- deviance(lm.target)
	    s2 <- SSE/n
	    p <- lm.target$rank
	    rest.se <- (summary(lm.target)$coefficients[,2])*sqrt((n-p)/n)
	    r <- as.vector(residuals(lm.target))
	    fit <- as.vector(y - r)
            vcov <- vcov(lm.target)
        } else {
            fit <- as.vector(x %*% coef.lambda)
            r <- as.vector(y - fit)
            e <- residuals(ols)
            et <- e - lambda*lag.listw(listw, e, zero.policy=zero.policy)
            SSE <- c(crossprod(et))
            s2 <- SSE/n
            Bx <- x - lambda*WX
            Qr <- qr(Bx/(sqrt(s2)))
            invxpx <- chol2inv(Qr$qr)
            rest.se <- sqrt(diag(invxpx))
            vcov <- invxpx
        }

        W <- as(listw, "CsparseMatrix")
        lambda.se <- NULL
        if (!arnoldWied && se.lambda) {
# produce an std for "rho" following Kelejian-Prucha (2004)
# implemented following sem_gmm.m in the Matlab Spatial Econometrics
# toolbox, written by Shawn Bucholtz, modified extensively by J.P. LeSage
# after http://econweb.umd.edu/~prucha/STATPROG/OLS/desols.pdf
          
          KP04a <- (1/n) * vvo$trwpw
          KP04c <- sqrt(1/(1+(KP04a*KP04a)))
          KP04se <- vvo$wu
          KP04de <- vvo$wwu
          KP04eo <- residuals(ols)

          J <- matrix(0.0, ncol=2, nrow=2)
          J[1,1] <- 2*KP04c*(crossprod(KP04de, KP04se) - 
            KP04a*crossprod(KP04se, KP04eo))
          J[2,1] <- crossprod(KP04de, KP04eo) + crossprod(KP04se)
          J[1,2] <- - KP04c*(crossprod(KP04de) - KP04a*crossprod(KP04se))
          J[2,2] <- - crossprod(KP04de, KP04se)

          J <- (1/n)*J

          J1 <- J %*% matrix(c(1, 2*lambda), ncol=1)
          A2N <- crossprod(W)
          A1N <- KP04c*(A2N - KP04a*as_dsCMatrix_I(n))
          A1NA1Np <- A1N+t(A1N)
          A2NA2Np <- A2N+t(A2N)

          trA1A1 <- sum(colSums(t(A1NA1Np)*A1NA1Np))
          trA1A2 <- sum(colSums(crossprod(A2NA2Np, A1NA1Np)))
          trA2A2 <- sum(colSums(crossprod(A2NA2Np, A2NA2Np)))
          sigh <- s2*s2

          phihat <- matrix(0.0, ncol=2, nrow=2)
          phihat[1,1] <- (sigh)*trA1A1/(2*n)
          phihat[1,2] <- (sigh)*trA1A2/(2*n)
          phihat[2,1] <- (sigh)*trA1A2/(2*n)
          phihat[2,2] <- (sigh)*trA2A2/(2*n)

          JJI <- 1/crossprod(J1)
          omega <- JJI * t(J1) %*% phihat %*% J1 * JJI
          lambda.se <- sqrt(omega/n)
        }

	call <- match.call()
	names(r) <- names(y)
	names(fit) <- names(y)
        Hcov <- NULL
        if (returnHcov) {
	    pp <- ols$rank
            p1 <- 1L:pp
            R <- chol2inv(ols$qr$qr[p1, p1, drop = FALSE])
            B <- tcrossprod(R, x)
            B <- as(powerWeights(W=W, rho=lambda, order=pWOrder,
                X=B, tol=tol.Hcov), "matrix")
            C <- x %*% R
            C <- as(powerWeights(W=t(W), rho=lambda, order=pWOrder,
                X=C, tol=tol.Hcov), "matrix")
            Hcov <- B %*% C
            attr(Hcov, "method") <- "Matrix"
        }

	ret <- structure(list(type= "ERROR", lambda=lambda,
		coefficients=coef.lambda, rest.se=rest.se, 
		s2=s2, SSE=SSE, parameters=(m+2), lm.model=ols, 
		call=call, residuals=r, lm.target=lm.target,
		fitted.values=fit, formula=formula, aliased=aliased,
		zero.policy=zero.policy, vv=vv, optres=optres,
                pars=pars, Hcov=Hcov, legacy=legacy, lambda.se=lambda.se,
                arnoldWied=arnoldWied, GMs2=GMs2, scaleU=scaleU, vcov=vcov),
                class=c("gmsar"))
	if (zero.policy) {
		zero.regs <- attr(listw$neighbours, 
			"region.id")[which(card(listw$neighbours) == 0)]
		if (length(zero.regs) > 0L)
			attr(ret, "zero.regs") <- zero.regs
	}
	if (!is.null(na.act))
		ret$na.action <- na.act
	ret
}

# Copyright 2005 by Roger Bivand

residuals.gmsar <- function(object, ...) {
	if (is.null(object$na.action))
		object$residuals
	else napredict(object$na.action, object$residuals)
}

deviance.gmsar <- function(object, ...) {
	deviance(object$lm.target)
}


coef.gmsar <- function(object, ...) {
	ret <- c(object$coefficients, object$lambda)
	ret
}

fitted.gmsar <- function(object, ...) {
	if (is.null(object$na.action))
		object$fitted.values
	else napredict(object$na.action, object$fitted.values)
}


print.gmsar <- function(x, ...)
{
	cat("\nCall:\n")
	print(x$call)
	cat("\n")
	cat("\nCoefficients:\n")
	print(coef(x))
	invisible(x)
}

summary.gmsar <- function(object, correlation = FALSE, Hausman=FALSE, ...)
{
	object$coeftitle <- "(GM standard errors)"
	object$Coef <- cbind(object$coefficients, object$rest.se, 
		object$coefficients/object$rest.se,
		2*(1-pnorm(abs(object$coefficients/object$rest.se))))
	colnames(object$Coef) <- c("Estimate", "Std. Error", 
		"z value", "Pr(>|z|)")
	rownames(object$Coef) <- names(object$coefficients)
        if (Hausman && !is.null(object$Hcov)) {
                object$Haus <- Hausman.test(object)
        }

	structure(object, class=c("summary.gmsar", class(object)))
}





###modified to acomodate the SARAR model
print.summary.gmsar<-function (x, digits = max(5, .Options$digits - 3), signif.stars = FALSE, 
    ...) 
{
    cat("\nCall:", deparse(x$call), sep = "", fill = TRUE)
    cat("\nResiduals:\n")
    resid <- residuals(x)
    nam <- c("Min", "1Q", "Median", "3Q", "Max")
    rq <- if (length(dim(resid)) == 2L) 
        structure(apply(t(resid), 1, quantile), dimnames = list(nam, 
            dimnames(resid)[[2]]))
    else structure(quantile(resid), names = nam)
    print(rq, digits = digits, ...)

    if(x$type=="SARAR") cat("\nType: GM SARAR estimator")
    else  cat("\nType: GM SAR estimator")
    if (x$arnoldWied) cat(" (Arnold and Wied (2010) moment definitions)\n")
    else cat("\n")
    if (x$zero.policy) {
        zero.regs <- attr(x, "zero.regs")
        if (!is.null(zero.regs)) 
            cat("Regions with no neighbours included:\n", zero.regs, 
                "\n")
    }
    cat("Coefficients:", x$coeftitle, "\n")
    coefs <- x$Coef
    if (!is.null(aliased <- x$aliased) && any(x$aliased)) {
        cat("    (", table(aliased)["TRUE"], " not defined because of singularities)\n", 
            sep = "")
        cn <- names(aliased)
        coefs <- matrix(NA, length(aliased), 4, dimnames = list(cn, 
            colnames(x$Coef)))
        coefs[!aliased, ] <- x$Coef
    }
    printCoefmat(coefs, signif.stars = signif.stars, digits = digits, 
        na.print = "NA")
    cat("\nLambda:", format(signif(x$lambda, digits)))
    if (!is.null(x$lambda.se)) {
      cat(" (standard error):", format(signif(x$lambda.se, digits)))
      cat(" (z-value):", format(signif(x$lambda/x$lambda.se, digits)))
    }
    cat("\n")
    cat("Residual variance (sigma squared): ", format(signif(x$s2, 
        digits)), ", (sigma: ", format(signif(sqrt(x$s2), digits)), 
        ")\n", sep = "")
    if (x$scaleU) cat("(scaled) ")
    cat("GM argmin sigma squared: ", format(signif(x$GMs2, 
        digits)), "\n", sep = "")
    cat("Number of observations:", length(x$residuals), "\n")
    cat("Number of parameters estimated:", x$parameters, "\n")
    if (!is.null(x$Haus)) {
        cat("Hausman test: ", format(signif(x$Haus$statistic, 
            digits)), ", df: ", format(x$Haus$parameter), ", p-value: ", 
            format.pval(x$Haus$p.value, digits), "\n", sep = "")
    }
    cat("\n")
    invisible(x)
}

# Copyright 2004 by Luc Anselin
# Kelejian-Prucha generalized moments equations
# helper function to provide function to nonlinear optimizer
# must have parameter vector first for nlm
# Usage:
#    kpgm(par,v)
# Arguments:
#    par: 2x1 parameter vector rho,sig2
#    v: list containing bigG and litg as computed by kpwuwu
# Details:
#    sets up the equation as squared residuals
# Value:
#    value: evaluated nonlinear least squares for parameter value

.kpgm <- function(rhopar,v,verbose=FALSE) {
  vv <- v$bigG %*% c(rhopar[1],rhopar[1]^2,rhopar[2]) - v$litg
  value <- sum(vv^2)
  if (verbose)
    cat("function:", value, "lambda:", rhopar[1], "sig2:", rhopar[2], "\n")
  value
  
}


# Copyright 2004 by Luc Anselin
# Kelejian-Prucha generalized moments equations
# helper function
# Usage:
#    kpwuwu(listw,u)
# Arguments:
#    listw: spatial weights file as listw object
#    u: OLS residual vector
#    zero.policy: allow no-neighbour observations if TRUE
# Details:
#    sets up the bigG matrix and littleg vector needed
#    for the nonlinear least squares in the GM estimator
#    see Kelejian-Prucha(1999) p. 515
# Value:
# a list with two elements
#    bigG: the 3x3 G matrix
#    litg: the 3x1 g vector

.kpwuwu <- function(listw, u, zero.policy=FALSE, arnoldWied=FALSE, X=NULL) {
        if (arnoldWied) {
            stopifnot(!is.null(X))
            invXtX <- chol2inv(qr.R(qr(X)))
            W <- as(listw, "CsparseMatrix")
            WX <- W %*% X
        }
	n <- length(u)
# Gianfranco Piras 081119 
        trwpw <- sum(unlist(listw$weights)^2)
#	tt <- matrix(0,n,1)
#	for (i in 1:n) {tt[i] <- sum(W$weights[[i]]^2) }
#	trwpw <- sum(tt)
	wu <- lag.listw(listw, u, zero.policy=zero.policy)
	wwu <- lag.listw(listw, wu, zero.policy=zero.policy)
    	uu <- crossprod(u,u)
    	uwu <- crossprod(u,wu)
    	uwpuw <- crossprod(wu,wu)
    	uwwu <- crossprod(u,wwu)
    	wwupwu <- crossprod(wwu,wu)
    	wwupwwu <- crossprod(wwu,wwu)
    	bigG <- matrix(0,3,3)
        if (arnoldWied) {
            k <- ncol(X)
            uwpX <- crossprod(wu, X)
            upWX <- crossprod(u, WX)
            uwpWX <- crossprod(wu, WX)
            iXpXXpwu <- invXtX %*% t(uwpX)
            c22 <- wwu - WX %*% iXpXXpwu
            XiXpX <- X %*% invXtX
            WXpW <- t(WX) %*% W
            c23 <- trwpw - sum(diag(WXpW %*% XiXpX))
            c32 <- crossprod(wu, wwu) - t(wu) %*% WX %*% iXpXXpwu
            c32 <- c32 - ((t(wu) %*% XiXpX %*% crossprod(X, wwu)) - 
                (t(wu) %*% XiXpX %*% crossprod(X, WX) %*% iXpXXpwu))
            bigG[,1] <- c(2*uwu, 2*as.vector(wwupwu - (uwpWX %*% iXpXXpwu)),
                as.vector(uwwu - (upWX %*% iXpXXpwu)) +
                (uwpuw - (uwpX %*% iXpXXpwu)))/n
    	    bigG[,2] <- - c(as.vector(uwpuw - (uwpX %*% iXpXXpwu)), as.vector(crossprod(c22)), as.vector(c32))/n
    	    bigG[,3] <- c(n-k, as.vector(c23), -as.vector(sum(diag(t(X) %*% (WX %*% invXtX)))))/n
# M <- diag(length(u)) - X %*% invXtX %*% t(X)
# BGc1 <- c(2*as.vector(crossprod(u, W %*% u)), 2*as.vector(t(u) %*% t(W) %*% W %*% M %*% W %*% u), as.vector(u %*% (W+t(W)) %*% M %*% W %*% u))/n
# BGc2 <- - c(as.vector(t(u) %*% t(W) %*% M %*% W %*% u), as.vector(t(u) %*% t(W) %*% M %*% t(W) %*% W %*% M %*% W %*% u), as.vector(t(u) %*% t(W) %*% M %*% W %*% M %*% W %*% u))/n
# BGc3 <- c((n-k), sum(diag(M %*% t(W) %*% W)), sum(diag(W %*% M)))/n
        } else {
    	    bigG[,1] <- c(2*uwu, 2*wwupwu, (uwwu+uwpuw))/n
    	    bigG[,2] <- - c(uwpuw, wwupwwu, wwupwu) / n
    	    bigG[,3] <- c(1, trwpw/n, 0)
        }
    	litg <- c(uu,uwpuw,uwu) / n
    	list(bigG=bigG, litg=litg, trwpw=trwpw, wu=wu, wwu=wwu)
}

####SARAR model

gstsls<-function (formula, data = list(), listw, listw2=NULL,
 na.action = na.fail, zero.policy = NULL, pars, scaleU=FALSE,
 control = list(), verbose = NULL, method = "nlminb", robust = FALSE,
 legacy = FALSE, W2X = TRUE ) 
{
	
	
	 if (is.null(verbose)) 
        verbose <- get("verbose", envir = .spdepOptions)
    stopifnot(is.logical(verbose))

    if (is.null(zero.policy)) 
        zero.policy <- get.ZeroPolicyOption()
    stopifnot(is.logical(zero.policy))

    if (!inherits(listw, "listw")) 
        stop("The weights matrix is not a listw object")

    if (is.null(listw2)) 
        listw2 <- listw
    else if (!inherits(listw2, "listw")) 
        stop("No 2nd neighbourhood list")

    if (class(formula) != "formula") formula <- as.formula(formula)
    mt <- terms(formula, data = data)
    mf <- lm(formula, data, na.action = na.fail, method = "model.frame")
    na.act <- attr(mf, "na.action")
    cl <- match.call()
    if (!is.null(na.act)) {
        subset <- !(1:length(listw$neighbours) %in% na.act)
        subset2 <- !(1:length(listw2$neighbours) %in% na.act)
        listw <- subset(listw, subset, zero.policy = zero.policy)
        listw2 <- subset(listw2, subset2, zero.policy = zero.policy)    
    }
    
    

    y <- model.extract(mf, "response")
    x <- model.matrix(mt, mf)
    if (length(y) != nrow(x)) 
        stop("x and y have different length")
    if (nrow(x) != length(listw$neighbours)) 
        stop("Input data and weights have different dimension")
    if (any(is.na(y))) 
        stop("NAs in dependent variable")
    if (any(is.na(x))) 
        stop("NAs in independent variable")
    n <- nrow(x)
    k <- ncol(x)
    xcolnames <- colnames(x)
    K <- ifelse(xcolnames[1] == "(Intercept)" || all(x[, 1] == 
        1), 2, 1)


        wy <- lag.listw(listw, y, zero.policy = zero.policy)
        wy <- array(wy, c(length(y), 1L))
        colnames(wy) <- ("Wy")
        if (any(is.na(wy))) 
            stop("NAs in spatially lagged dependent variable")
        if (k > 1) {
            WX <- matrix(nrow = n, ncol = (k - (K - 1)))
            WWX <- matrix(nrow = n, ncol = (k - (K - 1)))
            for (i in K:k) {
                wx <- lag.listw(listw, x[, i], zero.policy = zero.policy)
                wwx <- lag.listw(listw, wx, zero.policy = zero.policy)
                if (any(is.na(wx))) 
                  stop("NAs in lagged independent variable")
                WX[, (i - (K - 1))] <- wx
                WWX[, (i - (K - 1))] <- wwx
            }
        }

        instr <- cbind(WX, WWX)
        firststep <- tsls(y = y, yend = wy, X = x, Zinst = instr, robust = robust, legacy = legacy)

        ukp <- residuals(firststep)

    if (missing(pars)) {
    	
        scorr <- c(crossprod(lag.listw(listw2, ukp,
            zero.policy=zero.policy), ukp)/crossprod(ukp, ukp))
        scorr <- scorr/(sum(unlist(listw2$weights))/length(ukp))
        if (scaleU) ukp <- scale(ukp)
        pars <- c(scorr, var(ukp))
    }
    if (length(pars) != 2L || !is.numeric(pars)) 
        stop("invalid starting parameter values")
    vv <- .kpwuwu(listw2, ukp, zero.policy = zero.policy,
        arnoldWied=FALSE, X=x)
    if (method == "nlminb") 
        optres <- nlminb(pars, .kpgm, v = vv, verbose = verbose, 
            control = control)
    else optres <- optim(pars, .kpgm, v = vv, verbose = verbose, 
        method = method, control = control)
    if (optres$convergence != 0) 
        warning(paste("convergence failure:", optres$message))
    lambda <- optres$par[1]
    names(lambda) <- "lambda"
    GMs2 <- optres$par[2]

#        Gn <- vv$bigG
#        Gn2 <- vv$litg
#        pars <- c(lambda, lambda^2, GMs2)
#        Hfun <- function(pars, Gn, Gn2) {
#            val <- Gn2 - Gn %*% pars
#            sum(val^2)
#        }
#        e1 <- Gn2 - Gn %*% pars
#        vare1 <- sd(e1)^2

#        Hess <- fdHess(pars=pars, fun=Hfun, Gn=Gn, Gn2=Gn2)$Hessian
#        res <- solve(Hess)
#        lambda.se <- sqrt(vare1*diag(res))[1]

       lambda.se <- NULL

        w2y <- lag.listw(listw2, y)
        yt <- y - lambda * w2y
        xt <- x - lambda * lag.listw(listw2, x)
        wyt <- wy - lambda * lag.listw(listw2, wy)

        colnames(xt) <- xcolnames
        colnames(wyt) <- c("Rho_Wy")
        secstep <- tsls(y = yt, yend = wyt, X = xt, Zinst = instr,
            robust = robust, legacy = legacy)
	rho<-secstep$coefficients[1]
	coef.sac<-secstep$coefficients
	rest.se <- sqrt(diag(secstep$var))
	rho.se <- sqrt(diag(secstep$var))[1]
	s2<-secstep$sse / secstep$df
	r<- secstep$residuals
	fit<- y - r
	SSE<- crossprod(r)

    	call <- match.call()

    	ret <- structure(list(type= "SARAR", lambda = lambda,
	coefficients = coef.sac, 
        rest.se = rest.se, s2 = s2, SSE = SSE, parameters = (k + 
            3), lm.model = NULL, call = call, residuals = r, lm.target = NULL, 
            fitted.values = fit, formula = formula, aliased = NULL, 
            zero.policy = zero.policy, vv = vv, optres = optres, 
            pars = pars, Hcov = NULL, lambda.se=lambda.se,
            arnoldWied=FALSE, GMs2=GMs2, scaleU=scaleU,
            secstep_var=secstep$var), class = c("gmsar"))
        if (zero.policy) {
        zero.regs <- attr(listw$neighbours,
            "region.id")[which(card(listw$neighbours) == 0)]
        if (length(zero.regs) > 0L) 
            attr(ret, "zero.regs") <- zero.regs
    }
        
    if (!is.null(na.act)) ret$na.action <- na.act
    ret
}

GMargminImage <- function(obj, lambdaseq, s2seq) {
    if (missing(lambdaseq)) {
        lamin <- obj$lambda-0.5
        lamin <- ifelse(lamin < -1, -1, lamin)
        lamax <- obj$lambda+0.5
        lamax <- ifelse(lamax >= 1, (1-.Machine$double.eps), lamax)
        lambdaseq <- seq(lamin, lamax, length.out=40)
    }
    if (missing(s2seq)) 
        s2seq <- seq(0.5*obj$GMs2, 1.5*obj$GMs2, length.out=40)
    xy <- as.matrix(expand.grid(lambdaseq, s2seq))
    vres <- apply(xy, 1, function(x) .kpgm(rhopar=x, v=obj$vv))
    res <- matrix(vres, ncol=length(lambdaseq))
    list(x=lambdaseq, y=s2seq, z=res)
}

