# Copyright 2023-4 by Roger Bivand 
#

is.formula <- function(x){
   inherits(x,"formula")
}

create_X0 <- function(X, listw, Durbin=TRUE, data=NULL, na.act=NULL, have_factor_preds=FALSE) {
        if (isTRUE(Durbin)) {
            if (have_factor_preds) warn_factor_preds(have_factor_preds)
            n <- NROW(X)
	    m <- NCOL(X)
	    # check if there are enough regressors
	    xcolnames <- colnames(X)
            stopifnot(!is.null(xcolnames))
	    K <- ifelse(xcolnames[1] == "(Intercept)", 2, 1)
            vars <- NULL
            xI <- NULL
            X0 <- NULL
	    if (K == 2) {
            # unnormalized weight matrices
               	if (!(listw$style == "W")) {
 		    xI <- as.double(rep(1, n))
                    vars <-"X0.(Intercept)"
               	} 
            }   
	    if (m > 1 || (m == 1 && K == 1)) {
                X0 <- matrix(as.numeric(NA), nrow=n,
                    ncol=ifelse(m==1, 1, (m-(K-1))))
		for (k in K:m) {
                        j <- ifelse(k==1, 1, k-(K-1))
			X0[,j] <- X[,xcolnames[k]]
                        vars <- c(vars, xcolnames[k])
		}
	    }
            if (!is.null(xI)) X0 <- cbind(xI, X0)
            colnames(X0) <- vars
            rownames(X0) <- rownames(X)
        } else if (is.formula(Durbin)) {
            data1 <- data
            if (!is.null(na.act) && (inherits(na.act, "omit") ||
                inherits(na.act, "exclude"))) {
                data1 <- data1[-c(na.act),]
            }
            dmf <- lm(Durbin, data1, na.action=na.fail, 
	        method="model.frame")
	    formula_durbin_factors <- have_factor_preds_mf(dmf)
            if (formula_durbin_factors) 
                warn_factor_preds(formula_durbin_factors)
#	    dmf <- lm(Durbin, data, na.action=na.action, 
#	         method="model.frame")
            X0 <- try(model.matrix(Durbin, dmf), silent=TRUE)
            if (inherits(X0, "try-error")) 
                 stop("Durbin variable mis-match")
            
            inds <- match(colnames(X0), colnames(X))
            if (anyNA(inds)) {
              wna <- which(is.na(inds)) #TR: continue if Durbin has intercept, but formula has not
              if (length(wna) == 1 && grepl("Intercept", colnames(X0)[wna])
                 && attr(terms(Durbin), "intercept") == 1) {
                inds <- inds[-wna]
              } else {
                stop("X0 variables not in X: ",
                     paste(colnames(X0)[is.na(inds)], collapse=" "))
              }
            } 
            icept <- grep("(Intercept)", colnames(X0))
            if (length(icept) == 1L && listw$style == "W") 
                X0 <- X0[, -icept, drop=FALSE]
        } else stop("Durbin argument neither TRUE nor formula")
        X0
}

SD.RStests <- function(model, listw, zero.policy=attr(listw, "zero.policy"), test="SDM", Durbin=TRUE, data=NULL) {

	if (inherits(model, "lm")) na.act <- model$na.action
	else na.act <- attr(model, "na.action")

	listw_name <- deparse(substitute(listw))

	SDM.tests <- c("SDM_RSlag", "SDM_adjRSlag", "SDM_RSWX", "SDM_adjRSWX", "SDM_Joint")
	SDEM.tests <- c("SDEM_RSerr", "SDEM_RSWX", "SDEM_Joint")
        all.tests <- c(SDM.tests, SDEM.tests)
	if (test[1] == "SDM") test <- SDM.tests
	if (test[1] == "SDEM") test <- SDEM.tests
        if (test[1] == "all") test <- all.tests
	if (!all(test %in% all.tests))
	  stop("Invalid test selected - must be either \"all\", \"SDM\", \"SDEM\" or a vector of tests")		
	nt <- length(test)
	if (nt < 1) stop("non-positive number of tests")

	if (!inherits(listw, "listw")) stop(paste(listw_name,
		"is not a listw object"))
        if (is.null(zero.policy))
            zero.policy <- get.ZeroPolicyOption()
        stopifnot(is.logical(zero.policy))
	if (!is.null(na.act)) {
	    subset <- !(1:length(listw$neighbours) %in% na.act)
	    listw <- subset(listw, subset, zero.policy=zero.policy)
	}

	if(!inherits(model, "lm")) stop(paste(deparse(substitute(model)),
		"not an lm object"))
	N <- length(listw$neighbours)
	u <- resid(model)
	if (N != length(u)) stop("objects of different length")
	u <- as.vector(u)

	if (is.null(attr(listw$weights, "W")) || !attr(listw$weights, "W"))
		warning("Spatial weights matrix not row standardized")

        if (is.formula(Durbin)) {
            if (is.null(data)) stop("Original data object from lm() call required for formula Durbin terms")
        }

	mf <- model.frame(model)
        y <- model.response(mf)
	X <- model.matrix(terms(model), mf)
        have_factor_preds <- have_factor_preds_mf(mf)
        X0 <- create_X0(X=X, listw=listw, Durbin=Durbin, data=data, na.act=na.act,
            have_factor_preds=have_factor_preds)
	yhat <- as.vector(fitted(model))
	p <- model$rank
	p1 <- 1:p
	nacoefs <- which(is.na(coefficients(model)))
# fixed after looking at TOWN dummy in Boston data
	if (length(nacoefs) > 0L) X <- X[,-nacoefs]
	XtXinv <- chol2inv(model$qr$qr[p1, p1, drop = FALSE])
	sigma2 <- c(t(u) %*% u) / N
	TrW <- tracew(listw)
	Wu <- lag.listw(listw, u, zero.policy)
	Wy <- lag.listw(listw, y, zero.policy)
        dr <- (t(Wy) %*% u)/sigma2 # lagged y
        dl <- (t(Wu) %*% u)/sigma2 # lagged residuals
	Wyhat <- lag.listw(listw, yhat, zero.policy)
        WX0 <- lag.listw(listw, X0, zero.policy)
        dg <- c(t(WX0) %*% u)/sigma2
        k <- ncol(X)
        k0 <- ncol(X0)
        J_11 <- rbind(cbind((crossprod(X)/(N*sigma2)), rep(0, k)), 
                      cbind(t(rep(0, k)), (1/(2*(sigma2^2)))))
        invJ_11 <- solve(J_11)
        Jrp <- rbind((t(X) %*% Wyhat)/(N*sigma2), t(rep(0, 1)))
        Jgb <- (t(X) %*% WX0)/(N*sigma2)
        Jgp <- rbind(Jgb, t(rep(0, k0)))
        J_12 <- cbind(Jrp, Jgp)
        Jrr <- (c(crossprod(Wyhat)) + TrW*sigma2)/(N*sigma2)
        Jgg <- crossprod(WX0)/(N*sigma2)
        Jrg <-  (t(WX0) %*% Wyhat)/(N*sigma2)
        J_22 <- rbind(cbind(Jrr, t(Jrg)), cbind(Jrg, Jgg))
        Jrg.p <- t(Jrg) - c(t(Jrp) %*% invJ_11 %*% Jgp)
        Jr.p <- Jrr - c(t(Jrp) %*% invJ_11 %*% Jrp)
        Jg.p <- Jgg - (t(Jgp) %*% invJ_11 %*% Jgp)
        invJg.p <- solve(Jg.p)
        dr_adj <- dr - (Jrg.p %*% invJg.p %*% dg)
        Jr.p_adj <- Jr.p - (Jrg.p %*% invJg.p %*% t(Jrg.p))
        dg_adj <- dg - c(dr * (1/Jr.p)) * Jrg.p
        Jg.p_adj <- Jg.p - ((1/Jr.p) * crossprod(Jrg.p))
        J.22 <- solve(J_22 - t(J_12) %*% invJ_11 %*% J_12)
        invJg.b <- solve(Jgg - t(Jgb) %*% solve(crossprod(X)/(N*sigma2)) %*%
            Jgb)
	tres <- vector(mode="list", length=nt)
	names(tres) <- test
	for (i in 1:nt) {
		testi <- test[i]
		zz <- switch(testi,
		SDM_RSlag = vec <- c((1/N) * ((dr^2) * 1/Jr.p), 1),
		SDM_adjRSlag = vec <- c((1/N)*((dr_adj^2)*(1/Jr.p_adj)), 1),
		SDM_RSWX = vec <- c((1/N) * (t(dg) %*% invJg.p %*% dg),
                    ncol(X0)),
		SDM_adjRSWX = vec <- c((1/N) * (dg_adj %*% solve(Jg.p_adj) %*% 
                    t(dg_adj)), ncol(X0)),
		SDM_Joint = vec <- c(((1/N) * (t(c(dr, dg)) %*% 
                    J.22 %*% c(dr, dg))), ncol(X0)+1),
                SDEM_RSerr = vec <- c((dl^2) / TrW, 1),
                SDEM_RSWX = vec <- c(((t(dg) %*% invJg.b %*% dg) / N),
                    ncol(X0)),
                SDEM_Joint = vec <- c(((t(dg) %*% invJg.b %*% dg) / N) + 
                    ((dl^2) / TrW), ncol(X0)+1)
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
		method <- "Rao's score test spatial Durbin diagnostics"
                Durf <- ""
                if (is.formula(Durbin))
                    Durf <- paste0("Durbin: ", paste(as.character(Durbin),
                    collapse=" "), "\n")
		data.name <- paste("\n", paste(strwrap(paste("model: ",
		    gsub("[ ]+", " ", paste(deparse(model$call), 
		    sep="", collapse="")))), collapse="\n"),
    	            "\nweights: ", listw_name, "\n", Durf, sep="")
		tres[[i]] <- list(statistic=statistic, parameter=parameter,
			p.value=p.value, method=method, data.name=data.name)
		class(tres[[i]]) <- "htest"
	}
	class(tres) <- "RStestlist"
	tres
}


