SpatialFiltering <- function (formula, lagformula, data=list(), nb,
 glist=NULL, style="C", zero.policy=NULL, tol=0.1, zerovalue = 0.0001,
 ExactEV=FALSE, symmetric=TRUE, alpha=NULL, alternative="two.sided",
 verbose=NULL) {
#
# tol: tolerance value for convergence of spatial filtering (Moran's I).
# The search for eigenvector terminates, once the residual
# autocorrelation falls below abs(Moran's I) < tol. For positive
# spatial autocorrelation in the residuals of the basic unfiltered model,
# only those eigenvectors associated with positive autocorrelation are in
# the selection set. Vice versa, for negative autocorrelation in the
# regression residuals.
#
# zerovalue: eigenvectors with eigenvalues smaller than zerovalue
# will be excluded in eigenvector search. Allows to restrict the
# search set of eigenvectors to those with extreme autocorrelation levels.
#
# ExactEV: In some incidences the approximation of using the expectation
# and variance of Moran's I from the previous iteration will lead
# to inversions. Set ExactEV=TRUE in this situation to use exact
# expectations and variances
# alpha: Added for Pedro Peres-Neto to explore its consequences as
# compared to tol= as a stopping rule.
#
#           Authors: Yongwan Chun and Michael Tiefelsdorf
#                    Dept. of Geography - The Ohio State University
#                    Columbus, Ohio 43210
#                    emails: chun.49@osu.edu and tiefelsdorf.1@osu.edu
#		Modified by Roger Bivand
#
# Reference: Tiefelsdorf M, Griffith DA. Semiparametric Filtering of Spatial
# Autocorrelation: The Eigenvector Approach. Environment and Planning A
# 2007, 39 (5) 1193 - 1221
#
#  Version 0.9.1 - September 11, 2004
# Adaptation to formula format Roger Bivand December 2005
    
    if (missing(nb)) stop("Neighbour list argument missing")
    if (missing(formula)) stop("Formula argument missing")
    if (is.null(verbose)) verbose <- get("verbose", envir = .spdepOptions)
    stopifnot(is.logical(verbose))
        if (is.null(zero.policy))
            zero.policy <- get("zeroPolicy", envir = .spdepOptions)
        stopifnot(is.logical(zero.policy))
    lw <- nb2listw(nb, glist=glist, style=style, zero.policy=zero.policy)
    if (symmetric) lw <- listw2U(lw)
    S <- listw2mat(lw)
    a <- sum(S)
    S <- nrow(S)/a*S

    nofreg <- nrow(S)           # number of observations


# Generate Eigenvectors if eigen vectors are not given
# (M1 for no SAR, MX for SAR)
    if (class(formula) != "formula") formula <- as.formula(formula)
    mt <- terms(formula, data = data)
    mf <- lm(formula, data, method="model.frame")
    y <- model.extract(mf, "response")
    if (any(is.na(y))) stop("NAs in dependent variable")
    xsar <- model.matrix(mt, mf)
    if (any(is.na(xsar))) stop("NAs in independent variable")
    if (NROW(xsar) != length(nb))
        stop("Input data and neighbourhood list have different dimensions")
    
    mx <- diag(1,nofreg) - xsar %*% qr.solve(crossprod(xsar), t(xsar))
    S <- mx %*% S %*% mx
                                                     
#Get EigenVectors and EigenValues
    eigens <- eigen(S,symmetric=symmetric)
    val <- as.matrix(eigens$values)
    vec <- as.matrix(eigens$vectors)    
    
    if (missing(lagformula)) X <- xsar
    else {
	xlag <- model.matrix(lagformula, data=data)
	isIntercept <- match("(Intercept)", colnames(xlag))
	if (!is.na(isIntercept)) xlag <- xlag[,-(isIntercept), drop=FALSE]
	X <- cbind(xsar, xlag)
    }
    coll_test <- lm(y ~ X - 1)
    if (any(is.na(coefficients(coll_test)))) stop("Collinear RHS variable detected")

    y <- as.matrix(y)
#    Xorg <- X                                 
# X will be augmented by the selected eigenvectors
    
#Total sum of squares for R2
    TSS <- sum((y - mean(y))^2)

#Compute first Moran Expectation and Variance
    nofexo <- ncol(X)
# Number of exogenous variables (incl. const)
    degfree <- nofreg - nofexo
    M <- diag(1,nofreg) - X %*% solve(crossprod(X),t(X))
    MSM <- M %*% S %*% M     
    MStat <- GetMoranStat(MSM, degfree) 
    E <- MStat$Mean
    V <- MStat$Var    

#Matrix storing the iteration history:
#   [1] Step counter of the selection procedure
#   [2] number of selected eigenvector (sorted descending)
#   [3] its associated eigenvalue
#   [4] value Moran's I for residual autocorrelation
#   [5] standardized value of Moran's I assuming a normal approximation
#   [6] p-value of [5] for given alternative
#   [7] R^2 of the model including exogenous variables and eigenvectors
#   c("Step","SelEvec","Eval","MinMi","ZMinMi","R2","gamma")
#Store the results at Step 0 (i.e., no eigenvector selected yet)
    cyMy <- crossprod(y, M) %*% y
    cyMSMy <- crossprod(y, MSM) %*% y
    IthisTime <- (cyMSMy) / (cyMy)
    zIthisTime <- (IthisTime - E) / sqrt(V)
    altfunc <- function(ZI, alternative="two.sided") {
        if (alternative == "two.sided") 
	    PrI <- 2 * pnorm(abs(ZI), lower.tail=FALSE)
        else if (alternative == "greater")
            PrI <- pnorm(ZI, lower.tail=FALSE)
        else PrI <- pnorm(ZI)
        PrI
    }

    out <- c(	0,
		0,
		0,
		IthisTime,
            	zIthisTime,
		altfunc(zIthisTime, alternative=alternative),
            	1 - ((cyMy) / TSS)
	    )
    if (verbose) cat("Step", out[1], "SelEvec", out[2], "MinMi", out[4], 
	"ZMinMi", out[5],"Pr(ZI)", out[6], "\n")
    Aout <- out
    
#Define search eigenvalue range
#The search range is restricted into a sign range based on Moran's I
#Put a sign for eigenvectors associated with their eigenvalues
#if val > zerovalue (e.g. if val > 0.0001), then 1
#if val < zerovalue (e.g. if val < -0.0001), then -1
#otherwise 0

    sel <- cbind(row(y)[,1],val,matrix(0,nofreg,1))
#    sel[,3] <- (val > zerovalue) - (val < -zerovalue)
    sel[,3] <- (val > abs(zerovalue)) - (val < -abs(zerovalue))

#Compute the Moran's I of the aspatial model (without any eigenvector)
#i.e., the sign of autocorrelation
#if MI is positive, then acsign = 1
#if MI is negative, then acsign = -1 
    
    res <- y - X %*% solve(crossprod(X), crossprod(X, y))
    acsign <- 1
    if (((crossprod(res, S) %*% res) / crossprod(res)) < 0) acsign <- -1

#If only sar model is applied or just the intercept,
#Compute and save coefficients for all eigenvectors 
    is.onlysar <- FALSE
# if (missing(xlag) & !missing(xsar))  # changed by MT
    if (missing(lagformula)) {
        is.onlysar <- TRUE
        Xcoeffs <- solve(crossprod(X), crossprod(X, y))
        gamma4eigenvec <- cbind(row(y)[,1],matrix(0,nofreg,1))
        
# Only SAR the first parameter estimation for all eigenvectors
# Due to orthogonality each coefficient can be estimate individually
        for (j in 1:nofreg) { #Loop
            if (sel[j,3] == acsign ) { #Use only feasible unselected evecs
                gamma4eigenvec[j,2] <- solve(crossprod(vec[,j]), 
                    crossprod(vec[,j], y))  
            }  
        }      
    }

# Here the actual search starts - The inner loop check each candidate -
# The outer loop selects eigenvectors until the residual autocorrelation
# falls below 'tol'
# Loop over all eigenvectors with positive or negative eigenvalue 

    oldZMinMi <- Inf
    for (i in 1:nofreg) { #Outer Loop
        z <- Inf
        idx <- 0
        
        for (j in 1:nofreg) { #Inner Loop - Find next eigenvector
            if (sel[j,3] == acsign ) { #Use only feasible unselected evecs
                xe <- cbind(X, vec[,j])  #Add test eigenvector
                
                #Based on whether it is an only SAR model or not
                if (is.onlysar) 
                    res <- y - xe %*% as.matrix(rbind(Xcoeffs,
                        gamma4eigenvec[j,2]))
                else 
                    res <- y - xe %*% solve(crossprod(xe), crossprod(xe, y))
                
                mi <- (crossprod(res, S) %*% res) / crossprod(res)

                if (ExactEV) {
                    M <- diag(1,nofreg) - xe %*% solve(crossprod(xe),t(xe))
                    degfree <- nofreg - ncol(xe)
                    MSM <- M %*% S %*% M
                    MStat <- GetMoranStat(MSM, degfree) 
                    E <- MStat$Mean
                    V <- MStat$Var          
                  }

                if (abs((mi - E) / sqrt(V)) < z) { #Identify min z(Moran)
                    MinMi = mi
                    z <- (MinMi - E) / sqrt(V)
                    idx =j
                }
            }
        }  #End inner loop
        
#Update design matrix permanently by selected eigenvector
        X <- cbind(X,vec[,idx])
        if (is.onlysar) Xcoeffs <- (rbind(Xcoeffs,gamma4eigenvec[idx,2]))        
                
        M <- diag(1,nofreg) - X %*% solve(crossprod(X),t(X))
        degfree <- nofreg - ncol(X)
        
#Update Expectation and Variance
        MSM <- M %*% S %*% M
        MStat <- GetMoranStat(MSM, degfree) 
        E <- MStat$Mean
        V <- MStat$Var          
        ZMinMi <- ((MinMi - E) / sqrt(V))
                
#Add results of i-th step
	out <- c(i, idx, val[idx],MinMi,ZMinMi, altfunc(ZMinMi, 
	    alternative=alternative), (1 - (crossprod(y, M) %*% y / TSS)))
	if (verbose) cat("Step", out[1], "SelEvec", out[2], "MinMi", 
	    out[4], "ZMinMi", out[5],"Pr(ZI)", out[6], "\n")

        Aout <- rbind(Aout, out)

#To exclude the selected eigenvector in the next loop
        sel[idx,3] <- 0 

        if (is.null(alpha)) {
	    if (abs(ZMinMi) < tol) {
		break
	    } else if (abs(ZMinMi) > abs(oldZMinMi)) {
		if (!ExactEV) {
		   cat("   An inversion has been detected. The procedure will terminate now.\n")
           	   cat("   It is suggested to use the exact expectation and variance of Moran's I\n")
           	   cat("   by setting the option ExactEV to TRUE.\n")
		}
		break
	    }
	} else {
	    if (altfunc(ZMinMi, alternative=alternative) >= alpha) break
	}
        if (!ExactEV) {
           if (abs(ZMinMi) > abs(oldZMinMi)) {
		cat("   An inversion has been detected. The procedure will terminate now.\n")
           	cat("   It is suggested to use the exact expectation and variance of Moran's I\n")
           	cat("   by setting the option ExactEV to TRUE.\n")
                break
           }
        }
        oldZMinMi <- ZMinMi
    } # End Outer Loop
    
# Regression coefficients of selected eigenvectors
    betagam <- solve(crossprod(X),crossprod(X,y))
    gammas <- as.matrix(betagam[(nofexo+1):(nrow(betagam)),1])
        
#Formatting the output
    gammas <- rbind(0, gammas)  # Add 0 for iteration zero
    out <- cbind(Aout,gammas)    
    colnames(out) <- c("Step","SelEvec","Eval","MinMi","ZMinMi","Pr(ZI)","R2","gamma")
    rownames(out) <- out[,1]
    
    selVec <- vec[,out[,2], drop=FALSE]
    colnames(selVec) <- c(paste("vec",out[2:nrow(out),2],sep=""))
    
#Generating a result object 
    SFResult <- list(selection=out, dataset=selVec)
    class(SFResult) <- "SFResult"
    return(SFResult)
}

print.SFResult <- function(x, ...) {
	print(x$selection, ...)
}

fitted.SFResult <- function(object, ...) {
	object$dataset
}

GetMoranStat <- function(MSM, degfree) {
    #MSM    : M %*% S %*% M matrix
    #         M : projection matrix
    #         S : coded symmetric spatial link matrix
    #degfree: degrees of freedom
    
    MSM <- as.matrix(MSM)
    t1 <- sum(diag(MSM))
    t2 <- sum(diag(MSM %*% MSM))
    
    E <- t1 / degfree
    V <- 2 * (degfree * t2 - t1 * t1)/(degfree * degfree * (degfree + 2))
    return(list(Mean=E,Var=V))     
}



