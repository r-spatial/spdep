# translated from Matlab code sar_g.m in the Spatial Econometrics toolbox by
# James LeSage and R. Kelley Pace (http://www.spatial-econometrics.com/).
# GSoc 2011 project by Abhirup Mallik mentored by Virgilio Gómez-Rubio

spBreg_lag <- function(formula, data = list(), listw, na.action, Durbin, type,
    zero.policy=NULL, control=list()) {
    timings <- list()
    .ptime_start <- proc.time()
#control
    con <- list(tol.opt=.Machine$double.eps^0.5, ldet_method="SE_classic",
        Imult=2, cheb_q=5, MC_p=16L, MC_m=30L, super=NULL, spamPivot="MMD",
        in_coef=0.1, type="MC", correct=TRUE, trunc=TRUE,
        SE_method="LU", nrho=200, interpn=2000, SElndet=NULL, LU_order=FALSE,
        pre_eig=NULL, interval=c(-1, 1), ndraw=2500L, nomit=500L, thin=1L,
        verbose=FALSE, detval=NULL, prior=list(rhoMH=FALSE, Tbeta=NULL,
        c_beta=NULL, rho=0.5, sige=1, nu=0, d0=0, a1 = 1.01, a2 = 1.01,
        cc = 0.2, c=NULL, T=NULL))
    priors <- con$prior
    nmsP <- names(priors)
    priors[(namp <- names(control$prior))] <- control$prior
    if (length(noNms <- namp[!namp %in% nmsP])) 
        warning("unknown names in control$prior: ",
            paste(noNms, collapse = ", "))
    control$prior <- NULL
    con$prior <- NULL
    nmsC <- names(con)
    con[(namc <- names(control))] <- control
    if (length(noNms <- namc[!namc %in% nmsC])) 
        warning("unknown names in control: ", paste(noNms, collapse = ", "))
    stopifnot(is.logical(con$verbose))
    stopifnot(is.integer(con$ndraw))
    stopifnot(is.integer(con$nomit))
    stopifnot(is.integer(con$thin))

    if (is.null(zero.policy))
        zero.policy <- get.ZeroPolicyOption()
    stopifnot(is.logical(zero.policy))
    if (class(formula) != "formula") formula <- as.formula(formula)
    mt <- terms(formula, data = data)
    mf <- lm(formula, data, na.action=na.action,  method="model.frame")
    na.act <- attr(mf, "na.action")
    if (!inherits(listw, "listw")) stop("No neighbourhood list")
    can.sim <- FALSE
    if (listw$style %in% c("W", "S")) can.sim <- can.be.simmed(listw)
    if (!is.null(na.act)) {
        subset <- !(1:length(listw$neighbours) %in% na.act)
        listw <- subset(listw, subset, zero.policy=zero.policy)
    }
    y <- model.extract(mf, "response")
#MatrixModels::model.Matrix()
#    x <- Matrix::sparse.model.matrix(mt, mf)
    x <- model.matrix(mt, mf)
    n <- nrow(x)
    if (n != length(listw$neighbours))
        stop("Input data and weights have different dimensions")
    xcolnames <- colnames(x)
    K <- ifelse(xcolnames[1] == "(Intercept)", 2, 1)
    wy <- lag.listw(listw, y, zero.policy=zero.policy)
    if (anyNA(wy)) stop("NAs in lagged dependent variable")
#create_WX
# check for dgCMatrix
    if (missing(type)) type <- "lag"
    if (type == "mixed") {
        type <- "Durbin"
    }
    if (missing(Durbin)) Durbin <- ifelse(type == "lag", FALSE, TRUE)
    if (listw$style != "W" && is.formula(Durbin)) {
        Durbin <- TRUE
        warning("formula Durbin requires row-standardised weights; set TRUE")
    }
    if (is.logical(Durbin) && isTRUE(Durbin)) type <- "Durbin"
    if (is.formula(Durbin)) type <- "Durbin"
    if (is.logical(Durbin) && !isTRUE(Durbin)) type <- "lag"

    m <- ncol(x)
    dvars <- c(NCOL(x), 0L)

#    if (type == "Durbin") {
#        WX <- create_WX(x, listw, zero.policy=zero.policy, prefix="lag")
#FIXME
    if (is.formula(Durbin) || isTRUE(Durbin)) {
        prefix <- "lag"
        if (isTRUE(Durbin)) {
            WX <- create_WX(x, listw, zero.policy=zero.policy,
                prefix=prefix)
        } else {
            dmf <- lm(Durbin, data, na.action=na.action, 
	        method="model.frame")
            fx <- try(model.matrix(Durbin, dmf), silent=TRUE)
            if (class(fx) == "try-error") 
                stop("Durbin variable mis-match")
            WX <- create_WX(fx, listw, zero.policy=zero.policy,
                prefix=prefix)
            inds <- match(substring(colnames(WX), 5,
	        nchar(colnames(WX))), colnames(x))
            if (anyNA(inds)) stop("WX variables not in X: ",
                paste(substring(colnames(WX), 5,
                nchar(colnames(WX)))[is.na(inds)], collapse=" "))
            icept <- grep("(Intercept)", colnames(x))
            iicept <- length(icept) > 0L
            if (iicept) {
                xn <- colnames(x)[-1]
            } else {
                xn <- colnames(x)
            }
            wxn <- substring(colnames(WX), nchar(prefix)+2,
                nchar(colnames(WX)))
            zero_fill <- NULL
            if (length((which(!(xn %in% wxn)))) > 0L)
                zero_fill <- length(xn) + (which(!(xn %in% wxn)))
        }
        dvars <- c(NCOL(x), NCOL(WX))
        if (is.formula(Durbin)) {
            attr(dvars, "f") <- Durbin
            attr(dvars, "inds") <- inds
            attr(dvars, "zero_fill") <- zero_fill
        }
	x <- cbind(x, WX)
	m <- NCOL(x)
	rm(WX)
    }
#        x <- cbind(x, WX)
#        rm(WX)
#    } else if (type != "lag") stop("No such type:", type)
    lm.base <- lm(y ~ x - 1) # doesn't like dgCMatrix
    aliased <- is.na(coefficients(lm.base))
    cn <- names(aliased)
    names(aliased) <- substr(cn, 2, nchar(cn))
    if (any(aliased)) {
          if (is.formula(Durbin)) {
	    stop("Aliased variables found: ",
                paste(names(aliased)[aliased], collapse=" "))
          } else {
	    warning("Aliased variables found: ",
                paste(names(aliased)[aliased], collapse=" "))
	    nacoef <- which(aliased)
		x <- x[,-nacoef]
          }
    }
    m <- ncol(x)
    timings[["set_up"]] <- proc.time() - .ptime_start
    .ptime_start <- proc.time()

    env <- new.env()
    assign("can.sim", can.sim, envir=env)
    assign("listw", listw, envir=env)
    assign("similar", FALSE, envir=env)
    assign("n", n, envir=env)
    assign("verbose", con$verbose, envir=env)
    assign("family", "SAR", envir=env)
    assign("method", con$ldet_method, envir=env)
    W <- as(listw, "CsparseMatrix")
    assign("W", W, envir=env)

    con$interval <- jacobianSetup(con$ldet_method, env, con,
        pre_eig=con$pre_eig, interval=con$interval)
    detval1 <- get("detval1", envir=env)[,1]
    detval2 <- get("detval1", envir=env)[,2]
    bprior <-  dbeta(detval1, priors$a1, priors$a2)

    nm <- paste(con$ldet_method, "set_up", sep="_")
    timings[[nm]] <- proc.time() - .ptime_start
    .ptime_start <- proc.time()

    k <- m 

    if (is.null(priors$c_beta))
        priors$c_beta <- rep(0, k)
    else 
        stopifnot(length(priors$c_beta) == k)

    if (is.null(priors$Tbeta))
        priors$Tbeta <- diag(k)*1e+12
    else
        stopifnot(nrow(priors$Tbeta) == k && ncol(priors$Tbeta) == k)

    sige <- priors$sige
    rho <- priors$rho

#% storage for draws
    bsave <- matrix(0, nrow=con$ndraw, ncol=k)
    psave <- numeric(con$ndraw)
    ssave <- numeric(con$ndraw)
    lsave <- numeric(con$ndraw)
    acc_rate <- NULL
    if (priors$rhoMH) acc_rate <- numeric(con$ndraw)

#% ====== initializations
#% compute this stuff once to save time

    TI = solve(priors$Tbeta); # see eq 5.29, Lesage & Pace (2009) p. 140
    TIc = TI%*%priors$c_beta;
           
    xpx = crossprod(x)
    xpy = crossprod(x, y)
    xpWy = crossprod(x, wy)
    nu1 = n + 2*priors$nu
    nrho = length(detval1)
    rho_out = 0
    gsize = detval1[2] - detval1[1]
    acc = 0
    noninf <- TRUE
    if (!is.null(priors$c) && !is.null(priors$T)) {
        if (length(priors$c) == 1L && is.numeric(priors$c) && 
            length(priors$T) == 1L && is.numeric(priors$T)) noninf <- FALSE
    }
    nmk = (n-k)/2
    cc <- priors$cc
#    nano_1 = 0
#    nano_2 = 0
#    nano_3 = 0
 
    timings[["complete_setup"]] <- proc.time() - .ptime_start
    .ptime_start <- proc.time()


    for (iter in 1:con$ndraw) { #% start sampling;
                  
          ##% update beta   
#        nano_p <- microbenchmark::get_nanotime()
        AI = solve((xpx + sige*TI))#,diag(rep(1,k)));        
        ys = y - rho*wy          
        b = crossprod(x, ys) + sige*TIc
        b0 = AI %*% b # see eq 5.29, p. 140
        bhat = MASS::mvrnorm(1, b0, sige*AI) #norm_rnd(sige*AI) + b0;  
        bsave[iter, 1:k] = as.vector(bhat)
#        nano_1 <- nano_1 + microbenchmark::get_nanotime() - nano_p
          
          ##% update sige
#        nano_p <- microbenchmark::get_nanotime()
        xb = x %*% bhat
        e = (ys - xb)
        d1 = 2*priors$d0 + crossprod(e)
        chi = rchisq(1, nu1) #chi = chis_rnd(1,nu1);
        sige = as.numeric(d1/chi) # see eq 5.30, p. 141
        ssave[iter] = as.vector(sige)

        if (priors$rhoMH) {

##         % metropolis step to get rho update
            i1 = max(which(detval1 <= (rho + gsize)))
	    i2 = max(which(detval1 <= (rho - gsize)))
            index = round((i1+i2)/2)
            if (!is.finite(index)) index = 1 
	    detm = detval2[index]
            if (noninf) {
                epe = (crossprod(e))/(2*sige)
            } else {
                epe = (crossprod(e))/(2*sige) + 0.5*(((rho-priors$c)^2)/(priors$T*sige))
            }
            rhox = detm - epe
            accept = 0L;
	    rho2 = rho + cc*rnorm(1);
	    while(accept == 0L) {
	        if((rho2 > con$interval[1]) & (rho2 < con$interval[2])) {
                    accept=1
	        } else { 
		    rho2 = rho + cc*rnorm(1)
	        }
            }
            i1 = max(which(detval1 <= (rho2 + gsize)))
	    i2 = max(which(detval1 <= (rho2 - gsize)))
            index = round((i1+i2)/2)
            if (!is.finite(index)) index = 1 
	    detm = detval2[index]
            yss = y - rho2*wy
            e = yss - x %*% bhat
            if (noninf) {
                epe = (crossprod(e))/(2*sige)
            } else {
                epe = (crossprod(e))/(2*sige) +
                    0.5*(((rho-priors$c)^2)/(priors$T*sige))
            }
            rhoy = detm - epe
            if ((rhoy - rhox) > exp(1)) {
                p = 1
            } else {
	        ratio = exp(rhoy-rhox)
	        p = min(1, ratio)
            }
	    ru = runif(1)
	    if(ru < p) {
  	        rho = rho2
	        acc = acc + 1
	    }
	    acc_rate[iter] = acc/iter
	    if(acc_rate[iter] < 0.4) cc=cc/1.1
	    if(acc_rate[iter] > 0.6) cc=cc*1.1	


            AI = solve((xpx + sige*TI))
            b0 = AI %*% (xpy + sige*TIc)
            bd = AI %*% (xpWy + sige*TIc)
            e0 = y - x%*%b0
            ed = wy - x%*%bd
            epe0 = as.vector(crossprod(e0))
            eped = as.vector(crossprod(ed))
            epe0d = as.vector(crossprod(ed, e0))
        } else {
          ###% update rho using griddy Gibbs
            AI = solve((xpx + sige*TI))
            b0 = AI %*% (xpy + sige*TIc)
            bd = AI %*% (xpWy + sige*TIc)
            e0 = y - x%*%b0
            ed = wy - x%*%bd
            epe0 = as.vector(crossprod(e0))
            eped = as.vector(crossprod(ed))
            epe0d = as.vector(crossprod(ed, e0))
	    nmk = (n-k)/2
	    z = epe0 - 2*detval1*epe0d + detval1*detval1*eped
	    den = detval2 - nmk*log(z)
	    den = den + bprior

	    nd = length(den)
	    adj = max(den)
	    den = den - adj
	    xd = exp(den)

	    ## trapezoid rule
	    isum = sum((detval1[2:nd] + detval1[1:(nd-1)])*(xd[2:nd]
                - xd[1:(nd-1)])/2)#VIRGILIO:FIXED
	    zd = abs(xd/isum)
	    den = cumsum(zd)

	    rnd = runif(1)*sum(zd)
            cond <- den <= rnd
            if (any(cond)) {
#	ind = which(den <= rnd)
	        idraw = which.min(cond) - 1 #max(ind)
#	    if (idraw > 0 & idraw < nrho) 
                rho = detval1[idraw]#FIXME: This sometimes fail...
            } else {
                rho_out = rho_out+1
            }
        }

        psave[iter] = as.vector(rho)
#        nano_3 <- nano_3 + microbenchmark::get_nanotime() - nano_p

    }
### % end of sampling loop
    timings[["sampling"]] <- proc.time() - .ptime_start
    .ptime_start <- proc.time()

    mat <- cbind(bsave, psave, ssave)
    colnames(mat) <- c(colnames(x), "rho", "sige")
    res <- coda::mcmc(mat, start=con$nomit+1, end=con$ndraw, thin=con$thin)
    means <- summary(res)$statistics[,1]
    beta <- means[1:(length(means)-2)]
    rho <- means[(length(means)-1)]
    xb = x %*% beta
    ys = y - rho*wy
    e = (ys - xb)
    sse <- crossprod(e)
    s2 <- sse/n
    ldet <- do_ldet(rho, env)
    ll_mean <- (ldet - ((n/2)*log(2*pi)) - (n/2)*log(s2)
        - (1/(2*s2))*sse)

    timings[["finalise"]] <- proc.time() - .ptime_start
    attr(res, "timings") <- do.call("rbind", timings)
#    attr(res, "nano") <- c(nano_1, nano_2, nano_3)
    attr(res, "control") <- con
    attr(res, "type") <- type
    attr(res, "rho_out") <- rho_out
    attr(res, "listw_style") <- listw$style
    attr(res, "ll_mean") <- as.vector(ll_mean)
    attr(res, "aliased") <- aliased
    attr(res, "acc_rate") <- acc_rate
    attr(res, "dvars") <- dvars
    attr(res, "MH") <- priors$rhoMH
    class(res) <- c("MCMC_sar_g", class(res))
    res

#output mcmc object
}


impacts.MCMC_sar_g <- function(obj, ..., tr=NULL, listw=NULL, evalues=NULL,
    Q=NULL) {
    if (is.null(listw) && !is.null(attr(obj, "listw_style")) && 
        attr(obj, "listw_style") != "W")
        stop("Only row-standardised weights supported")
    means <- summary(obj)$statistics[,1]
    irho <- length(means)-1
    drop2beta <- irho:length(means)
    rho <- means[irho]
    s2 <- means[length(means)]
    beta <- means[1:(length(means)-2)]
    icept <- grep("(Intercept)", names(beta))
    iicept <- length(icept) > 0L
    zero_fill <- NULL
    dvars <- NULL
    samples <- as.matrix(obj)
    interval <- attr(obj, "control")$interval
    if (attr(obj, "type") == "lag") {
      type <- "lag"
      if (iicept) {
        P <- matrix(beta[-icept], ncol=1)
        bnames <- names(beta[-icept])
      } else {
        P <- matrix(beta, ncol=1)
        bnames <- names(beta)
      }
      p <- length(beta)
    } else if (attr(obj, "type") == "Durbin") {
      type <- "mixed"
      if (!is.null(attr(obj, "dvars"))) {
          dvars <- attr(obj, "dvars")
          zero_fill <- attr(dvars, "zero_fill")
      }
      if (iicept) {
        b1 <- beta[-icept]
      } else {
        b1 <- beta
      }
      if (!is.null(zero_fill)) {
        if (length(zero_fill) > 0L) {
          inds <- attr(dvars, "inds")
          b1_long <- rep(0, 2*(dvars[1]-1L))
          b1_long[1:(dvars[1]-1L)] <- b1[1:(dvars[1]-1L)]
          names(b1_long)[1:(dvars[1]-1L)] <- names(b1)[1:(dvars[1]-1)]
          for (i in seq(along=inds)) {
            b1_long[(dvars[1]-1L)+(inds[i]-1L)] <- b1[(dvars[1]-1L)+i]
          }
          b1 <- b1_long
#          s_zero_fill <- sort(zero_fill, decreasing=TRUE)
#          for (i in s_zero_fill) {
#            b1 <- append(b1, values=as.numeric(NA), after=i-1L)
#          }
        }
      }
      p <- length(b1)
      if (p %% 2 != 0) stop("non-matched coefficient pairs")
      P <- cbind(b1[1:(p/2)], b1[((p/2)+1):p])
      bnames <- names(b1[1:(p/2)])
    }
    if (is.null(evalues)) {
        if (is.null(listw) && is.null(tr)) {
            stop("either tr or listw must be given")
        } else {
            if (is.null(tr) && !is.null(listw)) n <- length(listw$neighbours)
            else if (!is.null(tr)) n <- attr(tr, "n")
        }
    } else {
        if (!is.null(listw)) {
            warning("evalues given: listw will be ignored")
            listw <-NULL
        }
        if (!is.null(tr)) {
            warning("evalues given: listw will be ignored")
            tr <- NULL
        }
        n <- length(evalues)
    }

    R <- nrow(samples)

    res <- intImpacts(rho=rho, beta=beta, P=P, n=n, mu=NULL, Sigma=NULL,
        irho=irho, drop2beta=drop2beta, bnames=bnames, interval=interval,
        type=type, tr=tr, R=R, listw=listw, evalues=evalues, tol=NULL,
        empirical=NULL, Q=Q, icept=icept, iicept=iicept, p=p, samples=samples,
        zero_fill=zero_fill, dvars=dvars)
    attr(res, "iClass") <- class(obj)
    res


}


# translated from Matlab code sem_g.m in the Spatial Econometrics toolbox by
# James LeSage and R. Kelley Pace (http://www.spatial-econometrics.com/).
# GSoc 2011 project by Abhirup Mallik mentored by Virgilio Gómez-Rubio

spBreg_err <- function(formula, data = list(), listw, na.action, Durbin, etype,
    zero.policy=NULL, control=list()) {
    timings <- list()
    .ptime_start <- proc.time()
#control
    con <- list(tol.opt=.Machine$double.eps^0.5, ldet_method="SE_classic",
        Imult=2, cheb_q=5, MC_p=16L, MC_m=30L, super=NULL, spamPivot="MMD",
        in_coef=0.1, type="MC", correct=TRUE, trunc=TRUE, LAPACK=FALSE,
        SE_method="LU", nrho=200, interpn=2000, SElndet=NULL, LU_order=FALSE,
        pre_eig=NULL, interval=c(-1, 1), ndraw=2500L, nomit=500L, thin=1L,
        verbose=FALSE, detval=NULL, prior=list(lambdaMH=FALSE, Tbeta=NULL,
        c_beta=NULL, lambda=0.5, sige=1, nu=0, d0=0, a1 = 1.01, a2 = 1.01,
        cc = 0.2, gG_sige=TRUE))
    priors <- con$prior
    nmsP <- names(priors)
    priors[(namp <- names(control$prior))] <- control$prior
    if (length(noNms <- namp[!namp %in% nmsP])) 
        warning("unknown names in control$prior: ",
            paste(noNms, collapse = ", "))
    control$prior <- NULL
    con$prior <- NULL
    nmsC <- names(con)
    con[(namc <- names(control))] <- control
    if (length(noNms <- namc[!namc %in% nmsC])) 
        warning("unknown names in control: ", paste(noNms, collapse = ", "))
    stopifnot(is.logical(con$verbose))
    stopifnot(is.integer(con$ndraw))
    stopifnot(is.integer(con$nomit))
    stopifnot(is.integer(con$thin))
    stopifnot(is.logical(con$LAPACK))

    if (is.null(zero.policy))
        zero.policy <- get.ZeroPolicyOption()
    stopifnot(is.logical(zero.policy))
    if (class(formula) != "formula") formula <- as.formula(formula)
    mt <- terms(formula, data = data)
    mf <- lm(formula, data, na.action=na.action,  method="model.frame")
    na.act <- attr(mf, "na.action")
    if (!inherits(listw, "listw")) stop("No neighbourhood list")
    can.sim <- FALSE
    if (listw$style %in% c("W", "S")) can.sim <- can.be.simmed(listw)
    if (!is.null(na.act)) {
        subset <- !(1:length(listw$neighbours) %in% na.act)
        listw <- subset(listw, subset, zero.policy=zero.policy)
    }
    y <- model.extract(mf, "response")
#MatrixModels::model.Matrix()
#    x <- Matrix::sparse.model.matrix(mt, mf)
    x <- model.matrix(mt, mf)
    n <- nrow(x)
    if (n != length(listw$neighbours))
        stop("Input data and weights have different dimensions")
    xcolnames <- colnames(x)
    K <- ifelse(xcolnames[1] == "(Intercept)", 2, 1)
    wy <- lag.listw(listw, y, zero.policy=zero.policy)
    if (anyNA(wy)) stop("NAs in lagged dependent variable")
#create_WX
# check for dgCMatrix
    if (missing(etype)) etype <- "error"
    if (missing(Durbin)) Durbin <- ifelse(etype == "error", FALSE, TRUE)
    if (listw$style != "W" && is.formula(Durbin)) {
        Durbin <- TRUE
        warning("formula Durbin requires row-standardised weights; set TRUE")
    }
    if (is.logical(Durbin) && isTRUE(Durbin)) etype <- "emixed"
    if (is.formula(Durbin)) etype <- "emixed"
    if (is.logical(Durbin) && !isTRUE(Durbin)) etype <- "error"

    m <- ncol(x)
    xxcolnames <- colnames(x)
    dvars <- c(NCOL(x), 0L)

    if (is.formula(Durbin) || isTRUE(Durbin)) {
        prefix <- "lag"
        if (isTRUE(Durbin)) {
            WX <- create_WX(x, listw, zero.policy=zero.policy,
                prefix=prefix)
        } else {
            dmf <- lm(Durbin, data, na.action=na.action, 
	        method="model.frame")
            fx <- try(model.matrix(Durbin, dmf), silent=TRUE)
            if (class(fx) == "try-error") 
                stop("Durbin variable mis-match")
            WX <- create_WX(fx, listw, zero.policy=zero.policy,
                prefix=prefix)
            inds <- match(substring(colnames(WX), 5,
	        nchar(colnames(WX))), colnames(x))
            if (anyNA(inds)) stop("WX variables not in X: ",
                paste(substring(colnames(WX), 5,
                nchar(colnames(WX)))[is.na(inds)], collapse=" "))
            icept <- grep("(Intercept)", colnames(x))
            iicept <- length(icept) > 0L
            if (iicept) {
                xn <- colnames(x)[-1]
            } else {
                xn <- colnames(x)
            }
            wxn <- substring(colnames(WX), nchar(prefix)+2,
                nchar(colnames(WX)))
            zero_fill <- NULL
            if (length((which(!(xn %in% wxn)))) > 0L)
                zero_fill <- length(xn) + (which(!(xn %in% wxn)))
        }
        dvars <- c(NCOL(x), NCOL(WX))
        if (is.formula(Durbin)) {
            attr(dvars, "f") <- Durbin
            attr(dvars, "inds") <- inds
            attr(dvars, "zero_fill") <- zero_fill
        }
	x <- cbind(x, WX)
        xcolnames <- colnames(x)
	m <- NCOL(x)
	rm(WX)
    }
#        x <- cbind(x, WX)
#        rm(WX)
#    } else if (type != "lag") stop("No such type:", type)
    lm.base <- lm(y ~ x - 1) # doesn't like dgCMatrix
    aliased <- is.na(coefficients(lm.base))
    cn <- names(aliased)
    names(aliased) <- substr(cn, 2, nchar(cn))
    if (any(aliased)) {
          if (is.formula(Durbin)) {
	    stop("Aliased variables found: ",
                paste(names(aliased)[aliased], collapse=" "))
          } else {
	    warning("Aliased variables found: ",
                paste(names(aliased)[aliased], collapse=" "))
	    nacoef <- which(aliased)
		x <- x[,-nacoef]
          }
    }
    m <- ncol(x)
    timings[["set_up"]] <- proc.time() - .ptime_start
    .ptime_start <- proc.time()

    env <- new.env()
    assign("can.sim", can.sim, envir=env)
    assign("listw", listw, envir=env)
    assign("similar", FALSE, envir=env)
    assign("n", n, envir=env)
    assign("verbose", con$verbose, envir=env)
    assign("family", "SAR", envir=env)
    assign("method", con$ldet_method, envir=env)
    assign("LAPACK", con$LAPACK, envir=env)
    W <- as(listw, "CsparseMatrix")
    assign("W", W, envir=env)

    con$interval <- jacobianSetup(con$ldet_method, env, con,
        pre_eig=con$pre_eig, interval=con$interval)
    detval1 <- get("detval1", envir=env)[,1]
    detval2 <- get("detval1", envir=env)[,2]
    bprior <-  dbeta(detval1, priors$a1, priors$a2)

    nm <- paste(con$ldet_method, "set_up", sep="_")
    timings[[nm]] <- proc.time() - .ptime_start
    .ptime_start <- proc.time()

    k <- m 

    if (is.null(priors$c_beta))
        priors$c_beta <- rep(0, k)
    else 
        stopifnot(length(priors$c_beta) == k)

    if (is.null(priors$Tbeta))
        priors$Tbeta <- diag(k)*1e+12
    else
        stopifnot(nrow(priors$Tbeta) == k && ncol(priors$Tbeta) == k)

    sige <- priors$sige
    lambda <- priors$lambda
    cc <- priors$cc

#% storage for draws
    bsave <- matrix(0, nrow=con$ndraw, ncol=k)
    psave <- numeric(con$ndraw)
    ssave <- numeric(con$ndraw)
    acc_rate <- NULL
    if (priors$lambdaMH) acc_rate <- numeric(con$ndraw)
#% ====== initializations
#% compute this stuff once to save time

    TI = solve(priors$Tbeta); # see eq 5.29, Lesage & Pace (2009) p. 140
    TIc = TI%*%priors$c_beta;
           
    if (m > 1 || (m == 1 && K == 1)) {
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


#    xpx = crossprod(x)
#    xpy = crossprod(x, y)
#    xpWy = crossprod(x, wy)
    nu1 = n + 2*priors$nu
    nrho = length(detval1)
    gsize = detval1[2] - detval1[1]
    acc = 0
    nmk = (n-k)/2

# RSB Lifted out of draw_rho_sem for homoscedastic model
    if (!priors$lambdaMH) {
        rgrid = seq(con$interval[1]+0.01, con$interval[2]-0.01, 0.01)
        ng = length(rgrid)
        epet = numeric(ng)
        detxt = numeric(ng)
        for(i in 1:ng) {
            xs = x - rgrid[i]*WX
            ys = y - rgrid[i]*wy
            xsxs <- crossprod(xs)
            AI = solve(xsxs)
            bs = AI %*% crossprod(xs, ys)
            e = ys - xs%*%bs
            epet[i] = crossprod(e)
            detxt[i] = det(xsxs)
        }
        EPE = exp((spline(x=rgrid, y=log(epet), xout=detval1))$y)
        DETX = exp(spline(x=rgrid, y=log(detxt), xout=detval1)$y)
    }
 
    timings[["complete_setup"]] <- proc.time() - .ptime_start
    .ptime_start <- proc.time()

    for (iter in 1:con$ndraw) { #% start sampling;

	#% update beta   
        xss = x - lambda*WX
        AI = solve(crossprod(xss) + sige*TI)
        yss = y - lambda*wy
        b = crossprod(xss, yss) + sige*TIc
        b0 = AI %*% b
        bhat = MASS::mvrnorm(1, b0, sige*AI)
        bsave[iter, 1:k] = as.vector(bhat)

	#update sige:
	e = yss - xss %*% bhat
	ed = e - lambda*lag.listw(listw, e, zero.policy=zero.policy)
	d1 = 2*priors$d0 + crossprod(ed)
	chi = rchisq(1, nu1)
	sige = as.numeric(d1/chi)
        ssave[iter] = as.vector(sige)

        if (priors$lambdaMH) {
        #update lambda using M-H
            i1 = max(which(detval1 <= (lambda + gsize)))
	    i2 = max(which(detval1 <= (lambda - gsize)))
            index = round((i1+i2)/2)
            if (!is.finite(index)) index = 1 #Fixed this
	    detm = detval2[index]
            epe = (crossprod(e))/(2*sige)
            detx = 0.5*log(det(crossprod(xss)))
            rhox = detm - detx - epe
            accept = 0L;
	    lambda2 = lambda + cc*rnorm(1);
	    while(accept == 0L) {
	        if((lambda2 > con$interval[1]) & (lambda2 < con$interval[2])) {
                    accept=1
	        } else { 
		    lambda2 = lambda + cc*rnorm(1)
	        }
            }
            i1 = max(which(detval1 <= (lambda2 + gsize)))
	    i2 = max(which(detval1 <= (lambda2 - gsize)))
            index = round((i1+i2)/2)
            if (!is.finite(index)) index = 1 #Fixed this
	    detm = detval2[index]
            xss = x - lambda2*WX
            yss = y - lambda2*wy
            detx = 0.5*log(det(crossprod(xss)))
            e = yss - xss %*% bhat
            epe = (crossprod(e))/(2*sige)
            rhoy = detm - detx - epe
            if ((rhoy - rhox) > exp(1)) {
                p = 1
            } else {
	        ratio = exp(rhoy-rhox)
	        p = min(1, ratio)
            }
	    ru = runif(1)
	    if(ru < p) {
  	        lambda = lambda2
	        acc = acc + 1
	    }
	    acc_rate[iter] = acc/iter
	    if(acc_rate[iter] < 0.4) cc=cc/1.1
	    if(acc_rate[iter] > 0.6) cc=cc*1.1	
        } else {
# RSB updated sige added
            if (priors$gG_sige) {
                den = detval2 - 0.5*log(DETX) - nmk*log(EPE/(2*sige))
            } else {
                den = detval2 - 0.5*log(DETX) - nmk*log(EPE)
            }
            adj = max(den)
            den = den - adj
            den = exp(den)
            nden = length(den)
	
            isum = sum((detval1[2:nden] + detval1[1:(nden-1)]) * (den[2:nden] -
                den[1:(nden-1)])/2)
            z = abs(den/isum)
            den = cumsum(z)
            rnd = runif(1) * sum(z)
	    ind = which(den <= rnd)
	    idraw = max(ind)
	    if((idraw > 0) && (idraw < nrho)) {
		lambda = detval1[idraw]
	    }
        }
        psave[iter] = as.vector(lambda)
    }
### % end of sampling loop
    timings[["sampling"]] <- proc.time() - .ptime_start
    .ptime_start <- proc.time()

    mat <- cbind(bsave, psave, ssave)
    colnames(mat) <- c(colnames(x), "lambda", "sige")
    res <- coda::mcmc(mat, start=con$nomit+1, end=con$ndraw, thin=con$thin)

    emixedImps <- NULL
    if (etype == "emixed") {
        sum_stats <- summary(res)$statistics
        if (isTRUE(Durbin)) {
            odd <- (m%/%2) > 0
            if (odd) {
                m2 <- (m-1)/2
            } else {
                m2 <- m/2
            }
            if (K == 1 && odd) {
                warning("model configuration issue: no total impacts")
            } else {
                cm <- matrix(0, ncol=m, nrow=m2)
                if (K == 2) {
                    if (odd) {
                        rownames(cm) <- xxcolnames[2:(m2+1)]
                    } else {
                        rownames(cm) <- xxcolnames[1:m2]
                    }
                    for (i in 1:m2) cm[i, c(i+1, i+(m2+1))] <- 1
# drop bug fix 2016-09-21 Philipp Hunziker
                    dirImps <- sum_stats[2:(m2+1), 1:2, drop=FALSE]
                    rownames(dirImps) <- rownames(cm)
                    indirImps <- sum_stats[(m2+2):m, 1:2, drop=FALSE]
                    rownames(indirImps) <- rownames(cm)
                    tI <- res[, 2:(m2+1), drop=FALSE] +
                        res[, (m2+2):m, drop=FALSE]
                    totImps <- t(apply(tI, 2, function(x) c(mean(x), sd(x))))
                    rownames(totImps) <- rownames(cm)
                } else {
                    rownames(cm) <- xxcolnames[1:m2]
                    for (i in 1:m2) cm[i, c(i, i+m2)] <- 1
                    dirImps <- sum_stats[1:m2, 1:2, drop=FALSE]
                    rownames(dirImps) <- rownames(cm)
                    indirImps <- sum_stats[(m2+1):m, 1:2, drop=FALSE]
                    rownames(indirImps) <- rownames(cm)
                    tI <- res[, 1:m2, drop=FALSE] +
                        res[, (m2+1):m, drop=FALSE]
                    totImps <- t(apply(tI, 2, function(x) c(mean(x), sd(x))))
                    rownames(totImps) <- rownames(cm)
                    colnames(totImps) <- colnames(dirImps)
                }
            }
        } else if (is.formula(Durbin)) {
#FIXME
            m <- sum(dvars)
            m2 <- dvars[2]
            cm <- matrix(0, ncol=m, nrow=m2)
            for (i in 1:m2) {
                cm[i, c(inds[i], i+dvars[1])] <- 1
            }
            rownames(cm) <- wxn
            dirImps <- sum_stats[2:dvars[1], 1:2, drop=FALSE]
            rownames(dirImps) <- xn
            indirImps <- sum_stats[(dvars[1]+1):m, 1:2, drop=FALSE]
            if (!is.null(zero_fill)) {
                if (length(zero_fill) > 0L) {
                    lres <- vector(mode="list", length=2L)
                    for (j in 1:2) {
                        jindirImps <- rep(as.numeric(NA), (dvars[1]-1))
                        for (i in seq(along=inds)) {
                            jindirImps[(inds[i]-1)] <- indirImps[i, j]
                        }
                        lres[[j]] <- jindirImps
                    }
                    indirImps <- do.call("cbind", lres)
                }
            }
            rownames(indirImps) <- xn
            tI <- res[, inds, drop=FALSE] +
                res[, (dvars[1]+1):m, drop=FALSE]
            totImps <- t(apply(tI, 2, function(x) c(mean(x), sd(x))))
            if (!is.null(zero_fill)) {
                if (length(zero_fill) > 0L) {
                    lres <- vector(mode="list", length=2L)
                    for (j in 1:2) {
                        jtotImps <- dirImps[, j]
                        for (i in seq(along=inds)) {
                            jtotImps[(inds[i]-1)] <- totImps[i, j]
                        }
                        lres[[j]] <- jtotImps
                    }
                    totImps <- do.call("cbind", lres)
                }
            }
            rownames(totImps) <- xn
            colnames(totImps) <- colnames(dirImps)
        } else stop("undefined emixed state")
        emixedImps <- list(dirImps=dirImps, indirImps=indirImps,
            totImps=totImps)
    }


    timings[["finalise"]] <- proc.time() - .ptime_start
    attr(res, "timings") <- do.call("rbind", timings)
#    attr(res, "nano") <- c(nano_1, nano_2, nano_3)
    attr(res, "control") <- con
    attr(res, "etype") <- etype
    attr(res, "listw_style") <- listw$style
#    attr(res, "ll_mean") <- as.vector(ll_mean)
    attr(res, "aliased") <- aliased
    attr(res, "dvars") <- dvars
    attr(res, "emixedImps") <- emixedImps
    attr(res, "acc_rate") <- acc_rate
    attr(res, "cc") <- cc
    attr(res, "n") <- n
    attr(res, "k") <- k
    attr(res, "MH") <- priors$lambdaMH
    class(res) <- c("MCMC_sem_g", class(res))
    res

#output mcmc object

}

impacts.MCMC_sem_g <- function(obj, ..., tr=NULL, listw=NULL, evalues=NULL,
    Q=NULL) {
    emixedImps <- attr(obj, "emixedImps")
    if (is.null(emixedImps)) {
        stop("No indirect impacts, use summary()")
    }
    n <- attr(obj, "n")
    k <- attr(obj, "k")
    impactsWX(emixedImps, n, k, type="SDEM", method="MCMC")
}


# translated from Matlab code sac_g.m in the Spatial Econometrics toolbox by
# James LeSage and R. Kelley Pace (http://www.spatial-econometrics.com/).
# GSoc 2011 project by Abhirup Mallik mentored by Virgilio Gómez-Rubio

spBreg_sac <- function(formula, data = list(), listw, listw2=NULL, na.action, 
    Durbin, type, zero.policy=NULL, control=list()) {
    timings <- list()
    .ptime_start <- proc.time()
#control
    con <- list(tol.opt=.Machine$double.eps^0.5, ldet_method="SE_classic",
        Imult=2, cheb_q=5, MC_p=16L, MC_m=30L, super=NULL, spamPivot="MMD",
        in_coef=0.1, type="MC", correct=TRUE, trunc=TRUE,
        SE_method="LU", nrho=200, interpn=2000, SElndet=NULL, LU_order=FALSE,
        pre_eig1=NULL, pre_eig2=NULL, interval1=c(-1, 1), interval2=c(-1, 1), 
        ndraw=2500L, nomit=500L, thin=1L, verbose=FALSE, detval1=NULL,
        detval2=NULL, prior=list(Tbeta=NULL, c_beta=NULL, lambda=0.5, 
        rho=0.5, sige=1, nu=0, d0=0, a1 = 1.01, a2 = 1.01, cc1 = 0.2, 
        cc2 = 0.2))
    priors <- con$prior
    nmsP <- names(priors)
    priors[(namp <- names(control$prior))] <- control$prior
    if (length(noNms <- namp[!namp %in% nmsP])) 
        warning("unknown names in control$prior: ",
            paste(noNms, collapse = ", "))
    control$prior <- NULL
    con$prior <- NULL
    nmsC <- names(con)
    con[(namc <- names(control))] <- control
    if (length(noNms <- namc[!namc %in% nmsC])) 
        warning("unknown names in control: ", paste(noNms, collapse = ", "))
    stopifnot(is.logical(con$verbose))
    stopifnot(is.integer(con$ndraw))
    stopifnot(is.integer(con$nomit))
    stopifnot(is.integer(con$thin))

    if (is.null(zero.policy))
        zero.policy <- get.ZeroPolicyOption()
    stopifnot(is.logical(zero.policy))
    if (class(formula) != "formula") formula <- as.formula(formula)
    mt <- terms(formula, data = data)
    mf <- lm(formula, data, na.action=na.action,  method="model.frame")
    na.act <- attr(mf, "na.action")
    if (!inherits(listw, "listw")) stop("No neighbourhood list")
    if (is.null(listw2)) listw2 <- listw
    if (!is.null(con$pre_eig1) && is.null(con$pre_eig2))
        con$pre_eig2 <- con$pre_eig1
    else if (!inherits(listw2, "listw")) stop("No 2nd neighbourhood list")
    can.sim <- FALSE
    if (listw$style %in% c("W", "S")) can.sim <- can.be.simmed(listw)
    can.sim2 <- FALSE
    if (listw2$style %in% c("W", "S")) can.sim2 <- can.be.simmed(listw2)
    if (!is.null(na.act)) {
        subset <- !(1:length(listw$neighbours) %in% na.act)
        listw <- subset(listw, subset, zero.policy=zero.policy)
    }
    y <- model.extract(mf, "response")
    x <- model.matrix(mt, mf)
    n <- nrow(x)
    if (n != length(listw$neighbours))
        stop("Input data and weights have different dimensions")
    xcolnames <- colnames(x)
    K <- ifelse(xcolnames[1] == "(Intercept)", 2, 1)
    wy <- lag.listw(listw, y, zero.policy=zero.policy)
    if (anyNA(wy)) stop("NAs in lagged dependent variable")
#create_WX
# check for dgCMatrix
    if (missing(type)) type <- "sac"
    if (type == "sacmixed") {
        type <- "Durbin"
    }
    if (missing(Durbin)) Durbin <- ifelse(type == "sac", FALSE, TRUE)
    if (listw$style != "W" && is.formula(Durbin)) {
        Durbin <- TRUE
        warning("formula Durbin requires row-standardised weights; set TRUE")
    }
    if (is.logical(Durbin) && isTRUE(Durbin)) type <- "Durbin"
    if (is.formula(Durbin)) type <- "Durbin"
    if (is.logical(Durbin) && !isTRUE(Durbin)) type <- "sac"

    m <- ncol(x)
    dvars <- c(NCOL(x), 0L)

#    if (type == "Durbin") {
#        WX <- create_WX(x, listw, zero.policy=zero.policy, prefix="lag")
#FIXME
    if (is.formula(Durbin) || isTRUE(Durbin)) {
        prefix <- "lag"
        if (isTRUE(Durbin)) {
            WX <- create_WX(x, listw, zero.policy=zero.policy,
                prefix=prefix)
        } else {
            dmf <- lm(Durbin, data, na.action=na.action, 
	        method="model.frame")
            fx <- try(model.matrix(Durbin, dmf), silent=TRUE)
            if (class(fx) == "try-error") 
                stop("Durbin variable mis-match")
            WX <- create_WX(fx, listw, zero.policy=zero.policy,
                prefix=prefix)
            inds <- match(substring(colnames(WX), 5,
	        nchar(colnames(WX))), colnames(x))
            if (anyNA(inds)) stop("WX variables not in X: ",
                paste(substring(colnames(WX), 5,
                nchar(colnames(WX)))[is.na(inds)], collapse=" "))
            icept <- grep("(Intercept)", colnames(x))
            iicept <- length(icept) > 0L
            if (iicept) {
                xn <- colnames(x)[-1]
            } else {
                xn <- colnames(x)
            }
            wxn <- substring(colnames(WX), nchar(prefix)+2,
                nchar(colnames(WX)))
            zero_fill <- NULL
            if (length((which(!(xn %in% wxn)))) > 0L)
                zero_fill <- length(xn) + (which(!(xn %in% wxn)))
        }
        dvars <- c(NCOL(x), NCOL(WX))
        if (is.formula(Durbin)) {
            attr(dvars, "f") <- Durbin
            attr(dvars, "inds") <- inds
            attr(dvars, "zero_fill") <- zero_fill
        }
	x <- cbind(x, WX)
	m <- NCOL(x)
	rm(WX)
    }
    if (NROW(x) != length(listw2$neighbours))
        stop("Input data and neighbourhood list2 have different dimensions")
    w2y <- lag.listw(listw2, y, zero.policy=zero.policy)
    w2wy <- lag.listw(listw2, wy, zero.policy=zero.policy)
    lm.base <- lm(y ~ x - 1) # doesn't like dgCMatrix
    aliased <- is.na(coefficients(lm.base))
    cn <- names(aliased)
    names(aliased) <- substr(cn, 2, nchar(cn))
    if (any(aliased)) {
          if (is.formula(Durbin)) {
	    stop("Aliased variables found: ",
                paste(names(aliased)[aliased], collapse=" "))
          } else {
	    warning("Aliased variables found: ",
                paste(names(aliased)[aliased], collapse=" "))
	    nacoef <- which(aliased)
		x <- x[,-nacoef]
          }
    }
    m <- ncol(x)
    xcolnames <- colnames(x)
    K <- ifelse(xcolnames[1] == "(Intercept)", 2, 1)
    if (any(is.na(wy)))
        stop("NAs in lagged dependent variable")
    if (m > 1 || (m == 1 && K == 1)) {
        W2X <- matrix(nrow=n,ncol=(m-(K-1)))
        for (k in K:m) {
    	    wx <- lag.listw(listw2, x[,k], zero.policy=zero.policy)
	    if (any(is.na(wx)))
	        stop("NAs in lagged independent variable")
	    W2X[,(k-(K-1))] <- wx
	}
    }
    if (K == 2) {
	wx1 <- as.double(rep(1, n))
	wx <- lag.listw(listw2, wx1, zero.policy=zero.policy)
	if (m > 1) W2X <- cbind(wx, W2X)
	else W2X <- matrix(wx, nrow=n, ncol=1)
    }
    colnames(W2X) <- xcolnames
    rm(wx)

    timings[["set_up"]] <- proc.time() - .ptime_start
    .ptime_start <- proc.time()

    env <- new.env()
    assign("can.sim", can.sim, envir=env)
    assign("can.sim2", can.sim2, envir=env)
    assign("listw", listw, envir=env)
    assign("listw2", listw2, envir=env)
    assign("similar", FALSE, envir=env)
    assign("similar2", FALSE, envir=env)
    assign("n", n, envir=env)
    assign("verbose", con$verbose, envir=env)
    assign("family", "SAR", envir=env)
    assign("method", con$ldet_method, envir=env)
    W <- as(listw, "CsparseMatrix")
    assign("W", W, envir=env)

    con$interval1 <- jacobianSetup(con$ldet_method, env, con,
        pre_eig=con$pre_eig1, interval=con$interval1, which=1)
    detval11 <- get("detval1", envir=env)[,1]
    detval12 <- get("detval1", envir=env)[,2]
    con$interval2 <- jacobianSetup(con$ldet_method, env, con,
        pre_eig=con$pre_eig2, interval=con$interval2, which=2)
    detval21 <- get("detval2", envir=env)[,1]
    detval22 <- get("detval2", envir=env)[,2]

    nm <- paste(con$ldet_method, "set_up", sep="_")
    timings[[nm]] <- proc.time() - .ptime_start
    .ptime_start <- proc.time()

    k <- m 

    if (is.null(priors$c_beta))
        priors$c_beta <- rep(0, k)
    else 
        stopifnot(length(priors$c_beta) == k)

    if (is.null(priors$Tbeta))
        priors$Tbeta <- diag(k)*1e+12
    else
        stopifnot(nrow(priors$Tbeta) == k && ncol(priors$Tbeta) == k)

    sige <- priors$sige
    rho <- priors$rho
    lambda <- priors$lambda

#% storage for draws
    bsave <- matrix(0, nrow=con$ndraw, ncol=k)
    psave <- numeric(con$ndraw)
    lsave <- numeric(con$ndraw)
    ssave <- numeric(con$ndraw)
    acc_rate1 <- numeric(con$ndraw)
    acc_rate2 <- numeric(con$ndraw)

#% ====== initializations
#% compute this stuff once to save time

    TI = solve(priors$Tbeta); # see eq 5.29, Lesage & Pace (2009) p. 140
    TIc = TI%*%priors$c_beta;
           
    xpx = crossprod(x)
    xpy = crossprod(x, y)
    xpWy = crossprod(x, wy)
    nu1 = n + 2*priors$nu
    nrho = length(detval11)
    gsize1 = detval11[2] - detval11[1]
    gsize2 = detval21[2] - detval21[1]
    acc1 = 0
    acc2 = 0
    nmk = (n-k)/2
    cc1 <- priors$cc1
    cc2 <- priors$cc2
#    nano_1 = 0
#    nano_2 = 0
#    nano_3 = 0
 
    timings[["complete_setup"]] <- proc.time() - .ptime_start
    .ptime_start <- proc.time()


    for (iter in 1:con$ndraw) { #% start sampling;
                  
          ##% update beta   
        xss = x - lambda*W2X
        AI = solve(crossprod(xss) + sige*TI)
        yss = y - rho*wy -lambda*w2y + rho*lambda*w2wy
        b = crossprod(xss, yss) + sige*TIc
        b0 = AI %*% b
        bhat = MASS::mvrnorm(1, b0, sige*AI)
        bsave[iter, 1:k] = as.vector(bhat)
          
          ##% update sige
# FIXME - sources do not use bhat
        lbhat <- solve(crossprod(xss)) %*% crossprod(xss, yss)
	e = yss - xss %*% lbhat
	ed = e - lambda*lag.listw(listw2, e, zero.policy=zero.policy)
	d1 = 2*priors$d0 + crossprod(ed)
	chi = rchisq(1, nu1)
	sige = as.numeric(d1/chi)
        ssave[iter] = as.vector(sige)

        #update lambda using M-H
        i1 = max(which(detval21 <= (lambda + gsize2)))
	i2 = max(which(detval21 <= (lambda - gsize2)))
        index = round((i1+i2)/2)
        if (!is.finite(index)) index = 1 #Fixed this
	detm = detval22[index]
        epe = (crossprod(e))/(2*sige)
        rhox = detm - epe
        accept = 0L;
	lambda2 = lambda + cc2*rnorm(1);
	while(accept == 0L) {
	    if((lambda2 > con$interval2[1]) & (lambda2 < con$interval2[2])) {
                accept=1
	    } else { 
	        lambda2 = lambda + cc2*rnorm(1)
	    }
        }
        i1 = max(which(detval21 <= (lambda2 + gsize2)))
	i2 = max(which(detval21 <= (lambda2 - gsize2)))
        index = round((i1+i2)/2)
        if (!is.finite(index)) index = 1 #Fixed this
	detm = detval22[index]
        xss = x - lambda2*W2X
        yss = y - rho*wy - lambda2*w2y + rho*lambda2*w2wy
        lbhat <- solve(crossprod(xss)) %*% crossprod(xss, yss)
        e = yss - xss %*% lbhat
        epe = (crossprod(e))/(2*sige)
        rhoy = detm - epe
        if ((rhoy - rhox) > exp(1)) {
            p = 1
        } else {
	    ratio = exp(rhoy-rhox)
	    p = min(1, ratio)
        }
	ru = runif(1)
	if(ru < p) {
  	    lambda = lambda2
	    acc2 = acc2 + 1
	}
	acc_rate2[iter] = acc2/iter
	if(acc_rate2[iter] < 0.4) cc2=cc2/1.1
	if(acc_rate2[iter] > 0.6) cc2=cc2*1.1	
        lsave[iter] = as.vector(lambda)

##      % metropolis step to get rho update
        i1 = max(which(detval11 <= (rho + gsize1)))
        i2 = max(which(detval11 <= (rho - gsize1)))
        index = round((i1+i2)/2)
        if (!is.finite(index)) index = 1 
	detm = detval12[index]
        xss = x - lambda*W2X
        yss = y - rho*wy - lambda*w2y + rho*lambda*w2wy
        lbhat <- solve(crossprod(xss)) %*% crossprod(xss, yss)
        e = yss - xss %*% lbhat
        epe = (crossprod(e))/(2*sige)
        rhox = detm - epe
        accept = 0L
	rho2 = rho + cc1*rnorm(1);
	while(accept == 0L) {
	    if((rho2 > con$interval1[1]) & (rho2 < con$interval1[2])) {
                accept=1
	    } else { 
	        rho2 = rho + cc1*rnorm(1)
	    }
        }
        i1 = max(which(detval11 <= (rho2 + gsize1)))
	i2 = max(which(detval11 <= (rho2 - gsize1)))
        index = round((i1+i2)/2)
        if (!is.finite(index)) index = 1 
	detm = detval12[index]
        yss = y - rho2*wy - lambda*w2y + rho2*lambda*w2wy
        lbhat <- solve(crossprod(xss)) %*% crossprod(xss, yss)
        e = yss - xss %*% lbhat
        epe = (crossprod(e))/(2*sige)
        rhoy = detm - epe
        if ((rhoy - rhox) > exp(1)) {
            p = 1
        } else {
	    ratio = exp(rhoy-rhox)
	    p = min(1, ratio)
        }
	ru = runif(1)
	if(ru < p) {
  	    rho = rho2
	    acc1 = acc1 + 1
	}
	acc_rate1[iter] = acc1/iter
	if(acc_rate1[iter] < 0.4) cc1=cc1/1.1
	if(acc_rate1[iter] > 0.6) cc1=cc1*1.1	

        psave[iter] = as.vector(rho)

    }
### % end of sampling loop
    timings[["sampling"]] <- proc.time() - .ptime_start
    .ptime_start <- proc.time()

    mat <- cbind(bsave, psave, lsave, ssave)
    colnames(mat) <- c(colnames(x), "rho", "lambda", "sige")
    res <- coda::mcmc(mat, start=con$nomit+1, end=con$ndraw, thin=con$thin)
    means <- summary(res)$statistics[,1]
    beta <- means[1:(length(means)-3)]
    rho <- means[(length(means)-2)]
    lambda <- means[(length(means)-1)]
    xb = (x - lambda*W2X) %*% beta
    ys = y - rho*wy - lambda*w2y + rho*lambda*w2wy
    e = (ys - xb)
    sse <- crossprod(e)
    s2 <- sse/n
    ldet <- do_ldet(rho, env)
    ll_mean <- (ldet - ((n/2)*log(2*pi)) - (n/2)*log(s2)
        - (1/(2*s2))*sse)

    timings[["finalise"]] <- proc.time() - .ptime_start
    attr(res, "timings") <- do.call("rbind", timings)
    attr(res, "control") <- con
    attr(res, "type") <- type
    attr(res, "listw_style") <- listw$style
    attr(res, "ll_mean") <- as.vector(ll_mean)
    attr(res, "aliased") <- aliased
    attr(res, "acc_rate1") <- acc_rate1
    attr(res, "acc_rate2") <- acc_rate2
    attr(res, "dvars") <- dvars
    class(res) <- c("MCMC_sac_g", class(res))
    res

#output mcmc object
}


impacts.MCMC_sac_g <- function(obj, ..., tr=NULL, listw=NULL, evalues=NULL,
    Q=NULL) {
    obj_lag <- obj[, -(which(colnames(obj) == "lambda"))]
    attributes(obj_lag) <- c(attributes(obj_lag),
        attributes(obj)[5:13])
    if (attr(obj, "type") == "sac") {
        attr(obj_lag, "type") <- "lag"
    } else {
        attr(obj_lag, "type") <- "Durbin"
    }
    class(obj_lag) <- c("MCMC_sar_g", class(obj_lag))
    res <- impacts.MCMC_sar_g(obj_lag, tr=tr, listw=listw, evalues=evalues,
        Q=Q)
    if (attr(obj, "type") == "sac") {
        attr(res, "type") <- "sac"
    } else {
        attr(res, "type") <- "sacmixed"
    }
    res
}
