# translated from Matlab code sar_g.m in the Spatial Econometrics toolbox by
# James LeSage and R. Kelley Pace (http://www.spatial-econometrics.com/).
# GSoc 2011 project by Abhirup Mallik mentored by Virgilio Gómez-Rubio

spBreg_lag <- function(formula, data = list(), listw, na.action, type="lag",
    zero.policy=NULL, control=list()) {
    timings <- list()
    .ptime_start <- proc.time()
#control
    con <- list(tol.opt=.Machine$double.eps^0.5, ldet_method="SE_classic",
        Imult=2, cheb_q=5, MC_p=16L, MC_m=30L, super=NULL, spamPivot="MMD",
        in_coef=0.1, type="MC", correct=TRUE, trunc=TRUE,
        SE_method="LU", nrho=200, interpn=2000, SElndet=NULL, LU_order=FALSE,
        pre_eig=NULL, interval=c(-1, 1), ndraw=2500L, nomit=500L, thin=1L,
        verbose=FALSE, detval=NULL, prior=list(Tbeta=NULL, c_beta=NULL,
        rho=0.5, sige=1, nu=0, d0=0, a1 = 1.01, a2 = 1.01))
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
    if (type == "mixed") {
        type <- "Durbin"
        warning("type \"mixed\" deprecated, changed to \"Durbin\"")
    }
    if (type == "Durbin") {
        WX <- create_WX(x, listw, zero.policy=zero.policy, prefix="lag")
        x <- cbind(x, WX)
        rm(WX)
    } else if (type != "lag") stop("No such type:", type)
    m <- ncol(x)
    lm.base <- lm(y ~ x - 1) # doesn't like dgCMatrix
    aliased <- is.na(coefficients(lm.base))
    cn <- names(aliased)
    names(aliased) <- substr(cn, 2, nchar(cn))
    if (any(aliased)) {
        nacoef <- which(aliased)
        x <- x[,-nacoef]
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
#        nano_2 <- nano_2 + microbenchmark::get_nanotime() - nano_p
          
          ###% update rho using griddy Gibbs
#        nano_p <- microbenchmark::get_nanotime()
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

        z = epe0 - 2*rho*epe0d + rho*rho*eped
	if (idraw > 0 & idraw < nrho) 
	    ldet = detval2[idraw]
        s2 <- z/n
        ll_iter <- (ldet - ((n/2)*log(2*pi)) - (n/2)*log(s2)
            - (1/(2*s2))*z)
        psave[iter] = as.vector(rho)
        lsave[iter] <- as.vector(ll_iter)
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
    lsave <- as.vector(coda::mcmc(matrix(lsave, ncol=1), start=con$nomit+1,
        end=con$ndraw, thin=con$thin)[,1])

    timings[["finalise"]] <- proc.time() - .ptime_start
    attr(res, "timings") <- do.call("rbind", timings)
#    attr(res, "nano") <- c(nano_1, nano_2, nano_3)
    attr(res, "control") <- con
    attr(res, "type") <- type
    attr(res, "rho_out") <- rho_out
    attr(res, "listw_style") <- listw$style
    attr(res, "lsave") <- lsave
    attr(res, "ll_mean") <- as.vector(ll_mean)
    class(res) <- c("MCMC_sar_g", class(res))
    res

#output mcmc object
}


impacts.MCMC_sar_g <- function(obj, ..., tr=NULL, listw=NULL, Q=NULL) {
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
      if (iicept) {
        b1 <- beta[-icept]
      } else {
        b1 <- beta
      }
      p <- length(b1)
      if (p %% 2 != 0) stop("non-matched coefficient pairs")
      P <- cbind(b1[1:(p/2)], b1[((p/2)+1):p])
      bnames <- names(b1[1:(p/2)])
    }
    if (is.null(tr) && !is.null(listw)) n <- length(listw$neighbours)
    else if (!is.null(tr)) n <- attr(tr, "n")
    else stop("either tr or listw must be given")
    R <- nrow(samples)

    res <- intImpacts(rho=rho, beta=beta, P=P, n=n, mu=NULL, Sigma=NULL,
        irho=irho, drop2beta=drop2beta, bnames=bnames, interval=interval,
        type=type, tr=tr, R=R, listw=listw, tol=NULL,
        empirical=NULL, Q=Q, icept=icept, iicept=iicept, p=p, samples=samples)
    attr(res, "iClass") <- class(obj)
    res


}
