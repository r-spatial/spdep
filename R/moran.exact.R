# Copyright (c) 2007-2008 Markus Reder and Roger Bivand

lm.morantest.exact <- function(model, listw, zero.policy = NULL, 
    alternative = "greater", spChk=NULL, resfun=weighted.residuals, 
    zero.tol=1.0e-7, Omega=NULL, save.M=NULL, save.U=NULL, useTP=FALSE,
    truncErr=1e-6, zeroTreat=0.1) 
{
    if (!inherits(listw, "listw")) 
        stop(paste(deparse(substitute(listw)), "is not a listw object"))
    if (!inherits(model, "lm")) 
        stop(paste(deparse(substitute(model)), "not an lm object"))
        if (is.null(zero.policy))
            zero.policy <- get("zeroPolicy", envir = .spdepOptions)
        stopifnot(is.logical(zero.policy))
    N <- length(listw$neighbours)
    u <- resfun(model)
    if (N != length(u)) 
        stop("objects of different length")
    if (is.null(spChk)) spChk <- get.spChkOption()
    if (spChk && !chkIDs(u, listw))
	stop("Check of data and weights ID integrity failed")
    if (!(alternative %in% c("greater", "less", "two.sided")))
	stop("alternative must be one of: \"greater\", \"less\", or \"two.sided\"")
    u <- as.vector(u)
    listw.U <- listw2U(listw)
    S0 <- sum(unlist(listw.U$weights))
    lu <- lag.listw(listw.U, u, zero.policy = zero.policy)
    Nnn <- N
    if (zero.policy) Nnn <- length(which(card(listw$neighbours) > 0))
    I <- (Nnn/S0) * (crossprod(u, lu) / crossprod(u))
    I_save <- I
    if (!isTRUE(all.equal((Nnn/S0), 1))) I <- I * (S0/Nnn)
    p <- model$rank
    p1 <- 1:p
    nacoefs <- which(is.na(coefficients(model)))
    XtXinv <- chol2inv(model$qr$qr[p1, p1, drop = FALSE])
    X <- model.matrix(terms(model), model.frame(model))
# fixed after looking at TOWN dummy in Boston data
    if (length(nacoefs) > 0L) X <- X[,-nacoefs]
    if (!is.null(wts <- weights(model))) {
	X <- sqrt(diag(wts)) %*% X
    }
    M <- diag(N) - X %*% tcrossprod(XtXinv, X)
    U <- listw2mat(listw.U)
    if (is.null(Omega)) {
        MVM <- M %*% U %*% M
        MVM <- 0.5 * (t(MVM) + MVM)
        evalue <- eigen(MVM, only.values=TRUE)$values
        idxpos <- which(abs(evalue) < zero.tol)
        if (length(idxpos) != p)
            warning("number of zero eigenvalues greater than number of variables")
        idxpos <- idxpos[1] - 1
	if (idxpos < 1) {
            warning("first eigenvalue index zero")
	    gamma <- c()
	} else gamma <- evalue[1:idxpos]
        gamma <- c(gamma, evalue[(idxpos+1+p):N])
        res <- exactMoran(I, gamma, alternative=alternative, type="Global",
            useTP=useTP, truncErr=truncErr, zeroTreat=zeroTreat)
    } else {
        if (dim(Omega)[1] != N) stop("Omega of different size than data")
        res <- exactMoranAlt(I, M, U, Omega, N, alternative=alternative,
            type="Alternative", useTP=useTP, truncErr=truncErr,
            zeroTreat=zeroTreat)
    }

    data.name <- paste("\nmodel:", paste(strwrap(gsub("[[:space:]]+", " ", 
	    paste(deparse(model$call), sep="", collapse=""))), collapse="\n"),
    	    "\nweights: ", deparse(substitute(listw)), "\n", sep="")
    res$estimate <- c(I_save)
    res$data.name <- data.name
    res$df <- (N-p)
    if (!is.null(save.M)) res$M <- M
    if (!is.null(save.U)) res$U <- U
    return(res)
}

# function contributed by Michael Tiefelsdorf 2008
# implements his 2000 book eq. 6.7, p.69

truncPoint <- function(SpecI, truncErr=1e-6, zeroTreat=0.1){
  m <- length(SpecI)
  absSpecI <- abs(SpecI)
  absSpecI[absSpecI < zeroTreat] <- zeroTreat
  TU <- (truncErr*pi*m/2)^(-2/m) * prod(absSpecI^(-1/m))  
# Break product up to reduce rounding errors and over- or under-flows
  return(TU)  
}


exactMoran <- function(I, gamma, alternative="greater", type="Global", np2=NULL, useTP=FALSE, truncErr=1e-6, zeroTreat=0.1) {
    if (!(alternative %in% c("greater", "less", "two.sided")))
	stop("alternative must be one of: \"greater\", \"less\", or \"two.sided\"")
    if (type == "Global") {
        SpecI <- gamma - c(I)
        integrand <- function(x) {
        sin(0.5 * colSums(atan(SpecI %*% t(x)))) /
        (x * apply((1 + SpecI^2 %*% t(x^2))^(1/4), 2, prod))}
    } else if (type == "Alternative") {
        if (useTP) SpecI <- gamma
        integrand <- function(x) {sin(0.5 * colSums(atan((gamma) %*% t(x)))) /
        (x * apply((1+(gamma)^2 %*% t(x^2))^(1/4), 2, prod))}
    } else if (type == "Local") {
        if (is.null(np2)) stop("Local requires np2")
        min <- gamma[1]
        max <- gamma[2]
        if (useTP) SpecI <- c(min, rep(0,np2), max) - I
        integrand <- function(x) {
            sin(0.5*(atan((min-I)*x)+(np2)*atan((-I)*x) + 
            atan((max-I)*x)))/(x*((1+(min-I)^2*x^2) *
            (1+I^2*x^2)^(np2) * (1+(max-I)^2*x^2))^(1/4))} 
    }
    if (useTP) upper <- truncPoint(SpecI, truncErr=truncErr,
        zeroTreat=zeroTreat)
    else upper <- Inf
    II <- integrate(integrand, lower=0, upper=upper)$value
# FIXME II > pi/2
    if (II > pi/2 && type == "Local") {
        tau <- gamma
	df <- np2 + 2
        if (length(tau) == 2L) tau <- c(tau[1], rep(0, df-2), tau[2])
        E.I <- sum(tau)/df
        tau <- tau - E.I
        V.I <- (2*sum(tau^2)) / (df*(df+2))
        sd.ex <- (I - E.I) / sqrt(V.I)
        warning("Normal approximation SD substituted", call.=FALSE)
        oType <- "N"
    } else {
        sd.ex <- qnorm(0.5-II/pi)
        oType <- "E"
    }
    if (alternative == "two.sided") p.v <- 2 * pnorm(abs(sd.ex), 
	lower.tail=FALSE)
    else if (alternative == "greater")
        p.v <- pnorm(sd.ex, lower.tail=FALSE)
    else p.v <- pnorm(sd.ex)
    if (!is.finite(p.v) || p.v < 0 || p.v > 1) 
	warning("Out-of-range p-value: reconsider test arguments")
    statistic <- sd.ex
    attr(statistic, "names") <- "Exact standard deviate"
    p.value <- p.v
    estimate <- c(I)
    attr(estimate, "names") <- "Observed Moran I"
    method <- paste(type, "Moran I statistic with exact p-value")

    res <- list(statistic = statistic, p.value = p.value,
        estimate = estimate, method = method,
	alternative = alternative, gamma=gamma, oType=oType)
    class(res) <- "moranex"
    return(res)
}

exactMoranAlt <- function(I, M, U, Omega, n, alternative="greater",
    type="Alternative", useTP=FALSE, truncErr=1e-6, zeroTreat=0.1) {
    Omega <- chol(Omega)
    A <- Omega %*% M %*% (U - diag(I, n)) %*% M %*% t(Omega)
    gamma <- sort(eigen(A)$values)
    obj <- exactMoran(I, gamma, alternative=alternative,
        type=type, useTP=useTP, truncErr=truncErr, zeroTreat=zeroTreat)
    obj
}


print.moranex <- function(x, ...) {
    class(x) <- c("htest", "moranex")
    print(x, ...)
    invisible(x)
}

H1_moments <- function(M, U, Omega, n) {
    B <- Omega %*% M %*%  t(Omega)
    eigen <- eigen(B)
    lambda <- eigen$values
    P <- eigen$vectors
    A <- Omega %*% M %*% U %*% M %*% t(Omega)
    H <- t(P) %*% A %*% P
    hh <- diag(H)
    integrand <- function(x) apply((1+2*lambda %*% t(x))^(-0.5),2,prod) *
        colSums(hh/(1+2*lambda %*% t(x)))
    mu <- integrate(integrand,lower=0, upper=Inf)$value
    integrand2 <- function(x) {
	res=0
	for (i in 1:n){
	    for(j in 1:n){
		res=res+(H[i,i]*H[j,j]+2*H[i,j]^2) /
                    ((1+2*lambda[i]*x)*(1+2*lambda[j]*x))*x
	    }
	}
	apply((1+2*lambda %*% t(x))^(-0.5),2,prod)*res
    }
    mu2 <- integrate(integrand2,lower=0, upper=Inf)$value
    res<-list(Ew=mu,Var=mu2-mu^2)
    res
}

moranExpect_H1 <- function(listw, rho, select=FALSE){
	if (!(select[1]))
		select=1:length(listw$neighbours)
	n <- length(select)
	V <- listw2mat(listw)[select,select]
	M <- diag(n)-matrix(rep(1/n,n*n),nrow=n)
	WOm <- invIrW(listw, rho=rho)[select,select]
	B <- WOm %*% M %*%  t(WOm)
	eigen <- eigen(B)
	lambda <- eigen$values
	P <- eigen$vectors
	A <- WOm %*% M %*% (0.5*(V + t(V))) %*% M %*% t(WOm)
	H <- t(P) %*% A %*% P
	hh <- diag(H)
	integrand <- function(x) apply((1+2*lambda %*% t(x))^(-0.5), 2, prod) *
		colSums(hh/(1+2*lambda %*% t(x)))
	mu <- integrate(integrand, lower=0, upper=Inf)$value
	mu
}

moranExpect_rho_H1 <- function(listw, I0, select=FALSE){
	if (!(select[1]))
		select=1:length(listw$neighbours)
	V <- listw2mat(listw)[select,select]
	eigenvalues<-as.double(eigen(V)$values)
	min=min(eigenvalues)
	max=max(eigenvalues)
	f <- function(x) abs(I0-moranExpect_H1(listw, x, select=select))
	rop <- optimise(f, c(1/min,1/max))$minimum
	rop
}

