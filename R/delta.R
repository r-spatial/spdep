# Copyright 2024 by Roger S. Bivand. 
# doi: 10.1111/gean.12390
# François Bavaud (2024) Measuring and Testing Multivariate Spatial
# Autocorrelation in a Weighted Setting: A Kernel Approach,
# Geographical Analysis (2024) 56, 573-599
spatialdelta <- function(D, w, f, alternative="greater", CF=FALSE) {
    alternative <- match.arg(alternative, c("greater", "less", "two.sided"))
    n <- length(f)
    stopifnot(nrow(D) == n)
    stopifnot(ncol(D) == n)
    stopifnot(nrow(w) == n)
    stopifnot(ncol(w) == n)
    E <- diag(f) %*% w
    H <- diag(n) - rep(1, times=n) %*% t(f) # above eq. 15
    B <- -0.5 * (H %*% D %*% t(H)) # eq. 15
#    di <- (1/Delta_Y)*diag(w %*% B) # eq. 30
    sf <- sqrt(f)
    dsf <- diag(sf)
    dsf1 <- diag(1/sf)
    Kx <- dsf %*% B %*% dsf # eq. 16
    Kw <- dsf %*% w %*% dsf1 - tcrossprod(sf) # eq. 21
    trKx <- sum(diag(Kx))
    trKxKx <- sum(diag(Kx %*% Kx))
    trKwKx <- sum(diag(Kw %*% Kx))
    trKwKw <- sum(diag(Kw %*% Kw))
    trKw <- sum(diag(Kw))
    mubar2 <- (trKwKw/(n-1)) - (trKw/(n-1))^2 # p. 591
    lbar2 <- (trKxKx/(n-1)) - (trKx/(n-1))^2 # p. 591
    trKwKwKw <- sum(diag(Kw %*% Kw %*% Kw))
    mubar3 <- (trKwKwKw/(n-1)) - 3*((trKwKw/(n-1))*(trKw/(n-1))) +
        2*((trKw/(n-1))^3) # p. 591
    trKxKxKx <- sum(diag(Kx %*% Kx %*% Kx))
    lbar3 <- trKxKxKx/(n-1) - 3*((trKxKx/(n-1))*(trKx/(n-1))) +
        2*((trKx/(n-1))^3) # p. 591
    trKwKwKwKw <- sum(diag(Kw %*% Kw %*% Kw %*% Kw))
    mubar4 <- (trKwKwKwKw/(n-1)) - 4*((trKwKwKw/(n-1))*(trKw/(n-1))) +
        6*((trKwKw/(n-1))*((trKw/(n-1))^2)) - 3*((trKw/(n-1))^4)
    trKxKxKxKx <- sum(diag(Kx %*% Kx %*% Kx %*% Kx)) # p. 591
    lbar4 <- (trKxKxKxKx/(n-1)) - 4*((trKxKxKx/(n-1))*(trKx/(n-1))) +
        6*((trKxKx/(n-1))*((trKx/(n-1))^2)) - 3*((trKx/(n-1))^4) # p. 591
    d <- trKwKx/trKx # eq. 12 and 26
    RV <- d*(trKx/(trKwKw*trKxKx)) # eq. 27
    trKw <- sum(diag(Kw))
    Ed <- trKw/(n-1) # eq. 36 
    Vd1 <- 2/((n-2)*((n-1)^2)*(n+1))
    Vd2 <- (n-1)*trKwKw - trKw^2
    Vd3 <- (((n-1)*trKxKx) - trKx^2)/(trKx^2)
    Vd <- Vd1*Vd2*Vd3 # eq. 37 
    VI <- (2/(n^2 - 1))*(sum(diag(w %*% w)) - 1 -((sum(diag(w)) - 1)^2/(n-1)))
    Vx <- ((1/(n-2))*(((n-1)/(trKx^2/trKxKx)) - 1))
    Vd0 <- VI*Vx # eq. 38
    skewd1 <- sqrt((8*(n-2)*(n+1)))/((n-3)*(n+3))
    alphamu <- (mubar3/(mubar2^(3/2)))
    alphalambda <- (lbar3/(lbar2^(3/2)))
    skewd2 <- alphamu * alphalambda
    skewd <- skewd1*skewd2 # eq. 39
    gammamu <- (mubar4/mubar2^2) - 3
    gammalambda <- (lbar4/lbar2^2) - 3
    kurtd1 <- ((3*(n-2)*(n+1)) / ((n-4)*(n-3)*(n-1)*n*(n+3)*(n+5)))
    kurtd2 <- (4*(n^2 - n +2)*gammamu*gammalambda) +
        (4*n^2 - 8*n + 52)*(gammamu + gammalambda)
    kurtd3 <- (4*(5*n^3 - 57*n^2 + 25*n + 169)/((n-2)*(n-1)))
    kurtd <- kurtd1 * (kurtd2 - kurtd3) # eq. 40
    if (!is.finite(kurtd)) CF <- FALSE
    std_d <- (d-Ed)/sqrt(Vd) # eq. 41
    

    if (CF) {
        cf_out <- cornish_fisher(std_d, skewd, kurtd)
        if (cf_out == 0) CF <- FALSE
        std_d <- ifelse(std_d <= 0, std_d + cf_out, std_d - cf_out)
        if (alternative == "two.sided") 
	    pv_d <- 2 * pnorm(abs(std_d), lower.tail=FALSE)
        else if (alternative == "greater")
            pv_d <- pnorm(std_d, lower.tail=FALSE)
        else pv_d <- pnorm(std_d)
    } else {
        if (alternative == "two.sided") 
	    pv_d <- 2 * pnorm(abs(std_d), lower.tail=FALSE)
        else if (alternative == "greater")
            pv_d <- pnorm(std_d, lower.tail=FALSE)
        else pv_d <- pnorm(std_d)
    }
    names(std_d) <- "Standard deviate"
    if (CF) names(std_d) <- paste0(names(std_d), " (Cornish-Fisher corrected)")
    vec <- c(d, Ed, Vd, skewd, kurtd)
    names(vec) <- c("delta", "Expectation", "Variance", "Skewness",
        "Excess Kurtosis")
    method <- ifelse(CF, "Bavaud delta under the Cornish-Fisher correction",
        "Bavaud delta normal approximation")
    data.name <- paste(deparse(substitute(D)), "\nweights:",
	deparse(substitute(w)))
    res <- list(statistic=std_d, p.value=pv_d, estimate=vec,
        alternative=alternative, method=method, data.name=data.name)
    attr(res, "Vd0") <- Vd0
    attr(res, "Kx") <- Kx
    attr(res, "Kw") <- Kw
    attr(res, "f") <- f
    attr(res, "w") <- w
    attr(res, "VI") <- VI
    attr(res, "Vx") <- Vx
    attr(res, "alphamu") <- alphamu
    attr(res, "alphalambda") <- alphalambda
    attr(res, "gammamu") <- gammamu
    attr(res, "gammalambda") <- gammalambda
    class(res) <- c("htest", "spatialdelta")
    res

}

summary.spatialdelta <- function(object, ...) {
# Tab. 2
    n <- length(attr(object, "f"))
    Kx <- attr(object, "Kx")
# features kernels
    lambdabar <- sum(diag(Kx))/(n - 1)
    nu <- sum(diag(Kx))^2 / sum(diag(Kx %*% Kx)) # Eq. 42
    kappa <- attr(object, "Vx") # Eq. 38
    alphalambda <- attr(object, "alphalambda")
    gammalambda <- attr(object, "gammalambda")
# spatial weights kernels
    mubar <- unname(object$estimate["Expectation"]) # Eq. 36
    varI <- attr(object, "VI") # Eq. 38
    alphamu <- attr(object, "alphamu")
    gammamu <- attr(object, "gammamu")
    res <- list(htest=object, lambdabar=lambdabar, nu=nu, kappa=kappa,
        alphalambda=alphalambda, gammalambda=gammalambda, mubar=mubar,
        varI=varI, alphamu=alphamu, gammamu=gammamu)
    class(res) <- c("summary.spatialdelta", class(object))
    res
}

print.summary.spatialdelta <- function(x, digits=getOption("digits"), ...) {
# Tab. 2
    htest <- x[[1]]
    print(htest, digits=digits)
    dn <- unlist(strsplit(htest$data.name, " \nweights: "))
    cat("\nFeatures kernels: ", dn[1], "\n")
    print(unlist(x[2:6]), digits=digits)
    cat("\nSpatial weights kernels: ", dn[2], "\n")
    print(unlist(x[7:10]), digits=digits)
    invisible(x)
}

linearised_diffusive_weights <- function(A, f, t_choice=2) { 
    n <- length(f)
    stopifnot(nrow(A) == n)
    stopifnot(ncol(A) == n)
    rS <- rowSums(A)
    n <- length(rS)
    LA <- (diag(rS) - A) # eq. 7
    dsf <- diag(1/sqrt(f))
    lA <- dsf %*% LA %*% dsf # eq. 7
    if (t_choice == 2) t <- 1/max(eigen(lA)$values) # below eq. 21
    else if (t_choice == 1) t <- min(f/rS) # below eq. 8
    else stop("t_choice 1 or 2")
    res <- diag(n) - t*(diag(1/f) %*% (diag(rS) - A)) # eq. 8
}

metropolis_hastings_weights <- function(A, f) {
    n <- length(f)
    stopifnot(nrow(A) == n)
    stopifnot(ncol(A) == n)
    rS <- rowSums(A)
    P <- diag(1/rS) %*% A
    fP <- diag(f) %*% P
    G <- pmin(fP, t(fP))
    E <- diag(f - rowSums(G)) + G
    diag(1/f) %*% E
}

iterative_proportional_fitting_weights <- function(A, f, g=0.001, iter=1000, tol=1e-10, tol.margins=1e-10, print=FALSE) {
    n <- length(f)
    stopifnot(nrow(A) == n)
    stopifnot(ncol(A) == n)
    if (!requireNamespace("mipfp", quietly=TRUE)) {
        warning("The mipfp package is required for this method,\n returning Metropolis-Hastings weights")
        return(metropolis_hastings_weights(A, f))
    }
    res <- mipfp::Ipfp(seed=(A+g), target.list=list(1, 2),
        target.data=list(f, f), iter=iter, tol=tol, tol.margins=tol.margins,
        print=print)
    diag(1/f) %*% res$x.hat
}

graph_distance_weights <- function(A, f) {
    n <- length(f)
    stopifnot(nrow(A) == n)
    stopifnot(ncol(A) == n)
    if (!requireNamespace("igraph", quietly=TRUE)) {
        warning("The igraph package is required for this method,\n returning Metropolis-Hastings weights")
        return(metropolis_hastings_weights(A, f))
    }
    D <- igraph::distances(igraph::graph_from_adjacency_matrix(A))
    H <- diag(n) - rep(1, times=n) %*% t(f) # above eq. 15
    B <- -0.5 * (H %*% D %*% t(H)) # eq. 15
    c1 <- -1/min(c(B))
    (1 + (c1*B)) %*% diag(f) # eq. 33
}

cornish_fisher <- function(z, s, k) {
    res0 <- (s/6)*((z^2) - 1)
    res1 <- (k/24)*((z^3) - 3*z)
    res2 <- ((s^2)/36)*(2*(z^3) - 5*z)
    res <- res0 + res1 - res2 # eq. 45
    if ((z > 0 && res > z) || (z <= 0 && res > abs(z))) {
        res <- 0
        warning("correction term set to zero")
    }
    res
}

plot_spatialcoords <- function(x, ...) {
  UseMethod("plot_spatialcoords")
}

plot_spatialcoords.default <- function(x, ...) {
  stop("x not a spatialdelta object")
}

plot_spatialcoords.spatialdelta <- function(x, cols=c(1L, 2L), mult=c(1, 1),
    power=1L, fmult=NULL, names=NULL, bg=1, pos=3, cex=0.6, ...) {
# Fig. 1
    cols <- as.integer(cols)
    stopifnot(length(cols) == 2L)
    stopifnot(length(mult) == 2L)
    f <- attr(x, "f")
    if (is.null(names)) names <- as.character(1:length(f))
    stopifnot(length(names) == length(f))
    dsf1 <- diag(1/sqrt(f))
    Kw <- attr(x, "Kw")
# Fig. 4
    power <- as.integer(power)
    stopifnot(is.integer(power))
    stopifnot(power > 0)
    if (power > 1L) {
        res <- Kw
        for (k in 2:power) res <- res %*% Kw
        Kw <- res
    }
    Kw_eig <- eigen(Kw) # eq. 23
    X_hat <- dsf1 %*% Kw_eig$vectors %*% diag(sqrt(abs(Kw_eig$values))) # eq. 24
    stopifnot(all(cols <= ncol(X_hat)))
    X <- mult[1]*X_hat[, cols[1]]
    Y <- mult[2]*X_hat[, cols[2]]
    plot(X, Y, type="n", xlab=paste0("spatial coordinate ", cols[1]),
        ylab=paste0("spatial coordinate ", cols[2]), ...)
    if (is.null(fmult)) {
        rX <- diff(range(X))
        rf <- diff(range(f))
        fmult <- (0.02*rX)/rf
    }
    mf <- f*fmult
    symbols(X, Y, circles=mf, bg=bg, add=TRUE, inches=FALSE)
    text(X, Y, names, pos=pos, cex=cex)
    abline(h=0, v=0, lty=2)
    invisible(list(X=X, Y=Y, circles=mf, names=names))
}

plot_moran <- function(x, y, ...) {
  UseMethod("plot_moran")
}

plot_moran.default <- function(x, y, ...) {
  stop("x not a spatialdelta object")
}

plot_moran.spatialdelta <- function(x, y, fmult=NULL, names=NULL, bg=1,
    pos=3, cex=0.6, ...) {
# Fig. 3 (+ fig. 2)
    f <- attr(x, "f")
    if (is.null(names)) names <- as.character(1:length(f))
    stopifnot(length(names) == length(f))
    stopifnot(is.numeric(y))
    stopifnot(length(y) == length(f))
    lw <- mat2listw(attr(x, "w"),  style="B", row.names=names)
    ly <- lag.listw(lw, y)
    a <- lm(ly ~ y)
    yn <- deparse(substitute(y))
    plot(y, ly, type="n", xlab=paste0(yn, " (",
        formatC(coef(a)[2], format="f", digits=4), ")"),
        ylab=paste0("lagged ", yn), ...)
    if (is.null(fmult)) {
        rX <- diff(range(y))
        rf <- diff(range(f))
        fmult <- (0.02*rX)/rf
    }
    mf <- f*fmult
    symbols(y, ly, circles=mf, bg=bg, add=TRUE, inches=FALSE)
    text(y, ly, names, pos=pos, cex=cex)
    abline(a=a)
    invisible(list(y=y, ly=ly, circles=mf, names=names))
}

factorial_coordinates <- function(x) {
  UseMethod("factorial_coordinates")
}

factorial_coordinates.default <- function(x) {
  stop("x not a spatialdelta object")
}

factorial_coordinates.spatialdelta <- function(x) {
# Fig. 2
    stopifnot(inherits(x, "spatialdelta"))
    dsf1 <- diag(1/sqrt(attr(x, "f")))
    Kx_eig <- eigen(attr(x, "Kx")) # eq. 18
    dsf1 %*% Kx_eig$vectors %*% diag(sqrt(abs(Kx_eig$values))) # eq. 19
}

plot_factorialcoords <- function(x, ...) {
  UseMethod("plot_factorialcoords")
}

plot_factorialcoords.default <- function(x, ...) {
  stop("x not a spatialdelta object")
}


plot_factorialcoords.spatialdelta <- function(x, cols=c(1L, 2L),
    mult=c(1, 1), fmult=NULL, names=NULL, bg=1, pos=3, cex=0.6, ...) {
# Fig. 5 left, not right
    cols <- as.integer(cols)
    stopifnot(length(cols) == 2L)
    stopifnot(length(mult) == 2L)
    f <- attr(x, "f")
    if (is.null(names)) names <- as.character(1:length(f))
    stopifnot(length(names) == length(f))
    x_tilde <- factorial_coordinates(x)
    stopifnot(all(cols <= ncol(x_tilde)))
    X <- mult[1]*x_tilde[, cols[1]]
    Y <- mult[2]*x_tilde[, cols[2]]
    plot(X, Y, type="n", xlab=paste0("factorial coordinate ", cols[1]),
        ylab=paste0("factorial coordinate ", cols[2]))
    if (is.null(fmult)) {
        rX <- diff(range(Re(X)))
        rf <- diff(range(f))
        fmult <- (0.02*rX)/rf
    }
    mf <- f*fmult
    symbols(X, Y, circles=mf, bg=bg, add=TRUE, inches=FALSE)
    text(X, Y, names, pos=pos, cex=cex)
    abline(h=0, v=0, lty=2)
    invisible(list(X=X, Y=Y, circles=mf, names=names))
}


plot_spatialscree <- function(x, ...) {
  UseMethod("plot_spatialscree")
}

plot_spatialscree.default <- function(x, ...) {
  stop("x not a spatialdelta object")
}


plot_spatialscree.spatialdelta <- function(x, ...) {
# Fig. 8
    e <- sort(Re(eigen(attr(x, "w"), only.values=TRUE)$values),
        decreasing=TRUE)
    plot(e, type="n", ylab="Eigenvalues", xlab="", ...)
    zero <- rep(0, length(e))
    n <- 1:length(e)
    segments(n, e, n, zero, lwd=2)
    invisible(e)
}

plot_factorialscree <- function(x, ...) {
  UseMethod("plot_factorialscree")
}

plot_factorialscree.default <- function(x, ...) {
  stop("x not a spatialdelta object")
}


plot_factorialscree.spatialdelta <- function(x, ...) {
# Fig. 7
    e <- sort(Re(eigen(attr(x, "Kx"), only.values=TRUE)$values),
        decreasing=TRUE)
    plot(e, type="n", ylab="Eigenvalues", xlab="", ...)
    zero <- rep(0, length(e))
    n <- 1:length(e)
    segments(n, e, n, zero, lwd=2)
    invisible(e)
}

