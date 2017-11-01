# Copyright 2002-2010 by Roger Bivand and Pierre Legendre
#

sp.correlogram <- function (neighbours, var, order = 1, method = "corr", 
    style = "W", randomisation = TRUE, zero.policy = NULL, spChk = NULL) {
    if (class(neighbours) != "nb") 
        stop("not a neighbours list")
    stopifnot(is.vector(var))
    if (any(is.na(var))) 
        stop("no NAs permitted in variable")
    if (is.null(zero.policy))
        zero.policy <- get("zeroPolicy", envir = .spdepOptions)
    stopifnot(is.logical(zero.policy))
    if (is.null(spChk)) 
        spChk <- get.spChkOption()
    if (spChk && !chkIDs(var, nb2listw(neighbours, zero.policy = zero.policy))) 
        stop("Check of data and weights ID integrity failed")
    if (order < 1) stop("order less than 1")
    nblags <- nblag(neighbours, maxlag = order)
    cardnos <- vector(mode = "list", length = order)
    for (i in 1:order) cardnos[[i]] <- table(card(nblags[[i]]))
    nobs <- sapply(cardnos, function(x) sum(x[names(x) > "0"]))
    if (any(nobs < 3))
        stop("sp.correlogram: too few included observations in higher lags:\n\treduce order.")
    if (method == "corr") {
        lags.x <- matrix(0, nrow = length(var), ncol = order)
        for (i in 1:order) lags.x[, i] <- lag.listw(nb2listw(nblags[[i]], 
            style = style, zero.policy = zero.policy), var, zero.policy = zero.policy)
        res <- cor(cbind(var, lags.x))[1, -1]
        names(res) <- 1:order
    }
    else if ((method == "I") || (method == "C")) {
        res <- matrix(NA, nrow = order, ncol = 3)
        for (i in 1:order) {
            listw <- nb2listw(nblags[[i]], style = style,
                zero.policy = zero.policy)
            if (method == "I") {
                res[i,] <- moran.test(var, listw,
                    randomisation = randomisation,
                    zero.policy = zero.policy)$estimate
            } else {                                          # Addition PL
                res[i, ] <- geary.test(var, listw,
                    randomisation = randomisation, 
                    zero.policy = zero.policy)$estimate           # Addition PL
            }
        }
        rownames(res) <- 1:order
    }
    else stop("method unknown")
    obj <- list(res=res, method=method, cardnos=cardnos, var=deparse(substitute(var)))
    class(obj) <- "spcor"
    obj
}


print.spcor <- function (x, p.adj.method="none", ...) 
{
    cat("Spatial correlogram for", x$var, "\nmethod: ") 
    if (x$method == "I") {
        cat( "Moran's I\n")       # Modif. PL
    } else if(x$method == "C") {
        cat( "Geary's C\n")   # Modif. PL
    } else {
        cat("Spatial autocorrelation\n")          # Modif. PL
    }
    if ((x$method == "I") || (x$method == "C")) {
        res <- as.matrix(x$res)
	ZI <- (res[,1]-res[,2])/sqrt(res[,3])
	pv <- p.adjust(2*pnorm(abs(ZI), lower.tail=FALSE), method=p.adj.method)
	res <- cbind(res, ZI, pv)
        rownames(res) <- paste(rownames(x$res), " (", sapply(x$cardnos,
          function(x) sum(x[names(x) > "0"])), ")", sep="")
        colnames(res) <- c("estimate", "expectation", "variance", 
	    "standard deviate", "Pr(I) two sided")
        printCoefmat(res, ...)
    } else {
        res <- as.vector(x$res)
        names(res) <- names(x$res)
	print(res)
    }
    invisible(res)
}


plot.spcor <- function (x, main, ylab, ylim, ...) 
{
    if (missing(main)) 
        main <- x$var
    if ((x$method == "I") || (x$method == "C")) {
        lags <- as.integer(rownames(x$res))
        to.plot <- which((x$res[,3] > 0) & (x$res[,3] != Inf))  # Addition PL
        sd2 <- rep(0, nrow(x$res))                              # Addition PL
        sd2[to.plot] <- 2 * sqrt(x$res[to.plot, 3])             # Modif. PL
#        sd2 <- 2*sqrt(x$res[,3])
        if (missing(ylim)) {
            ylim <- range(c(x$res[,1]+sd2, x$res[,1]-sd2))
	}
        if (missing(ylab)) 
            if(x$method == "I") ylab <- "Moran's  I"            # Addition PL
            if(x$method == "C") ylab <- "Geary's  C"            # Addition PL
#            ylab <- "Moran's I"
        plot(lags, x$res[,1], type="p", pch=18, ylim = ylim, main = main, ylab = ylab, xaxt = "n")
#        segments(lags, x$res[,1], lags, x$res[,2], lwd=4, col="grey")
        arrows(lags, x$res[,1]+sd2, lags, x$res[,1]-sd2, length=0.1, angle=90)
        arrows(lags, x$res[,1]-sd2, lags, x$res[,1]+sd2, length=0.1, angle=90)
        axis(1, at = lags)
        abline(h = x$res[1,2])
    }
    else {
        res <- as.vector(x$res)
        lags <- as.integer(names(x$res))
        if (missing(ylim)) 
            ylim <- c(-1, 1)
        if (missing(ylab)) 
            ylab <- "Spatial autocorrelation"
        plot(lags, res, type = "h", ylim = ylim, main = main, ylab = ylab, 
            lwd = 4, xaxt = "n")
        axis(1, at = lags)
        abline(h = 0)
    }
}

