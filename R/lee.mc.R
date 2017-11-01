# Copyright 20014 by Roger Bivand, Virgilio Gómez-Rubio
#


lee.mc <- function(x, y, listw, nsim, zero.policy=NULL,
	alternative="greater", na.action=na.fail, spChk=NULL,
        return_boot=FALSE) {
	alternative <- match.arg(alternative, c("greater", "less"))
	if(!inherits(listw, "listw")) stop(paste(deparse(substitute(listw)),
		"is not a listw object"))
	if(!is.numeric(x) | !is.numeric(y)) stop(paste(deparse(substitute(x)),
		"is not a numeric vector"))
        if (is.null(zero.policy))
            zero.policy <- get("zeroPolicy", envir = .spdepOptions)
        stopifnot(is.logical(zero.policy))
	if(missing(nsim)) stop("nsim must be given")
	if (is.null(spChk)) spChk <- get.spChkOption()
	if (spChk && !chkIDs(x, listw) && !chkIDs(y, listw))
		stop("Check of data and weights ID integrity failed")
	cards <- card(listw$neighbours)
	if (!zero.policy && any(cards == 0))
		stop("regions with no neighbours found")
#	if (any(is.na(x))) stop("NA in X")
#	if (any(is.na(y))) stop("NA in Y")
	xname <- deparse(substitute(x))
	yname <- deparse(substitute(y))
	wname <- deparse(substitute(listw))
	if (deparse(substitute(na.action)) == "na.pass")
	    stop("na.pass not permitted")

	#Check NA's in both vectors
	na.act <- attr(na.action(cbind(x,y)), "na.action")
	
	x[na.act]<-NA
	y[na.act]<-NA
	
	x<-na.action(x)
	y<-na.action(y)

	if (!is.null(na.act)) {
	    subset <- !(1:length(listw$neighbours) %in% na.act)
	    listw <- subset(listw, subset, zero.policy=zero.policy)
	}
	n <- length(listw$neighbours)
	if ((n != length(x)) |(n != length(y))) stop("objects of different length")
        gamres <- suppressWarnings(nsim > gamma(n + 1))
        if (gamres) stop("nsim too large for this number of observations")
	if (nsim < 1) stop("nsim too small")
	
#	S0 <- Szero(listw)
	S2<-sum ( (unlist(lapply(listw$weights, sum)))^2 )

	#Data frame with x, y
	xy<-data.frame(x,y)
        if (return_boot) {
            lee_boot <- function(var, i, ...) {
#                var <- var[i]
#                return(moran(x=var, ...)$I)
		return(lee(x=var[i,1], y=var[i,2], ...)$L)
            }
            cores <- get.coresOption()
            if (is.null(cores)) {
            parallel <- "no"
            } else {
                parallel <- ifelse (get.mcOption(), "multicore", "snow")
            }
            ncpus <- ifelse(is.null(cores), 1L, cores)
            cl <- NULL
            if (parallel == "snow") {
                cl <- get.ClusterOption()
                if (is.null(cl)) {
                    parallel <- "no"
                    warning("no cluster in ClusterOption, parallel set to no")
                }
            }
            res <- boot(xy, statistic=lee_boot, R=nsim,
                sim="permutation", listw=listw, n=n, S2=S2, 
                zero.policy=zero.policy, parallel=parallel, ncpus=ncpus, cl=cl)
            return(res)
        }
	res <- numeric(length=nsim+1)
	for (i in 1:nsim) 
	{
		idx<-sample(1:n)
		res[i] <- lee(x[idx], y[idx], listw, n, S2,
	    zero.policy)$L
	}
	res[nsim+1] <- lee(x, y, listw, n, S2, zero.policy)$L
	rankres <- rank(res)
	xrank <- rankres[length(res)]
	diff <- nsim - xrank
	diff <- ifelse(diff > 0, diff, 0)
	if (alternative == "less") 
        	pval <- punif((diff + 1)/(nsim + 1), lower.tail=FALSE)
    	else if (alternative == "greater") 
        	pval <- punif((diff + 1)/(nsim + 1))
	if (!is.finite(pval) || pval < 0 || pval > 1) 
		warning("Out-of-range p-value: reconsider test arguments")
	statistic <- res[nsim+1]
	names(statistic) <- "statistic"
	parameter <- xrank
	names(parameter) <- "observed rank"
	method <- "Monte-Carlo simulation of Lee's L"
	data.name <- paste(
	    xname, ", ", yname, "\nweights:",
	    wname, ifelse(is.null(na.act), "", paste("\nomitted:", 
	    paste(na.act, collapse=", "))), 
"\nnumber of simulations + 1:",
	    nsim+1, "\n")
	lres <- list(statistic=statistic, parameter=parameter,
	    p.value=pval, alternative=alternative, method=method, 
	    data.name=data.name, res=res)
	if (!is.null(na.act) ) attr(lres, "na.action") <- na.act
	class(lres) <- c("htest", "mc.sim")
	lres
}

