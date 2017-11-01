# Copyright 2003-2011 by Roger Bivand 
#

anova.sarlm <- function(object, ...) {
    if (length(list(object, ...)) > 1L) {
        getResponseFormula <- function (object) 
        {
            form <- formula(object$call)
            if (!(inherits(form, "formula") && (length(form) == 3L))) {
                stop("\"Form\" must be a two sided formula")
            }
            eval(parse(text = paste("~", deparse(form[[2]]))))
        }

        object <- list(object, ...)
        ancall <- sys.call()
        nmodels <- length(object)
        if (nmodels == 1) return(anova(object))
        termsClass <- unlist(lapply(object, data.class))
        if (!all(match(termsClass, c("lm", "sarlm"), 0))) {
            stop(paste("Objects must inherit from classes \"sarlm\" or \"lm\""))
        }
        resp <- unlist(lapply(object, 
	    function(el) deparse(getResponseFormula(el)[[2]])))
        subs <- as.logical(match(resp, resp[1], FALSE))
        if (!all(subs)) 
            warning(paste("Some fitted objects deleted because", 
                "response differs from the first model"))
        if (sum(subs) == 1) 
            stop("First model has a different response from the rest")
        object <- object[subs]
        aux <- lapply(object, logLik)
        if (length(unique(unlist(lapply(object, 
	    function(el) length(residuals(el)))))) > 1L) {
            stop("All fitted objects must use the same number of observations")
        }
        dfModel <- unlist(lapply(aux, function(el) attr(el, "df")))
        logLik <- unlist(lapply(aux, function(el) c(el)))
        AIC <- unlist(lapply(aux, AIC))
        aod <- data.frame(Model = (1:nmodels), df = dfModel, 
        AIC = AIC, logLik = logLik, check.names = FALSE)
        ddf <- diff(dfModel)
        if (sum(abs(ddf)) > 0) {
	   effects <- rep("", nmodels)
	   for (i in 2:nmodels) {
                if (ddf[i - 1] != 0) {
            	    effects[i] <- paste(i - 1, i, sep = " vs ")
                }
            }
            pval <- rep(NA, nmodels - 1)
            ldf <- as.logical(ddf)
            lratio <- 2 * abs(diff(logLik))
            lratio[!ldf] <- NA
            pval[ldf] <- 1 - pchisq(lratio[ldf], abs(ddf[ldf]))
	    aod <- data.frame(aod, Test = effects, L.Ratio = c(NA, 
                lratio), "p-value" = c(NA, pval), check.names = FALSE)
    	}
    	row.names(aod) <- unlist(lapply(as.list(ancall[-1]), 
            deparse))
    	attr(aod, "nmodels") <- nmodels
    	class(aod) <- c("anova", "data.frame")
    	return(aod)

    } else {
    	if (!inherits(object, "sarlm")) 
            stop("object not a fitted simultaneous autoregressive model")
        LL <- logLik(object)
        AIC <- AIC(LL)
        res <- data.frame("AIC"=AIC, "Log likelihood"=LL, "df"=attr(LL, "df"),
	    row.names=deparse(substitute(object)))
	class(res) <- c("anova", "data.frame")
        return(res)
    }
}


