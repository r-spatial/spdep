# Copyright 2025 by Roger Bivand 
#

have_factor_preds_mf <- function(mf) {
    if (!inherits(mf, "data.frame") || is.null(attr(mf, "terms")))
        stop("mf not a model.frame")
    dcs <- attr(attr(mf, "terms"), "dataClasses")
    dcfact <- which(!is.na(match(dcs, c("ordered", "factor"))))
    have_factor_preds <- FALSE
    if (length(dcfact) > 0) {
        have_factor_preds <- TRUE
        factnames <- names(dcs)[dcfact]
        xlevels <- lapply(factnames, function(xnms) levels(mf[[xnms]]))
        names(xlevels) <- factnames
        attr(have_factor_preds, "xlevels") <- xlevels 
        attr(have_factor_preds, "factnames") <- factnames
    }
    have_factor_preds
}

warn_factor_preds <- function(x) {
    warning("use of spatially lagged factors (categorical variables)\n", 
        paste(attr(x, "factnames"), collapse=", "),
        "\nis not well-understood")
}
