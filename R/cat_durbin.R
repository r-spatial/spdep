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
        pred_contrasts <- character(length(factnames))
        pred_ordered <- logical(length(factnames))
        for (pred in seq(along=factnames)) {
            contr <- C(mf[[factnames[pred]]])
            pred_contrasts[pred] <- attr(contr, "contrasts")
            if (pred_contrasts[pred] == "contr.poly" &&
                colnames(contrasts(mf[[factnames[pred]]]))[1] != ".L")
                pred_contrasts[pred] <- as.character(NA)
            pred_ordered[pred] <- names(attr(contr, "contrasts")) == "ordered"
        }
        names(pred_contrasts) <- names(pred_ordered) <- factnames
        attr(have_factor_preds, "pred_contrasts") <- pred_contrasts
        attr(have_factor_preds, "pred_ordered") <- pred_ordered
    }
    have_factor_preds
}

warn_factor_preds <- function(x) {
    plural <- length(attr(x, "factnames")) > 1L
    warning("use of spatially lagged ", ifelse(plural, "factors", "factor"),
        " (categorical ", ifelse(plural, "variables", "variable"), ")\n", 
        paste(attr(x, "factnames"), collapse=", "),
        "\nis not well-understood")
    pred_ordered <- attr(x, "pred_ordered")
    pred_contrasts <- attr(x, "pred_contrasts")
    if (any(pred_ordered & !is.na(pred_contrasts) &
        pred_contrasts == "contr.poly")) {
        ordered <- which(pred_ordered & !is.na(pred_contrasts) &
            pred_contrasts == "contr.poly")
        plural <- length(ordered) > 1L
        warning("In addition ", ifelse(plural, "variables", "variable"), ":\n",
            paste(names(pred_ordered)[ordered], collapse=", "), 
            "\n", ifelse(plural, "are", "is"), 
            " ordered (ordinal) with polynomial contrasts.")
    }
}
