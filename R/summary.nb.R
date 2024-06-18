# Copyright 2001-2019 by Roger Bivand
# Upgrade to sp classes February 2007
#


summary.nb <- function(object, coords=NULL, longlat=NULL, scale=1, ...) {
    nb <- object
    if (!inherits(nb, "nb")) stop("Not a neighbours list")
    c.nb <- card(nb)
    n.nb <- length(nb)
    regids <- attr(nb, "region.id")
    if(is.null(regids)) regids <- as.character(1:n.nb)
    attr(object, "c.nb") <-  c.nb
    attr(object, "n.nb") <-  n.nb
    attr(object, "regids") <- regids
    attr(object, "scale") <- scale
    if(!is.null(object$coords)) {
        dlist <- nbdists(object, coords, longlat=longlat)
        attr(object, "dlist") <- dlist
    }
    class(object) <- c("summary.nb", "nb")
    object
}

print.summary.nb <- function(x, ...) {
    print.nb(x)
    cat("Link number distribution:\n")
    print(table(attr(x, "c.nb"), deparse.level=0))
    if(any(attr(x, "c.nb") > 0)) {
        min.nb <- min(attr(x, "c.nb")[attr(x, "c.nb") > 0])
        cat(length(attr(x, "c.nb")[attr(x, "c.nb") == min.nb]), " least connected region",
	    ifelse(length(attr(x, "c.nb")[attr(x, "c.nb") == min.nb]) < 2L, "", "s"), ":\n",
	    paste(attr(x, "regids")[which(attr(x, "c.nb") == min.nb)], collapse=" "), " with ",
	    min.nb, " link", ifelse(min.nb < 2L, "", "s"), "\n", sep="")
        max.nb <- max(attr(x, "c.nb"))
	cat(length(attr(x, "c.nb")[attr(x, "c.nb") == max.nb]), " most connected region",
	    ifelse(length(attr(x, "c.nb")[attr(x, "c.nb") == max.nb]) < 2L, "", "s"), ":\n",
	    paste(attr(x, "regids")[which(attr(x, "c.nb") == max.nb)], collapse=" "), " with ",
	    max.nb, " link", ifelse(max.nb < 2L, "", "s"), "\n", sep="")
    }
    if(!is.null(attr(x, "dlist"))) {
	cat("Summary of link distances:\n")
	print(summary(unlist(attr(x, "dlist"))))
	stem(unlist(attr(x, "dlist")), scale=attr(x, "scale"))
    }
}


print.nb <- function(x, ...) {
    nb <- x
    if (!inherits(nb, "nb")) stop("Not a neighbours list")
    c.nb <- card(nb)
    n.nb <- length(nb)
    regids <- attr(nb, "region.id")
    if(is.null(regids)) regids <- as.character(1:n.nb)
    s.c.nb <- sum(c.nb)
    cat("Neighbour list object:\n")
    cat("Number of regions:", n.nb, "\n")
    cat("Number of nonzero links:", s.c.nb, "\n")
    cat("Percentage nonzero weights:", (100*sum(c.nb))/(n.nb^2), "\n")
    cat("Average number of links:", mean(c.nb), "\n")
    if(any(c.nb == 0)) cat(length(c.nb[c.nb == 0]), " region", 
        ifelse(length(c.nb[c.nb == 0]) < 2L, "", "s"), " with no links:\n",
	paste(strwrap(paste(regids[which(c.nb == 0)], collapse=" ")),
        collapse="\n"), "\n", sep="")
    nc <- 0
    if (!is.null(attr(x, "ncomp"))) {
        nc <- attr(x, "ncomp")$nc
    } else {
        if (get.SubgraphOption() && get.SubgraphCeiling() > s.c.nb+n.nb) {
            nc <- n.comp.nb(x)$nc
        }
    }
    if (nc > 1) cat(nc, " disjoint connected subgraphs\n", sep="")
    res <- is.symmetric.nb(nb, verbose=FALSE)
    if (!res) cat("Non-symmetric neighbours list\n")
    invisible(x)
}

summary.listw <- function(object, coords=NULL, longlat=FALSE, 
	zero.policy=attr(object, "zero.policy"), scale=1, adjust.n=TRUE, ...) {
        if (is.null(zero.policy))
            zero.policy <- get("zeroPolicy", envir = .spdepOptions)
        stopifnot(is.logical(zero.policy))
        if (any(card(object$neighbours) == 0) && !zero.policy)
            stop("regions with no neighbours found, use zero.policy=TRUE")
	nb_sum <- summary(object$neighbours, coords=coords, longlat=longlat, 
		scale=scale, ...)
        attr(object, "nb_sum") <- nb_sum
        style <- object$style
	if (is.na(style)) style = "NA"
	weights_const <- data.frame(rbind(unlist(spweights.constants(object,
		zero.policy=zero.policy, adjust.n=adjust.n))[c(1, 5:8)]), row.names=style)
        attr(object, "weights_const") <- weights_const
        class(object) <- c("summary.listw", "listw")
        object
}

print.summary.listw <- function(x, ...) {
	cat("Characteristics of weights list object:\n")
        print.summary.nb(attr(x, "nb_sum"))
	cat(paste("\nWeights style:", x$style, "\n"))
	cat("Weights constants summary:\n")
        print(attr(x, "weights_const"))
}

print.listw <- function(x, zero.policy=attr(x, "zero.policy"), ...) {
        if (is.null(zero.policy))
            zero.policy <- get("zeroPolicy", envir = .spdepOptions)
        stopifnot(is.logical(zero.policy))
        if (any(card(x$neighbours) == 0) && !zero.policy)
            stop("regions with no neighbours found, use zero.policy=TRUE")
	cat("Characteristics of weights list object:\n")
	print.nb(x$neighbours, ...)
	style <- x$style
	cat(paste("\nWeights style:", style, "\n"))
	if (is.na(style)) style = "NA"
	cat("Weights constants summary:\n")
        df <- data.frame(rbind(unlist(spweights.constants(x,
		zero.policy=zero.policy))[c(1, 5:8)]), row.names=style)
	print(df)
	invisible(x)

}

