# Copyright 2012 by Roger Bivand 
#

isCyclical <- function(nb) {
    stopifnot(inherits(nb, "nb"))
    cnb <- card(nb)
    if (any(cnb == 0)) stop("Neighbours must be connected")
    if (n.comp.nb(nb)$nc != 1) stop("Complete connection required")
    res <- 1L
      for (i in seq(along=nb)) {
        inbs <- nb[[i]]
        if (length(inbs) > 1) {
          for (j in 1:(length(inbs)-1)) {
            for (k in 2:length(inbs)) {
                hit <- (inbs[j] %in% nb[[inbs[k]]])
                if (hit) {
                  res <- 0L
                  break
                }
            }
            if (res == 0L) break
        }
        if (res == 0L) break
      }
    }
    res
}

find_q1_q2 <- function(lw) {
    stopifnot(lw$style == "W")
    nb <- lw$neighbours
    nc <- n.comp.nb(nb)
    members <- tapply(1:length(nb), nc$comp.id, c)
    q2 <- 0L
    q1 <- nc$nc
    t1 <- table(nc$comp.id)
    t2 <- table(t1)
    if ("1" %in% names(t2)) q1 <- unname(q1 - t2["1"])
    ids <- as.integer(names(t1[t1 > 1]))
    members1 <- members[ids]
    nums <- 1:length(nb)
    for (sub in seq(along=members1)) {
        subs <- members1[[sub]]
        nbsub <- subset(nb, (nums %in% subs))
        if (length(nbsub) > 2) q2sub <- isCyclical(nbsub)
        else q2sub <- 1L
        q2 <- q2 + q2sub
    }
    c(q1, q2)
}

