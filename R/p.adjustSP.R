# Copyright 2004 by Roger Bivand and Danlin Yu
#

p.adjustSP <- function(p, nb, method="none") {
# class to inherits Jari Oksanen 080603
  	if (!inherits(nb, "nb")) stop("not a neighbours list")
        n <- card(nb) + 1
        pn <- cbind(p, n)
        res <- apply(pn, 1, function(x) p.adjust(x[1], method=method, n=x[2]))
        res
}

