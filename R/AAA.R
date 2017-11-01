# Copyright 2001-15 by Roger Bivand 
#

.spdepOptions <- new.env(FALSE, globalenv())
assign("spChkID", FALSE, envir = .spdepOptions)
assign("zeroPolicy", FALSE, envir = .spdepOptions)
assign("verbose", FALSE, envir = .spdepOptions)
assign("mc", ifelse(.Platform$OS.type == "windows", FALSE, TRUE),
 envir = .spdepOptions)
assign("cores", NULL, envir = .spdepOptions)
assign("cluster", NULL, envir = .spdepOptions)
assign("rlecuyerSeed", rep(12345, 6), envir = .spdepOptions)
assign("listw_is_CsparseMatrix", FALSE, envir = .spdepOptions)

#.conflicts.OK <- TRUE

#.onLoad <- function(lib, pkg) {
#	require(methods)
#}

#.onLoad <- function(pkg, lib) {
#cat("spdep: a package for analysing spatial dependence\n")
#require(maptools)
#.First.lib <- function(lib, pkg) {
#	library.dynam("spdep", pkg, lib)
#}
#.noGenerics <- TRUE

