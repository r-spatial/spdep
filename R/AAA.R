# Copyright 2001-24 by Roger Bivand 
#

.spdepOptions <- new.env(TRUE, globalenv())
assign("spChkID", FALSE, envir = .spdepOptions)
assign("zeroPolicy", FALSE, envir = .spdepOptions)
assign("verbose", FALSE, envir = .spdepOptions)
assign("mc", ifelse(.Platform$OS.type == "windows", FALSE, TRUE),
 envir = .spdepOptions)
assign("cores", NULL, envir = .spdepOptions)
assign("cluster", NULL, envir = .spdepOptions)
assign("rlecuyerSeed", rep(12345, 6), envir = .spdepOptions)
assign("listw_is_CsparseMatrix", FALSE, envir = .spdepOptions)
assign("cluster", NULL, envir = .spdepOptions)
assign("report_nb_subgraphs", TRUE, envir = .spdepOptions)
assign("nb_subgraphs_N+E", 100000L, envir = .spdepOptions)
assign("report_nb_noneighs", TRUE, envir = .spdepOptions)
setOldClass(c("listw"))

.onLoad <- function(lib, pkg) {
  options(Matrix.warnDeprecatedCoerce = 1L)
}

#.conflicts.OK <- TRUE

#.onLoad <- function(lib, pkg) {
#	require(methods)
#}

#.onAttach <- function(lib, pkg) {
#packageStartupMessage("spdep: a package for analysing spatial dependence\nDEPRECATED: from 1.1-1, spatial regression functions moved to the spatialreg package\n", appendLF = FALSE)
#require(maptools)
#.First.lib <- function(lib, pkg) {
#	library.dynam("spdep", pkg, lib)
#}
#.noGenerics <- TRUE

