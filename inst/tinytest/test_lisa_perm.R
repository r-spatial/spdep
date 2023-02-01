library(spdep)
data(boston, package="spData")
lw <- nb2listw(boston.soi)
lws <- nb2listw(include.self(boston.soi))
x <- boston.c$NOX
y <- boston.c$LSTAT
xx <- cbind(boston.c$NOX, boston.c$LSTAT, boston.c$RM)
nsim <- 499L
iseed=1L
expect_silent(no <- system.time(localmoran_perm(x, lw, nsim=nsim, iseed=iseed))["elapsed"])
expect_silent(no <- system.time(localmoran_perm(x, lw, nsim=nsim, iseed=iseed, no_repeat_in_row=TRUE))["elapsed"])
expect_silent(no <- system.time(G <- localG_perm(x, lw, nsim=nsim, iseed=iseed))["elapsed"])
expect_equal(c(G), c(localG(x, lw)))
expect_silent(no <- system.time(Gs <- localG_perm(x, lws, nsim=nsim, iseed=iseed))["elapsed"])
expect_false(isTRUE(all.equal(attr(G, "internals")[,6], attr(Gs, "internals")[,6])))
expect_equal(c(Gs), c(localG(x, lws)))
expect_silent(no <- system.time(Gs_no_rep <- localG_perm(x, lws, nsim=nsim, iseed=iseed, no_repeat_in_row=TRUE))["elapsed"])
expect_silent(no <- system.time(localC_perm(x, lw, nsim=nsim, iseed=iseed))["elapsed"])
expect_silent(no <- system.time(localC_perm(xx, lw, nsim=nsim, iseed=iseed))["elapsed"])
expect_silent(no <- system.time(localmoran_bv(x, y, lw, nsim=nsim, iseed=iseed))["elapsed"])
if (require(parallel, quietly=TRUE)) {
 coresOpt <- get.coresOption()
 nc <- detectCores(logical=FALSE)-1L
 nc
 mcOpt <- get.mcOption()
# set nc to 1L here
 if (nc > 1L) nc <- 1L
 multicore <- snow <- NULL
 if (!is.na(nc)) {
  invisible(set.coresOption(nc))
  if (mcOpt) {
   expect_silent(multicore <- system.time(localmoran_perm(x, lw, nsim=nsim, iseed=iseed))["elapsed"])
   expect_silent(multicore <- system.time(localG_perm(x, lw, nsim=nsim, iseed=iseed))["elapsed"])
   expect_silent(multicore <- system.time(localC_perm(x, lw, nsim=nsim, iseed=iseed))["elapsed"])
   expect_silent(multicore <- system.time(localC_perm(xx, lw, nsim=nsim, iseed=iseed))["elapsed"])
   expect_silent(multicore <- system.time(localmoran_bv(x, y, lw, nsim=nsim, iseed=iseed))["elapsed"])
   invisible(set.mcOption(FALSE))
   cl <- makeCluster(get.coresOption())
   set.ClusterOption(cl)
   expect_silent(snow <- system.time(localmoran_perm(x, lw, nsim=nsim, iseed=iseed))["elapsed"])
   expect_silent(snow <- system.time(localG_perm(x, lw, nsim=nsim, iseed=iseed))["elapsed"])
   expect_silent(snow <- system.time(localC_perm(x, lw, nsim=nsim, iseed=iseed))["elapsed"])
   expect_silent(snow <- system.time(localC_perm(xx, lw, nsim=nsim, iseed=iseed))["elapsed"])
   expect_silent(snow <- system.time(localmoran_bv(x, y, lw, nsim=nsim, iseed=iseed))["elapsed"])
   invisible(stopCluster(cl))
   invisible(set.mcOption(mcOpt))
  } else {
   cl <- makeCluster(get.coresOption())
   set.ClusterOption(cl)
   expect_silent(snow <- system.time(localmoran_perm(x, lw, nsim=nsim, iseed=iseed))["elapsed"])
   expect_silent(snow <- system.time(localG_perm(x, lw, nsim=nsim, iseed=iseed))["elapsed"])
   expect_silent(snow <- system.time(localC_perm(x, lw, nsim=nsim, iseed=iseed))["elapsed"])
   expect_silent(snow <- system.time(localC_perm(xx, lw, nsim=nsim, iseed=iseed))["elapsed"])
   expect_silent(snow <- system.time(localmoran_bv(x, y, lw, nsim=nsim, iseed=iseed))["elapsed"])
   invisible(stopCluster(cl))
  }
  invisible(set.coresOption(coresOpt))
 }
}

