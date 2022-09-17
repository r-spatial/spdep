library(spdep)
data(boston, package="spData")
lw <- nb2listw(boston.soi)
x <- boston.c$NOX
xx <- cbind(boston.c$NOX, boston.c$LSTAT, boston.c$RM)
nsim <- 499L
iseed=1L
expect_silent(no <- system.time(localmoran_perm(x, lw, nsim=nsim, iseed=iseed))["elapsed"])
expect_silent(no <- system.time(localG_perm(x, lw, nsim=nsim, iseed=iseed))["elapsed"])
expect_silent(no <- system.time(localC_perm(x, lw, nsim=nsim, iseed=iseed))["elapsed"])
expect_silent(no <- system.time(localC_perm(xx, lw, nsim=nsim, iseed=iseed))["elapsed"])
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
   invisible(set.mcOption(FALSE))
   cl <- makeCluster(get.coresOption())
   set.ClusterOption(cl)
   expect_silent(snow <- system.time(localmoran_perm(x, lw, nsim=nsim, iseed=iseed))["elapsed"])
   expect_silent(snow <- system.time(localG_perm(x, lw, nsim=nsim, iseed=iseed))["elapsed"])
   expect_silent(snow <- system.time(localC_perm(x, lw, nsim=nsim, iseed=iseed))["elapsed"])
   expect_silent(snow <- system.time(localC_perm(xx, lw, nsim=nsim, iseed=iseed))["elapsed"])
   invisible(stopCluster(cl))
   invisible(set.mcOption(mcOpt))
  } else {
   cl <- makeCluster(get.coresOption())
   set.ClusterOption(cl)
   expect_silent(snow <- system.time(localmoran_perm(x, lw, nsim=nsim, iseed=iseed))["elapsed"])
   expect_silent(snow <- system.time(localG_perm(x, lw, nsim=nsim, iseed=iseed))["elapsed"])
   expect_silent(snow <- system.time(localC_perm(x, lw, nsim=nsim, iseed=iseed))["elapsed"])
   expect_silent(snow <- system.time(localC_perm(xx, lw, nsim=nsim, iseed=iseed))["elapsed"])
   invisible(stopCluster(cl))
  }
  invisible(set.coresOption(coresOpt))
 }
}

