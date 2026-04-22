library(spdep)
data(boston)
boston.tr <- sf::st_read(system.file("shapes/boston_tracts.gpkg", package="spData")[1])
boston.nb <- poly2nb(boston.tr)
boston.listw <- nb2listw(boston.nb)
a <- moran.plot.drop(boston.c$CMEDV, boston.listw,
 locmoran=localmoran(boston.c$CMEDV, boston.listw), alpha=0.01,
 significant = TRUE, labels = NULL, return_df=TRUE)
b <- moran.plot.drop(boston.c$CMEDV, boston.listw,
 locmoran=NULL, alpha=0.01,
 significant = TRUE, labels = NULL, return_df=TRUE)
expect_true(isTRUE(all.equal(a, b)))
a1 <- moran.plot.drop(boston.c$CMEDV, boston.listw,
 locmoran=localmoran(boston.c$CMEDV, boston.listw), alpha=0.01,
 significant = FALSE, labels = NULL, return_df=TRUE)
b1 <- moran.plot.drop(boston.c$CMEDV, boston.listw,
 locmoran=NULL, alpha=0.01,
 significant = FALSE, labels = NULL, return_df=TRUE)
expect_true(isTRUE(all.equal(a1, b1)))
c <- moran.plot.seismogram(boston.c$CMEDV, boston.listw,
 locmoran=localmoran(boston.c$CMEDV, boston.listw), alpha=0.01,
 zero.policy = TRUE, return_df=TRUE)
d <- moran.plot.seismogram(boston.c$CMEDV, boston.listw,
 locmoran=NULL, alpha=0.01, zero.policy = TRUE, return_df=TRUE)
expect_true(isTRUE(all.equal(c, d)))
