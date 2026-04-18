library(spdep)
data(oldcol)
lw <- nb2listw(COL.nb, style="W")
is.na(COL.OLD$CRIME) <- 1
expect_error(o <- moran.test(COL.OLD$CRIME, lw))
expect_warning(o <- moran.test(COL.OLD$CRIME, lw, na.action=na.pass))
expect_silent(o <- moran.test(COL.OLD$CRIME, lw, na.action=na.omit))
