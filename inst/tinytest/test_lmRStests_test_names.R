library(spdep)
data(oldcol)
oldcrime.lm <- lm(CRIME ~ HOVAL + INC, data = COL.OLD)
lw <- nb2listw(COL.nb)
old.tests <- c("LMerr", "LMlag", "RLMerr", "RLMlag", "SARMA")
set.seed(1)
for (i in 1:4) {
  for (j in 1:5) {
    expect_silent(res <- lm.RStests(oldcrime.lm, listw=lw, test=sample(old.tests, size=i)))
    expect_message(res <- lm.LMtests(oldcrime.lm, listw=lw, test=sample(old.tests, size=i)))
  }
}
