library(spdep)
data(oldcol)
lw <- nb2listw(COL.nb)
COL.OLD$fEW <- factor(COL.OLD$EW)
COL.OLD$fDISCBD <- ordered(cut(COL.OLD$DISCBD, c(0, 1.5, 3, 4.5, 6)))
f <- formula(CRIME ~ INC + HOVAL + fDISCBD + fEW)
lm_obj <- lm(f, data=COL.OLD)
expect_warning(COL.SD0 <- SD.RStests(lm_obj, lw, test="SDM", Durbin=TRUE))
expect_warning(COL.SD1 <- SD.RStests(lm_obj, lw, test="SDM", Durbin=~ INC + HOVAL + fDISCBD + fEW))
expect_warning(COL.SD2 <- SD.RStests(lm_obj, lw, test="SDM", Durbin=~ INC + HOVAL + fEW))
expect_silent(COL.SD3 <- SD.RStests(lm_obj, lw, test="SDM", Durbin=~ INC + HOVAL))
expect_warning(COL.SDE0 <- SD.RStests(lm_obj, lw, test="SDEM", Durbin=TRUE))
expect_warning(COL.SDE1 <- SD.RStests(lm_obj, lw, test="SDEM", Durbin=~ INC + HOVAL + fDISCBD + fEW))
expect_warning(COL.SDE2 <- SD.RStests(lm_obj, lw, test="SDEM", Durbin=~ INC + HOVAL + fEW))
expect_silent(COL.SDE3 <- SD.RStests(lm_obj, lw, test="SDEM", Durbin=~ INC + HOVAL))
mf_obj0 <- lm(f, data=COL.OLD, method="model.frame")
hfp0 <- have_factor_preds_mf(mf_obj0)
expect_equal(unname(attr(hfp0, "pred_contrasts")),
    c("contr.poly", "contr.treatment"))
contrasts(COL.OLD$fDISCBD) <- "contr.treatment"
mf_obj0a <- lm(f, data=COL.OLD, method="model.frame")
hfp0a <- have_factor_preds_mf(mf_obj0a)
expect_equal(unname(attr(hfp0a, "pred_contrasts")),
    c("contr.treatment", "contr.treatment"))
if (require("codingMatrices", quietly=TRUE)) {
    contrasts(COL.OLD$fDISCBD) <- "code_diff"
    mf_obj1 <- lm(f, data=COL.OLD, method="model.frame")
    hfp1 <- have_factor_preds_mf(mf_obj1)
    expect_equal(unname(attr(hfp1, "pred_contrasts")),
        c("code_diff", "contr.treatment"))
    contrasts(COL.OLD$fEW) <- "code_control"
    mf_obj1 <- lm(f, data=COL.OLD, method="model.frame")
    hfp1 <- have_factor_preds_mf(mf_obj1)
    expect_equal(unname(attr(hfp1, "pred_contrasts")),
        c("code_diff", "code_control"))
}
