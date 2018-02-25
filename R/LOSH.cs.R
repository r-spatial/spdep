LOSH.cs <- function(x, listw, zero.policy = NULL, na.action = na.fail, 
                 p.adjust.method = "none", spChk = NULL) {
                 
  stopifnot(is.vector(x))
  if (!inherits(listw, "listw")) 
    stop(paste(deparse(substitute(listw)), "is not a listw object"))
  if (is.null(zero.policy)) 
    zero.policy <- get("zeroPolicy", envir = .spdepOptions)
  stopifnot(is.logical(zero.policy))
  if (!is.null(attr(listw$neighbours, "self.included")) && 
      attr(listw$neighbours, "self.included")) 
    stop("Self included among neighbours")
  if (is.null(spChk)) 
    spChk <- get.spChkOption()
  if (spChk && !chkIDs(x, listw)) 
    stop("Check of data and weights ID integrity failed")
  if (!is.numeric(x)) 
    stop(paste(deparse(substitute(x)), "is not a numeric vector"))
  # NAOK <- deparse(substitute(na.action)) == "na.pass"
  x <- na.action(x)
  na.act <- attr(x, "na.action")
  rn <- attr(listw, "region.id")
  if (!is.null(na.act)) {
    subset <- !(1:length(listw$neighbours) %in% na.act)
    listw <- subset(listw, subset, zero.policy = zero.policy)
    excl <- class(na.act) == "exclude"
  }
  ### calculate LOSH for all locations
  res <- LOSH(x, listw, 2, TRUE, zero.policy, na.action, spChk)   #### a = 2 is used, because the chi-square-based inference is only valid for a measure of variance.
  ### add a column for the chi-square-based p-values
  res <- cbind(res, "Pr()" = 1)
  ### p-values using the non-central chi-square distribution
  res[,"Pr()"] <- pchisq(res[,"Z.Hi"], df = 2/res[,"Var.Hi"], lower.tail = FALSE)
  ### correction for multiple hypothesis testing (if p.adjust.method != "none")
  res[,"Pr()"] <- p.adjustSP(res[,"Pr()"], listw$neighbours, method = p.adjust.method)
  if (!is.null(na.act) && excl) {
    res <- naresid(na.act, res)
  }
  if (!is.null(rn)) 
    rownames(res) <- rn
  attr(res, "call") <- match.call()
  if (!is.null(na.act)) 
    attr(res, "na.action") <- na.act
  class(res) <- c("LOSH", class(res))
  res
}
