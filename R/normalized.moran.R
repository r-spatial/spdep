normalized.moran <- function (x, listw, zero.policy = attr(listw, "zero.policy"),
                              NAOK = FALSE)
{
  if (is.null(zero.policy))
    zero.policy <- get.ZeroPolicyOption()
  stopifnot(is.logical(zero.policy))
  n1 <- length(listw$neighbours)
  x <- c(x)
  if (n1 != length(x))
    stop("objects of different length")
  xx <- mean(x, na.rm = NAOK)
  z <- x - xx
  lz <- spdep::lag.listw(listw, z, zero.policy = zero.policy, NAOK = NAOK)
  if(NAOK) {
    z_lz <- cbind(z,lz)
    z_lz <- na.omit(z_lz)
    z <- z_lz[,1]
    lz <- z_lz[,2]
  }
  I <- cor(z, lz)
  I
}
