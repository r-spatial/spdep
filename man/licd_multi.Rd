\name{licd_multi}
\alias{licd_multi}

\title{Local Indicators for Categorical Data}

\description{Local indicators for categorical data combine a measure of local composition in a window given by the per-observation set of neighbouring observations, with a local multi-category joincount test simplified to neighbours with the same or different categories compared to the focal observation}

\usage{
licd_multi(fx, listw, zero.policy = attr(listw, "zero.policy"), adjust.n = TRUE,
 nsim = 0L, iseed = NULL, no_repeat_in_row = FALSE, control = list())
}

\arguments{
  \item{fx}{a factor with two or more categories, of the same length as the neighbours and weights objects in listw}
  \item{listw}{a \code{listw} object created for example by \code{nb2listw}}
  \item{zero.policy}{default \code{attr(listw, "zero.policy")} as set when \code{listw} was created, if attribute not set, use global option value; if TRUE assign zero to the lagged value of zones without neighbours, if FALSE assign NA}
  \item{adjust.n}{default TRUE, if FALSE the number of observations is not adjusted for no-neighbour observations, if TRUE, the number of observations is adjusted}
  \item{nsim}{default 0, number of conditonal permutation simulations}
  \item{iseed}{default NULL, used to set the seed; the output will only be reproducible if the count of CPU cores across which computation is distributed is the same}
  \item{no_repeat_in_row}{default \code{FALSE}, if \code{TRUE}, sample conditionally in each row without replacements to avoid duplicate values, \url{https://github.com/r-spatial/spdep/issues/124}}
  \item{control}{comp_binary=TRUE, binomial_punif_alternative="greater",
        jcm_same_punif_alternative="less", jcm_diff_punif_alternative="greater"}
}
\details{The original code may be found at \doi{10.5281/zenodo.4283766}}

\value{
  \item{local_comp}{data.frame object with LICD local composition columns: ID,
    category_i, count_like_i, prop_i, count_nbs_i, pbinom_like_BW,
    pbinom_unlike_BW, pbinom_unlike_BW_alt, chi_BW_i, chi_K_i, anscombe_BW}
  \item{local_config}{data.frame object with LICD local configuration columns: ID, jcm_chi_obs, jcm_count_BB_obs, jcm_count_BW_obs, jcm_count_WW_obs, pval_jcm_obs_BB, pval_jcm_obs_WW, pval_jcm_obs_BW}
  \item{local_comp_sim}{data.frame object with permutation-based LICD local composition columns: ID, pbinom_like_BW, pbinom_unlike_BW, pbinom_unlike_BW_alt, rank_sim_chi_BW, rank_sim_chi_K, rank_sim_anscombe}
  \item{local_config_sim}{data.frame object with permutation-based LICD local configuration columns: ID, jcm_chi_sim_rank, pval_jcm_obs_BB_sim, pval_jcm_obs_BW_sim, pval_jcm_obs_WW_sim}
}
\references{
Cliff, A. D., Ord, J. K. 1981 Spatial processes, Pion, p. 20; 

Upton, G., Fingleton, B. 1985 Spatial data analysis by example: point pattern and qualitative data, Wiley, pp. 158--170; 

Boots, B., 2003. Developing local measures of spatial association for categorical data. Journal of Geographical Systems 5, 139–160; 

Boots, B., 2006. Local configuration measures for categorical spatial data: binary regular lattices. Journal of Geographical Systems 8 (1), 1–24; 

Pietrzak, M.B., Wilk, J., Kossowski, T., Bivand, R.S., 2014. The application of local indicators for categorical data (LICD) in the spatial analysis of economic development. Comparative Economic Research 17 (4), 203–220 \doi{10.2478/cer-2014-0041}; 

Bivand, R.S., Wilk, J., Kossowski, T., 2017. Spatial association of population pyramids across Europe: The application of symbolic data, cluster analysis and join-count tests. Spatial Statistics 21 (B), 339–361 \doi{10.1016/j.spasta.2017.03.003}; 

Francesco Carrer, Tomasz M. Kossowski, Justyna Wilk, Michał B. Pietrzak, Roger S. Bivand, The application of Local Indicators for Categorical Data (LICD) to explore spatial dependence in archaeological spaces, Journal of Archaeological Science, 126, 2021, \doi{10.1016/j.jas.2020.105306}
}
\author{Roger Bivand \email{Roger.Bivand@nhh.no} based on earlier code by Tomasz M. Kossowski, Justyna Wilk and Michał B. Pietrzak}

\note{In order to increase the numbers of neighbours using \code{\link{nblag}} and \code{\link{nblag_cumul}} is advisable; use of binary weights is advised and are in any case used for the composition measure}

\seealso{\code{\link{joincount.multi}}}

\examples{
columbus <- st_read(system.file("shapes/columbus.gpkg", package="spData")[1], quiet=TRUE)
HICRIME <- cut(columbus$CRIME, breaks=c(0,35,80), labels=c("low","high"))
(nb <- poly2nb(columbus))
lw <- nb2listw(nblag_cumul(nblag(nb, 2)), style="B")
obj <- licd_multi(HICRIME, lw)
str(obj)
h_obj <- hotspot(obj)
str(h_obj)
table(h_obj$both_recode)
columbus$both <- h_obj$both_recode
plot(columbus[, "both"])
}

\keyword{spatial}
