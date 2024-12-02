# spdep

[![Actions Status](https://github.com/r-spatial/spdep/workflows/R-CMD-check/badge.svg)](https://github.com/r-spatial/spdep/actions)
[![CRAN](http://www.r-pkg.org/badges/version/spdep)](https://cran.r-project.org/package=spdep)

### Spatial Dependence: Weighting Schemes and Statistics

A collection of functions to create spatial weights matrix objects from polygon contiguities, from point patterns by distance and tessellations, for summarizing these objects, and for permitting their use in spatial data analysis, including regional aggregation by minimum spanning tree; a collection of tests for spatial autocorrelation, including global Morans I and Gearys C proposed by Cliff and Ord (1973, ISBN: 0850860369) and (1981, ISBN: 0850860814), Hubert/Mantel general cross product statistic, Empirical Bayes estimates and Assunção/Reis (1999) (https://doi.org/10.1002/(SICI)1097-0258(19990830)18:16%3C2147%3A%3AAID-SIM179%3E3.0.CO%3B2-I) Index, Getis/Ord G (Getis and Ord 1992) (https://doi.org/10.1111/j.1538-4632.1992.tb00261.x) and multicoloured join count statistics, APLE (Li et al. ) (https://doi.org/10.1111/j.1538-4632.2007.00708.x), local Morans I (Anselin 1995) (https://doi.org/10.1111/j.1538-4632.1995.tb00338.x) and Getis/Ord G (Ord and Getis 1995) (https://doi.org/10.1111/j.1538-4632.1995.tb00912.x), saddlepoint approximations (Tiefelsdorf 2002) (https://doi.org/10.1111/j.1538-4632.2002.tb01084.x) and exact tests for global and local Morans I (Bivand et al. 2009) (https://doi.org/10.1016/j.csda.2008.07.021) and LOSH local indicators of spatial heteroscedasticity (Ord and Getis) (https://doi.org/10.1007/s00168-011-0492-y), with further extensions in Bivand (2022) <doi:10.1111/gean.12319>. The implementation of most of the measures is described in Bivand and Wong (2018) (https://doi.org/10.1007/s11749-018-0599-x). Extra measures contributed by Josiah Parry. Lagrange multiplier tests for spatial dependence in linear models are provided (Anselin et al. 1996) (https://doi.org/10.1016/0166-0462(95)02111-6), as are Rao's score tests for hypothesised spatial Durbin models based in fitted linear models (Koley and Bera 2024) (https://doi.org/10.1080/17421772.2023.2256810). A local indicators for categorical data (LICD) implementation based on Carrer et al. (2021) (https://doi.org/10.1016/j.jas.2020.105306) and Bivand et al. (2017) (https://doi.org/10.1016/j.spasta.2017.03.003) was added in 1.3-7.

**sfdep** (https://cran.r-project.org/package=sfdep) provides a piped interface to **spdep**.

From **spdep** and **spatialreg** versions >= 1.2-1, the model fitting functions previously present in this package are defunct in **spdep** and may be found in **spatialreg**.

Default branch now `main`

Moved from [R-Forge](https://r-forge.r-project.org/projects/spdep/)
