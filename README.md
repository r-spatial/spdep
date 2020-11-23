# spdep

[![Actions Status](https://github.com/r-spatial/spdep/workflows/R-CMD-check/badge.svg)](https://github.com/r-spatial/spdep/actions)
[![CRAN](http://www.r-pkg.org/badges/version/spdep)](https://cran.r-project.org/package=spdep)

### Spatial Dependence: Weighting Schemes and Statistics

A collection of functions to create spatial weights matrix objects from polygon 'contiguities', from point patterns by distance and tessellations, for summarizing these objects, and for permitting their use in spatial data analysis, including regional aggregation by minimum spanning tree; a collection of tests for spatial 'autocorrelation', including global 'Morans I' and 'Gearys C' proposed by 'Cliff' and 'Ord' (1973, ISBN: 0850860369) and (1981, ISBN: 0850860814), 'Hubert/Mantel' general cross product statistic, Empirical Bayes estimates and 'Assunção/Reis' (1999) https://doi.org/10.1002/(SICI)1097-0258(19990830)18:16%3C2147::AID-SIM179%3E3.0.CO;2-I Index, 'Getis/Ord' G ('Getis' and 'Ord' 1992) https://doi.org/10.1111/j.1538-4632.1992.tb00261.x and multicoloured join count statistics, 'APLE' ('Li 'et al.' ) https://doi.org/10.1111/j.1538-4632.2007.00708.x, local 'Moran's I' (Anselin 1995) https://doi.org/10.1111/j.1538-4632.1995.tb00338.x and 'Getis/Ord' G ('Ord' and 'Getis' 1995) https://doi.org/10.1111/j.1538-4632.1995.tb00912.x, 'saddlepoint' approximations ('Tiefelsdorf' 2002) https://doi.org/10.1111/j.1538-4632.2002.tb01084.x and exact tests for global and local 'Moran's I' ('Bivand et al.' 2009) https://doi.org/10.1016/j.csda.2008.07.021 and 'LOSH' local indicators of spatial heteroscedasticity ('Ord' and 'Getis') https://doi.org/10.1007/s00168-011-0492-y. The implementation of most of the measures is described in 'Bivand' and 'Wong' (2018) https://doi.org/10.1007/s11749-018-0599-x.

 **spdep** >= 1.1-1 corresponds to **spatialreg** >= 1.1-1, in which the model fitting functions are deprecated and pass through to **spatialreg**, but will mask those in **spatialreg**. From versions 1.2-1, the functions will be made defunct in **spdep**.

For now **spatialreg** only has functions from **spdep**, where they are shown as deprecated. **spatialreg** only loads the namespace of **spdep**; if you attach it, the same functions in the other package will be masked. Some feed through adequately, others do not (mostly where `stats::model.matrix()` facilities do not like the extra level of passing arguments).

Moved from [R-Forge](https://r-forge.r-project.org/projects/spdep/)
