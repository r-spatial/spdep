library(spdep)
library(sf)
columbus <- st_read(system.file("shapes/columbus.shp", package="spData")[1], quiet=TRUE)
col_geoms <- st_geometry(columbus)
col_geoms[21] <- st_buffer(col_geoms[21], dist=-0.05)
st_geometry(columbus) <- col_geoms
expect_false(get.SubgraphOption())
expect_silent(nb <- poly2nb(columbus))
set.SubgraphOption(TRUE)
expect_true(get.SubgraphOption())
expect_warning(nb <- poly2nb(columbus))

