library(spdep)
# https://github.com/r-spatial/spdep/issues/65
p <- c(st_as_sfc("POLYGON ((2 2, 2 8, 8 8, 8 2, 2 2))"), st_as_sfc("POLYGON ((12 2, 12 8, 18 8, 18 2, 12 2))"))
p_sp <- as(p, "Spatial")
expect_true(isTRUE(all.equal(sum(card(poly2nb(p_sp))), 0L)))
expect_true(isTRUE(all.equal(sum(card(poly2nb(p))), 0L)))
expect_true(isTRUE(all.equal(sum(card(poly2nb(p_sp, snap=5))), 2L)))
expect_true(isTRUE(all.equal(sum(card(poly2nb(p, snap=5))), 2L)))
