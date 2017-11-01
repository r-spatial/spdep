# Copyright 2001 by Roger Bivand 
#

airdist <- function(ann=FALSE) {
	usr <- diff(par("usr"))[c(1,3)]
	plt <- diff(par("plt"))[c(1,3)]
	if (abs(diff(plt/usr)) > 0.005)
		warning("plot x and y scales may differ: use plot(..., asp=1)")
	coords <- locator(2)
	res <- sqrt(diff(coords$x)^2 + diff(coords$y)^2)
	if (ann) {
		lines(coords)
		text(mean(coords$x), mean(coords$y), format(res, digits=3),
			pos=4, offset=0.2, cex=0.7)
	}
	if (.Platform$OS.type == "windows") bringToTop(-1)
	list(dist=res, coords=coords)
}
