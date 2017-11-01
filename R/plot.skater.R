`plot.skater` <-
function(x, coords, label.areas=NULL,
                        groups.colors, cex.circles=1, cex.labels=1, ...){
  n <- nrow(coords)
  if (is.null(label.areas))
    label.areas <- as.character(1:n)
  gr.lab <- unique(x$groups)
  if (missing(groups.colors))
    groups.colors <- rainbow(length(gr.lab))
  symbols(coords[,1], coords[,2], circles=rep(cex.circles,n),
          inches=FALSE, xlab=" ", ylab=" ", xaxt="n", yaxt="n",
          fg=groups.colors[x$groups], ...)
  id.edgp <- which(sapply(x$edges.groups, function(x)
                          length(x$node))>1L)
  if (length(id.edgp)>0L)
    for (i in 1:length(id.edgp)) {
      id1 <- x$edges.groups[[id.edgp[i]]]$edge[,1]
      id2 <- x$edges.groups[[id.edgp[i]]]$edge[,2]
      segments(coords[id1,1], coords[id1,2],
               coords[id2,1], coords[id2,2],
               col=groups.colors[id.edgp[i]], ...)
    }
  text(coords[,1], coords[,2], label.areas, cex=cex.labels)
  invisible()
}

