hotspot <- function(obj, ...) {
  UseMethod("hotspot")
}

hotspot.default <- function(obj, ...) {
  stop("obj not a localmoran, localG or localC object")
}

hotspot.localmoran <- function(obj, Prname, cutoff=0.005, quadrant.type="mean", p.adjust="fdr", droplevels=TRUE, ...) {
  quadrant.type <- match.arg(quadrant.type, c("mean", "median", "pysal"))
  if (is.null(attr(obj, "quadr"))) stop("object has no quadr attribute")
  res <- attr(obj, "quadr")[[quadrant.type]]
  if (!is.factor(res)) stop("chosen quadrant type is not a factor")
  nms <- colnames(obj)
  i <- match(Prname, nms)
  if (is.na(i)) stop("Prname not in column names: ", paste(nms, collapse=" "))
  if (is.null(attr(obj, "quadr"))) stop("object has no quadr attribute")
  is.na(res) <- p.adjust(obj[,i], p.adjust) >= cutoff
  if (droplevels) res <- droplevels(res)
  res
}

hotspot.summary.localmoransad <- function(obj, Prname, cutoff=0.005, quadrant.type="mean", p.adjust="fdr", droplevels=TRUE, ...) {
  quadrant.type <- match.arg(quadrant.type, c("mean", "median", "pysal"))
  if (is.null(attr(obj, "quadr"))) stop("object has no quadr attribute")
  res <- attr(obj, "quadr")[[quadrant.type]]
  if (!is.factor(res)) stop("chosen quadrant type is not a factor")
  nms <- colnames(obj)
  i <- match(Prname, nms)
  if (is.na(i)) stop("Prname not in column names: ", paste(nms, collapse=" "))
  is.na(res) <- p.adjust(obj[,i], p.adjust) >= cutoff
  if (droplevels) res <- droplevels(res)
  res
}

hotspot.data.frame.localmoranex <- function(obj, Prname, cutoff=0.005, quadrant.type="mean", p.adjust="fdr", droplevels=TRUE, ...) {
  quadrant.type <- match.arg(quadrant.type, c("mean", "median", "pysal"))
  if (is.null(attr(obj, "quadr"))) stop("object has no quadr attribute")
  res <- attr(obj, "quadr")[[quadrant.type]]
  if (!is.factor(res)) stop("chosen quadrant type is not a factor")
  nms <- colnames(obj)
  i <- match(Prname, nms)
  if (is.na(i)) stop("Prname not in column names: ", paste(nms, collapse=" "))
  is.na(res) <- p.adjust(obj[,i], p.adjust) >= cutoff
  if (droplevels) res <- droplevels(res)
  res
}

hotspot.localC <- function(obj, Prname, cutoff=0.005, p.adjust="fdr", droplevels=TRUE, ...) {
  stopifnot(!is.null(attr(obj, "pseudo-p")))
  if (is.null(attr(obj, "cluster"))) stop("object has no cluster attribute")
  res <- attr(obj, "cluster")
  if (!is.factor(res)) stop("cluster is not a factor")
  x <- attr(obj, "pseudo-p")
  nms <- colnames(x)
  i <- match(Prname, nms)
  if (is.na(i)) stop("Prname not in column names: ", paste(nms, collapse=" "))
  is.na(res) <- p.adjust(x[,i], p.adjust) >= cutoff
  if (droplevels) res <- droplevels(res)
  res
}

hotspot.localG <- function(obj, Prname, cutoff=0.005, p.adjust="fdr", droplevels=TRUE, ...) {
  stopifnot(!is.null(attr(obj, "internals")))
  if (is.null(attr(obj, "cluster"))) stop("object has no cluster attribute")
  res <- attr(obj, "cluster")
  if (!is.factor(res)) stop("cluster is not a factor")
  x <- attr(obj, "internals")
  nms <- colnames(x)
  i <- match(Prname, nms)
  if (is.na(i)) stop("Prname not in column names: ", paste(nms, collapse=" "))
  is.na(res) <- p.adjust(x[,i], p.adjust) >= cutoff
  if (droplevels) res <- droplevels(res)
  res
}
