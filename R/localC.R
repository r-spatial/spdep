localC <- function(x, ..., zero.policy=NULL) {
  UseMethod("localC")
}


localC.default <- function(x, listw, ..., zero.policy=NULL) {
  # check listw object
  if (!inherits(listw, "listw"))
    stop(paste(deparse(substitute(listw)), "is not a listw object"))

  # check missing values
  if (any(is.na(x))) stop(paste("NA in ", deparse(substitute(x))))

  localC_calc(scale(x), listw, zero.policy=zero.policy)
}

localC.formula <- function(formula, listw, data, ..., zero.policy=NULL) {
  # check listw object
  if (!inherits(listw, "listw"))
    stop(paste(deparse(substitute(listw)), "is not a listw object."))

  # check if sf object in data
  if (inherits(data, "sf")) {
    data[[attr(data, "sf_column")]] <- NULL
    data <- as.data.frame(data)
  }

  if (any(sapply(data, class) == "list"))
    stop("`data` cannot contain list columns or elements.")

  form_terms <- terms(formula, data = data)

  if (!attr(form_terms, "response") == 0)
    stop("`formula` must be one-sided with no response variable.")

  df <- model.frame(formula, data = data)

  char_cols <- colnames(df)[sapply(df, class) == "character"]

  if (length(char_cols) > 0)
    stop(paste("Formula contains character vectors:", char_cols))

  rowSums(apply(scale(df), 2, localC_calc, listw, zero.policy=zero.policy)) / ncol(df)

}

localC.list <- function(x, listw, ..., zero.policy=NULL) {

  if (!inherits(listw, "listw"))
    stop(paste(deparse(substitute(listw)), "is not a listw object,"))

  if (!length(unique(lengths(x))) == 1) {
    stop("Elements of x must be of equal length.")
  }

  x <- scale(Reduce(cbind, x))
  rowSums(apply(x, 2, localC_calc, listw, zero.policy=zero.policy)) / ncol(x)

}


localC.matrix <- function(x, listw, ..., zero.policy=NULL) {

  if (!inherits(listw, "listw"))
    stop(paste(deparse(substitute(listw)), "is not a listw object"))

  if (inherits(x, "character")) stop("x must be a numeric matrix.")

  rowSums(apply(scale(x), 2, localC_calc, listw, zero.policy=zero.policy)) / ncol(x)
}

localC.data.frame <- function(x, listw, ..., zero.policy=NULL) {

  if (inherits(x, "sf")) {
    x[[attr(x, "sf_column")]] <- NULL
    x <- as.data.frame(x)
  }

  if (!all(sapply(x, class) %in% c("numeric", "integer")))
    stop("Columns of x must be numeric.")

  rowSums(apply(scale(x), 2, localC_calc, listw, zero.policy=zero.policy)) / ncol(x)

}


localC_perm <- function(x, ..., zero.policy=NULL) {
  UseMethod("localC_perm")
}

localC_perm.default <- function(x, listw, nsim = 499, alternative = "less", ..., zero.policy=NULL) {

  # checks are inherited from localC no need to implement
  obs <- localC(x, listw)

  # if sf object remove geometry & cast as df
  if (inherits(x, "sf")) {
    x[[attr(x, "sf_column")]] <- NULL
    x <- as.data.frame(x)
  }

  if (inherits(x, "list")) {
    x <- scale(as.matrix(Reduce(cbind, x)))
    reps <- replicate(nsim, localC(x[sample.int(nrow(x)),], listw, zero.policy=zero.policy))
  }

  if (inherits(x, c("matrix", "data.frame"))) {
    reps <- replicate(nsim, localC(x[sample.int(nrow(x)),], listw, zero.policy=zero.policy))
  }

  if (is.vector(x) & is.numeric(x)) {
    reps <- replicate(nsim, localC(x[sample.int(length(x))], listw, zero.policy=zero.policy))
  }


  pseudo_p <- localC_p(reps, obs, alternative, nsim)

  attr(obs, "call") <- match.call()
  attr(obs, "pseudo-p") <- pseudo_p
  class(obs) <- c("localC", "numeric")

  obs

}

localC_perm.formula <- function(formula, listw, data,
                                nsim = 499, alternative = "less", ..., zero.policy=NULL) {

  # if any data issues the localC formula method will catch it
  obs <- localC(formula, listw, data)

  # check if sf object in data
  if (inherits(data, "sf")) {
    data[[attr(data, "sf_column")]] <- NULL
    data <- as.data.frame(data)
  }

  x <- scale(model.frame(formula, data = data))

  reps <- replicate(nsim, localC(x[sample.int(nrow(x)),], listw, zero.policy=zero.policy))

  pseudo_p <- localC_p(reps, obs, alternative, nsim)


  attr(obs, "call") <- match.call()
  attr(obs, "pseudo-p") <- pseudo_p
  class(obs) <- c("localC", "numeric")

  obs

}




# Local Geary Utils -------------------------------------------------------
localC_calc <- function(x, listw, zero.policy=NULL) {
  if (any(card(listw$neighbours) == 0L)) {
    res <- geary.intern(x, listw, n=length(listw$neighbours), zero.policy=zero.policy)
  } else {
    xij <- lapply(listw$neighbours, FUN = function(nbs_i) x[nbs_i])
# xij: list of vectors: for each i, x[j] values of its n_i neighbours
    res <- mapply(function(x, j, wi) sum(wi * (j - x)^2),
                x, xij, listw$weights,
                USE.NAMES = FALSE)
# res: numeric vector: for each i, the sum over j of w_{ij} * (x[j] - x[i])^2
  }
  res
}

localC_p <- function(reps, obs, alternative, nsim) {

  alternative <- match.arg(alternative, c("less", "greater"))

  switch(alternative,
         less = (rowSums(reps <= obs) + 1)/ (nsim + 1),
         greater = (rowSums(reps >= obs) + 1)/ (nsim + 1))

}

