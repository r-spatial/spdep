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

localC.formula <- function(formula, data, listw, ..., zero.policy=NULL) {
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

localC_perm.default <- function(x, listw, nsim = 499, alternative = "two.sided",
             ..., zero.policy=NULL) {

  alternative <- match.arg(alternative, c("two.sided", "less", "greater"))
  # checks are inherited from localC no need to implement
  obs <- localC(x, listw, zero.policy=zero.policy)

  # if sf object remove geometry & cast as df
  if (inherits(x, "sf")) {
    x[[attr(x, "sf_column")]] <- NULL
    x <- as.data.frame(x)
  }
  
  if (inherits(x, "list")) {
    xorig <- as.matrix(Reduce(cbind, x))
    x <- scale(xorig)
    reps <- localC_perm_calc(x, listw, obs, nsim,
      alternative=alternative, zero.policy=zero.policy)
  }

  if (inherits(x, c("matrix", "data.frame"))) {
    xorig <- as.matrix(x)
    x <- scale(xorig)
    reps <- localC_perm_calc(x, listw, obs, nsim,
      alternative=alternative, zero.policy=zero.policy)
  }

  if (is.vector(x) & is.numeric(x)) {
    xorig <- as.matrix(x)
    x <- scale(xorig)
    reps <- localC_perm_calc(x, listw, obs, nsim, alternative=alternative,
      zero.policy=zero.policy)
  }
  if (ncol(xorig) > 1L) {
    cluster <- rep(3L, length(obs))
    cluster[obs <= reps[, 1]] <- 1L
    cluster[obs > reps[, 1]] <- 2L
    cluster <- factor(cluster, levels=1:3, labels=c("Positive", "Negative",
      "Undefined"))
  } else {
    a <- scale(c(xorig), scale=FALSE)
    b <- lag(listw, a)
    q <- rep(4L, length(a))
    q[a > 0 & b > 0] <- 1L
    q[a <= 0 & b > 0] <- 3L
    q[a <= 0 & b <= 0] <- 2L
    cluster <- factor(q, levels=1:4, labels=c("High-High", "Low-Low",
      "Other Positive", "Negative"))
  }


  attr(obs, "call") <- match.call()
  attr(obs, "pseudo-p") <- reps
  attr(obs, "cluster") <- cluster
  class(obs) <- c("localC", "numeric")

  obs

}

localC_perm.formula <- function(formula, data, listw,
                                nsim = 499, alternative = "two.sided", ...,
                                zero.policy=NULL) {

  alternative <- match.arg(alternative, c("less", "two.sided", "greater"))
  # if any data issues the localC formula method will catch it
  obs <- localC(formula, listw, data, zero.policy=zero.policy)

  # check if sf object in data
  if (inherits(data, "sf")) {
    data[[attr(data, "sf_column")]] <- NULL
    data <- as.data.frame(data)
  }

  xorig <- model.frame(formula, data = data)
  x <- scale()

  reps <- localC_perm_calc(x, listw, obs, nsim, alternative=alternative,
    zero.policy=zero.policy)
  if (ncol(xorig) > 1L) {
    cluster <- rep(3L, length(obs))
    cluster[obs <= reps[, 1]] <- 1L
    cluster[obs > reps[, 1]] <- 2L
    cluster <- factor(cluster, levels=1:3, labels=c("Positive", "Negative",
      "Undefined"))
  } else {
    a <- scale(c(xorig), scale=FALSE)
    b <- lag(listw, a)
    q <- rep(4L, length(a))
    q[a > 0 & b > 0] <- 1L
    q[a <= 0 & b > 0] <- 3L
    q[a <= 0 & b <= 0] <- 2L
    cluster <- factor(q, levels=1:4, labels=c("High-High", "Low-Low",
      "Other Positive", "Negative"))
  }

  attr(obs, "call") <- match.call()
  attr(obs, "pseudo-p") <- reps
  attr(obs, "cluster") <- cluster
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

localC_perm_calc <- function(x, listw, obs, nsim, alternative="two.sided",
  zero.policy=NULL) {
    nc <- ncol(x)
    stopifnot(nc > 0L)
    gr <- punif((1:(nsim+1))/(nsim+1), 0, 1)
    ls <- rev(gr)
    ts <- (ifelse(gr > ls, ls, gr))*2
    if (alternative == "two.sided") {
        probs <- ts
        Prname <- "Pr(z != E(Ci))"
    } else if (alternative == "greater") {
        Prname <- "Pr(z > E(Ci))"
        probs <- gr
    } else {
        Prname <- "Pr(z < E(Ci))"
        probs <- ls
    }
    n <- length(listw$neighbours)
    if (n != nrow(x))stop("Different numbers of observations")
    crd <- card(listw$neighbours)
    permC_int <- function(i, zi, z_i, crdi, wtsi, nsim, Ci, alternative,
      probs, nc) {
      res_i <- rep(as.numeric(NA), 8)
      if (crdi > 0) {
        if (nc == 1L) {
          sz_i <- matrix(sample(c(z_i), size=crdi*nsim, replace=TRUE),
            ncol=crdi, nrow=nsim)
          diffs <- (c(zi) - sz_i)^2
          res_p <- c(diffs %*% wtsi)
        } else {
          res_ps <- matrix(NA, ncol=nc, nrow=nsim)
          for (j in 1:nc) {
            sz_i <- matrix(sample(z_i[,j], size=crdi*nsim, replace=TRUE),
              ncol=crdi, nrow=nsim)
            diffs <- (zi[,j] - sz_i)^2
            res_ps[,j] <- c(diffs %*% wtsi)
          }
          res_p <- apply(res_ps, 1, mean)
        }
# res_p length nsim for obs i conditional draws
        res_i[1] <- mean(res_p)
        res_i[2] <- var(res_p)
        res_i[3] <- (Ci - res_i[1])/sqrt(res_i[2])
        if (alternative == "two.sided")
          res_i[4] <- 2 * pnorm(abs(res_i[3]), lower.tail=FALSE)
        else if (alternative == "greater")
          res_i[4] <- pnorm(res_i[3], lower.tail=FALSE)
        else res_i[4] <- pnorm(res_i[3])
        res_i[5] <- probs[rank(c(res_p, Ci))[(nsim + 1L)]]
        rnk0 <- as.integer(sum(res_p >= Ci))
        drnk0 <- nsim - rnk0
        rnk <- ifelse(drnk0 < rnk0, drnk0, rnk0)
        res_i[6] <- (rnk + 1.0) / (nsim + 1.0)
        res_i[7] <- e1071::skewness(res_p)
        res_i[8] <- e1071::kurtosis(res_p)
      }
      res_i
    }
    z <- scale(x)
    oo <- lapply(1:n, function(i) permC_int(i, z[i,,drop=FALSE], z[-i,], crd[i],
      listw$weights[[i]], nsim, obs[i], alternative, probs, nc))
    res <- do.call("rbind", oo)
    colnames(res) <- c("E.Ci", "Var.Ci", "Z.Ci", Prname,
      paste0(Prname, " Sim"), "Pr(folded) Sim", "Skewness", "Kurtosis")
    res
}


