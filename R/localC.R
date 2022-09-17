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


localC_perm <- function(x, ..., zero.policy=NULL, iseed=NULL) {
  UseMethod("localC_perm")
}

localC_perm.default <- function(x, listw, nsim = 499, alternative = "two.sided",
             ..., zero.policy=NULL, iseed=NULL) {

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
      zero.policy=zero.policy, iseed=iseed)
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
                                zero.policy=NULL, iseed=NULL) {

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
    zero.policy=zero.policy, iseed=iseed)
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
  zero.policy=NULL, iseed=NULL) {
    nc <- ncol(x)
    stopifnot(nc > 0L)
    n <- length(listw$neighbours)

    if (n != nrow(x))stop("Different numbers of observations")
    probs <- probs_lut(nsim=nsim, alternative=alternative)
    Prname <- attr(probs, "Prname")

    crd <- card(listw$neighbours)
    z <- scale(x)
    lww <- listw$weights
    env <- new.env()
    assign("z", z, envir=env)
    assign("crd", crd, envir=env)
    assign("lww", lww, envir=env)
    assign("nsim", nsim, envir=env)
    assign("obs", obs, envir=env)
    assign("nc", nc, envir=env)
    varlist <- ls(envir = env)
    permC_int <- function(i, env) {
#    permC_int <- function(i, zi, z_i, crdi, wtsi, nsim, Ci, nc) {
      res_i <- rep(as.numeric(NA), 8)
      crdi <- get("crd", envir=env)[i]
      if (crdi > 0) { # if i has neighbours
        nsim <- get("nsim", envir=env)
        zi <- get("z", envir=env)[i,,drop=FALSE]
        z_i <- get("z", envir=env)[-i,]
        wtsi <- get("lww", envir=env)[[i]]
        nc <- get("nc", envir=env)
        if (nc == 1L) { # if univariate
          sz_i <- matrix(sample(c(z_i), size=crdi*nsim, replace=TRUE),
            ncol=crdi, nrow=nsim) # permute nsim*#neighbours from z[-i]
          diffs <- (c(zi) - sz_i)^2
          res_p <- c(diffs %*% wtsi)
        } else { # else multivariate
          res_p <- numeric(length=nsim) # for cumulation across columns
          sii <- sample.int(nrow(z_i), size=crdi*nsim, replace=TRUE)
          for (j in 1:nc) { # permute nsim*#neighbours row indices from z[-i]
            # create nsim by crdi matrix of z[-i, j] for j-th column
            sz_i <- matrix(z_i[sii, j], ncol=crdi, nrow=nsim) 
            diffs <- (zi[, j] - sz_i)^2
            res_p <- res_p + c(diffs %*% wtsi) # cumulate across columns
          }
          res_p <- res_p/nc # nsim simulated local Gi
        }
        # res_p length nsim for obs i conditional draws
        res_i[1] <- mean(res_p)
        res_i[2] <- var(res_p)
        Ci <- get("obs", envir=env)[i]
        res_i[5] <- rank(c(res_p, Ci))[(nsim + 1L)]
        res_i[6] <- as.integer(sum(res_p >= Ci))
        res_i[7] <- e1071::skewness(res_p)
        res_i[8] <- e1071::kurtosis(res_p)
      }
      res_i
    }

    res <- run_perm(fun=permC_int, n=n, env=env, iseed=iseed, varlist=varlist)

    res[,3] <- (obs - res[,1])/sqrt(res[,2])
    if (alternative == "two.sided")
      res[,4] <- 2 * pnorm(abs(res[,3]), lower.tail=FALSE)
    else if (alternative == "greater")
      res[,4] <- pnorm(res[,3], lower.tail=FALSE)
    else res[,4] <- pnorm(res[,3])
    res[,5] <- probs[as.integer(res[,5])]
    rnk0 <- as.integer(res[,6])
    drnk0 <- nsim - rnk0
    rnk <- ifelse(drnk0 < rnk0, drnk0, rnk0)
    res[,6] <- (rnk + 1.0) / (nsim + 1.0)
    
    colnames(res) <- c("E.Ci", "Var.Ci", "Z.Ci", Prname,
      paste0(Prname, " Sim"), "Pr(folded) Sim", "Skewness", "Kurtosis")
    res
}


