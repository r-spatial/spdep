# Copyright 2002-12 by Roger Bivand, 2015 Martin Gubri
#

residuals.sarlm <- function(object, ...) {
  if (is.null(object$na.action))
    object$residuals
  else napredict(object$na.action, object$residuals)
}

deviance.sarlm <- function(object, ...) {
  object$SSE
}

coef.sarlm <- function(object, ...) {
  ret <- NULL
  #	ret <- sqrt(object$s2)
  #	names(ret) <- "sigma"
  if (object$type == "error") ret <- c(ret, object$lambda)
  else if (object$type == "lag" || object$type == "mixed")
    ret <- c(ret, object$rho)
  else if (object$type == "sac" || object$type == "sacmixed")
    ret <- c(ret, object$rho, object$lambda)
  ret <- c(ret, object$coefficients)
  
  ret
}

vcov.sarlm <- function(object, ...) {
  if (object$ase) res <- object$resvar[-1,-1]
  else {
    if (!is.null(object$fdHess)) {
      if (object$insert) res <- object$resvar[-1,-1]
      else res <- object$resvar
    } else {
      stop("vcov not available for this model")
    }
  }
  res
}


fitted.sarlm <- function(object, ...) {
  if (is.null(object$na.action))
    object$fitted.values
  else napredict(object$na.action, object$fitted.values)
}


predict.sarlm <- function(object, newdata=NULL, listw=NULL, pred.type="TS", all.data=FALSE,
                          zero.policy=NULL, legacy=TRUE, legacy.mixed=FALSE, power=NULL, order=250, tol=.Machine$double.eps^(3/5), #pred.se=FALSE, lagImpact=NULL, 
                          spChk=NULL, ...) {
  if (is.null(zero.policy))
    zero.policy <- get("zeroPolicy", envir = .spdepOptions)
  stopifnot(is.logical(zero.policy))
  if (is.null(pred.type)) pred.type <- "TS"
  # check pred.type with model
  if (pred.type %in% c("TS") & object$type %in% c("sac", "sacmixed")) stop("no such predict method for sac model")
  if (pred.type %in% c("TC", "BP", "BPW", "BPN", "TS1", "BP1", "BPW1", "BPN1") & object$type == "error") stop("no such predict method for error model")
  if (pred.type %in% c("KP5") & object$type %in% c("lag", "lagmixed")) stop("no such predict method for lag model")
  
  if (pred.type %in% c("TC", "BP", "BPW", "BPN", "BP1", "BPW1", "BPN1") & object$type %in% c("sac", "sacmixed")) warning("predict method developed for lag model, use carefully")
  if (pred.type %in% c("KP5") & object$type %in% c("sac", "sacmixed")) warning("predict method developed for sem model, use carefully")
  
  
  if (is.null(power)) power <- object$method != "eigen"
  stopifnot(is.logical(all.data))
  stopifnot(is.logical(legacy))
  stopifnot(is.logical(legacy.mixed))
  stopifnot(is.logical(power))
  if (is.null(spChk)) spChk <- get.spChkOption()
  if (!is.null(newdata) && is.null(row.names(newdata))) stop("newdata should have region.id as row.names")
  
  #        if (pred.se && object$type == "error") {
  #            pred.se <- FALSE
  #            warning("standard error estimates not available for error models")
  #        }
  #        if (pred.se && is.null(lagImpact))
  #            stop("lagImpact object from impact method required for standard error estimate")
  
  Xs <- object$X
  B <- object$coefficients
  ys <- object$y
  if (class(ys) == "AsIs") ys <- c(ys)
  tarXs <- object$tarX
  tarys <- object$tary
  trends <- Xs %*% B
  
  # forecast case: newdata with the same names than data
  # use a sub-samble of in-sample predictors
  if (!is.null(newdata) && nrow(newdata) == length(ys) && row.names(newdata) == attr(ys, "names")) {
    if (!pred.type %in% c("trend", "TC", "TS"))
      warning("no such predictor type for prevision")
    # DATA
    frm <- formula(object$call)
    mt <- delete.response(terms(frm, data=newdata)) # returns a terms object for the same model but with no response variable
    mf <- model.frame(mt, newdata)
    # resolved problem of missing response column in newdata reported by
    # Christine N. Meynard, 060201
    if (dim(mf)[1] != nrow(newdata))
      stop("missing values in newdata")
    Xs <- model.matrix(mt, mf)
    
    if (object$type == "mixed" || (object$type == "error" && object$etype == "emixed")) { # mixed model: compute WXo
      if (is.null(listw) || !inherits(listw, "listw"))
        stop ("spatial weights list required")
      if (nrow(Xs) != length(listw$neighbours))
        stop("mismatch between data and spatial weights")
      if (spChk && !chkIDs(Xs, listw))
        stop("Check of data and weights ID integrity failed")
      if (pred.type == "TS" && legacy.mixed == FALSE) { # force legacy.mixed to allow the TS out-of-sample computation in the prevision case 
        warning("only legacy.mixed=TRUE is supported for pred.type='TS' and mixed models. legacy.mixed is forced")
        legacy.mixed <- TRUE
      }
      
      K <- ifelse(colnames(Xs)[1] == "(Intercept)", 2, 1)
      m <- ncol(Xs)
      # check if there are enough regressors
      if (m > 1) {
        WXs <- matrix(nrow=nrow(Xs),ncol=(m-(K-1)))
        for (k in K:m) {
          wx <- lag.listw(listw, Xs[,k], 
                          zero.policy=zero.policy)
          if (any(is.na(wx))) 
            stop("NAs in lagged independent variable")
          WXs[,(k-(K-1))] <- wx
        }
      }
      if (K == 2) {
        # unnormalized weight matrices
        if (!(listw$style == "W")) {
          intercept <- as.double(rep(1, nrow(Xs)))
          wx <- lag.listw(listw, intercept, 
                          zero.policy = zero.policy)
          if (m > 1) {
            WXs <- cbind(wx, WXs)
          } else {
            WXs <- matrix(wx, nrow = nrow(Xs), ncol = 1)
          }
        } 
      }
      if (any(object$aliased)) {
        if (K>1 && (listw$style == "W")) colnames(WXs) <- paste("lag.", colnames(Xs)[-1], sep="")
        else colnames(WXs) <- paste("lag.", colnames(Xs), sep="")
      }
      Xs <- cbind(Xs, WXs)
    }
    #  accommodate aliased coefficients 120314
    if (any(object$aliased))
      Xs <- Xs[,-which(object$aliased)]
    
    tarXs <- tarys <- NULL
    trends <- Xs %*% B
    if (pred.type != "TS") # use the out-of-sample computation of TS (which only depends on Woo)
      newdata <- NULL
  }
  
  
  if (is.null(newdata)) { # in-sample pred
    if (pred.type == "TS") { # defaut predictor
      res <- fitted.values(object)
      if (object$type == "error") { # We ou WX+We
        attr(res, "trend") <- as.vector(trends)
        attr(res, "signal") <- as.vector( -1 * (tarys - ys) -	-1 * (tarXs - Xs) %*% B)
      } else { # lag model or other
        attr(res, "trend") <- as.vector(trends)
        attr(res, "signal") <- as.vector( -1 * (tarys - ys))
      }
    } else { # new predictors
      if (pred.type != "trend") { # need listw
        if (is.null(listw) || !inherits(listw, "listw"))
          stop ("spatial weights list required")
        if (nrow(Xs) != length(listw$neighbours))
          stop("mismatch between data and spatial weights")
        if (spChk && !chkIDs(Xs, listw))
          stop("Check of data and weights ID integrity failed")
      }
      
      if (pred.type %in% c("TC", "BP")) { # need to compute TC
        if (power){
          W <- as(listw, "CsparseMatrix")
          TC <- powerWeights(W, rho = object$rho, X = Xs, order = order, tol = tol) %*% B
        } else {
          TC <- invIrW(listw, object$rho) %*% trends
        }
      }
      
      if (pred.type == "trend") {
        res <- as.vector(trends)
      } else if(pred.type == "TC") {
        res <- as.vector(TC)
      } else if(pred.type == "BP") {
        W <- as(listw, "CsparseMatrix")
        Qss <- 1/object$s2 * (Diagonal(dim(W)[1]) - (object$rho * t(W))) %*% (Diagonal(dim(W)[1]) - (object$rho * W)) # precision matrix for LAG model
        DiagQss <- Diagonal(x = diag(Qss))
        BP <- TC - solve(DiagQss) %*% (Qss - DiagQss) %*% (ys - TC)
        # TODO: Can BP also be applied to the SEM model? Cf LeSage and Pace (2004). Note: \hat{\mu_i} need to be adapted
        res <- as.vector(BP)
      } else {
        stop("no such in-sample predictor type")
      }
      if (pred.type != "trend") attr(res, "trend") <- as.vector(trends)
    }
    attr(res, "region.id") <- as.vector(attr(ys, "names"))
  } else { # out-of-sample
    #CHECK
    if (any(row.names(newdata) %in% attr(ys, "names")) && !(pred.type == "TS" && nrow(newdata) == length(ys) && row.names(newdata) == attr(ys, "names"))) # no warning in the computation of TS in forecast
      warning("some region.id are both in data and newdata")
    if (!(pred.type == "TS" && object$type == "error" && object$etype == "error") && !(pred.type == "trend" && (object$type != "mixed" && !(object$type == "error" && object$etype == "emixed")))) { # need of listw (ie. neither in the case of defaut predictor and SEM model, nor trend type without mixed models)
      if (is.null(listw) || !inherits(listw, "listw"))
        stop ("spatial weights list required")
      if (any(! row.names(newdata) %in% attr(listw, "region.id")))
        stop("mismatch between newdata and spatial weights. newdata should have region.id as row.names")
      listw.old <- NULL
      
      # wanted order of the listw
      if (pred.type == "TS") { # only need Woo
        region.id <- row.names(newdata)
        if (!legacy.mixed) listw.old <- listw # keep the old listw to allow the computation of lagged variable from the full WX
        # listw <- mat2listw(matrix(0), row.names = row.names(newdata)) # avoid crash of subset.listw
      } else {
        region.id <- c(attr(ys, "names"), row.names(newdata))
        if (legacy.mixed) listw.old <- listw # keep the old listw to allow the computation of lagged variable from Woo (as mat2listw() cannot compute Woo from the new listw, because it has general weights lists)
      }
      
      if (length(region.id) != length(attr(listw, "region.id")) || !all(region.id == attr(listw, "region.id"))) { # if listw is not directly ok
        if (all(subset(attr(listw, "region.id"), attr(listw, "region.id") %in% region.id) == region.id) && listw$style != "M") { # only need a subset.listw, ie. spatial units are in the right order and if weights style is not unknown
          listw <- subset.listw(listw, (attr(listw, "region.id") %in% region.id), zero.policy = zero.policy)
        } else { # we use a sparse matrix transformation to reorder a listw
          if (listw$style == "M") warning("unknown weight style: a subset of listw is used without re-normalization")
          W <- as(listw, "CsparseMatrix")
          W <- W[region.id, region.id]
          style <- listw$style
          listw <- mat2listw(W, row.names = region.id, style = style) # re-normalize to keep the style
          rm(W) # avoid the use of a wrong W
        }
      }
      #optional check
      temp <- rep(NA, length(region.id))
      names(temp) <- region.id
      if (spChk && !chkIDs(temp, listw))
        stop("Check of data and weights ID integrity failed")
    }
    
    
    # DATA
    frm <- formula(object$call)
    mt <- delete.response(terms(frm, data=newdata)) # returns a terms object for the same model but with no response variable
    mf <- model.frame(mt, newdata)
    # resolved problem of missing response column in newdata reported by
    # Christine N. Meynard, 060201
    if (dim(mf)[1] != nrow(newdata))
      stop("missing values in newdata")
    Xo <- model.matrix(mt, mf)
    
    #  accommodate aliased coefficients 120314
    if (any(object$aliased))
      Xo <- Xo[,-which(object$aliased)]
    
    if (object$type == "mixed" || (object$type == "error" && object$etype == "emixed")) { # mixed model: compute WXo
      K <- ifelse(colnames(Xo)[1] == "(Intercept)", 2, 1)
      # prepare listw for computation of lagged variables
      if (legacy.mixed) { # compute WooXo
        X <- Xo
        region.id.mixed <- rownames(Xo)
      } else { # compute [WX]o
        is.not.lagged <- 1:((K-1)+(ncol(Xs)-(K-1))/2) #TODO: change the way lag variables are detected
        Xs.not.lagged <- Xs[, is.not.lagged]
        if (any(colnames(Xs.not.lagged) != colnames(Xo)))
          stop("unknown mismatch. please report this bug")
        X <- rbind(Xs.not.lagged, Xo)
        region.id.mixed <- rownames(X)
      }
      if (is.null(listw.old)) listw.mixed <- listw
      else {
        listw.mixed <- listw.old
        rm(listw.old)
      }
      if (length(region.id.mixed) != length(attr(listw.mixed, "region.id")) || !all(region.id.mixed == attr(listw.mixed, "region.id"))) { # if listw is not directly ok
        if (all(subset(attr(listw.mixed, "region.id"), attr(listw.mixed, "region.id") %in% region.id.mixed) == region.id.mixed) && listw.mixed$style != "M") { # only need a subset.listw, ie. spatial units are in the right order and weights style is not unknown
          listw.mixed <- subset.listw(listw.mixed, attr(listw.mixed, "region.id") %in% region.id.mixed, zero.policy = zero.policy)
        } else { # we use a sparse matrix transformation to reorder a listw
          if (listw$style == "M") warning("unknown weight style: a subset of listw is used without re-normalization")
          W <- as(listw.mixed, "CsparseMatrix")
          W <- W[region.id.mixed, region.id.mixed]
          style <- listw.mixed$style
          listw.mixed <- mat2listw(W, row.names = region.id.mixed, style = style) # re-normalize to keep the style
          rm(W) # avoid the use of a wrong W
        }
      }
      #optional check
      temp <- rep(NA, length(region.id.mixed))
      names(temp) <- region.id.mixed
      if (spChk && !chkIDs(temp, listw.mixed))
        stop("Check of data and weights ID integrity failed")
      
      K <- ifelse(colnames(X)[1] == "(Intercept)", 2, 1)
      m <- ncol(X)
      # check if there are enough regressors
      if (m > 1) {
        WX <- matrix(nrow=nrow(X),ncol=(m-(K-1)))
        for (k in K:m) {
          wx <- lag.listw(listw.mixed, X[,k], 
                          zero.policy=zero.policy)
          if (any(is.na(wx)))
            stop("NAs in lagged independent variable")
          WX[,(k-(K-1))] <- wx
        }
      }
      if (K == 2) {
        # unnormalized weight matrices
        if (!(listw.mixed$style == "W")) {
          intercept <- as.double(rep(1, nrow(X)))
          wx <- lag.listw(listw.mixed, intercept, 
                          zero.policy = zero.policy)
          if (m > 1) {
            WX <- cbind(wx, WX)
          } else {
            WX <- matrix(wx, nrow = nrow(X), ncol = 1)
          }
        } 
      }
      if (any(object$aliased)) {
        if (K>1 && (listw.mixed$style == "W")) colnames(WX) <- paste("lag.", colnames(X)[-1], sep="")
        else colnames(WX) <- paste("lag.", colnames(X), sep="")
        WX <- WX[,!colnames(WX) %in% names(object$aliased[object$aliased])]
      }
      if (legacy.mixed) Xo <- cbind(Xo, WX)
      else {
        WXo <- WX[(length(ys)+1):nrow(WX),]
        Xo <- cbind(Xo, WXo)
        WXs <- WX[1:length(ys),]
        Xs <- cbind(Xs.not.lagged, WXs)
      }
      rm(list = c("listw.mixed", "region.id.mixed", "X"))
    }
    trendo <- Xo %*% B
    
    if (pred.type == "TS") { # defaut predictor
      if (object$type == "error") {
        if (object$etype == "error") { # We
          signal <- rep(0, length(trendo))
          res <- trendo + signal
          attr(res, "trend") <- trendo
          attr(res, "signal") <- signal
        } else if (object$etype == "emixed") { # WX + We
          signal <- rep(0, length(trendo))
          res <- trendo + signal
          attr(res, "trend") <- trendo
          attr(res, "signal") <- signal
        } else stop("unkown error model etype")
      } else if (object$type == "mixed") { # Wy+WX
        if (power) {
          W <- as(listw, "CsparseMatrix")
          res <- c(as(powerWeights(W, rho=object$rho,
                                   X=trendo, order=order, tol=tol), "matrix"))
        } else {
          res <- c(invIrW(listw, object$rho) %*% trendo)
        }
        if (legacy) {
          signal <- object$rho * lag.listw(listw, 
                                           res, zero.policy=zero.policy)
          res <- c(trendo + signal)
        } else {
          signal <- res - trendo
        }
        #                        if (pred.se) {
        #                            samples <- attr(lagImpact, "samples")$samples
        #                            irho <- attr(lagImpact, "samples")$irho
        #                            drop2beta <- attr(lagImpact, "samples")$drop2beta
        #                            nSim <- nrow(samples)
        #                            outmat <- matrix(NA, ncol=nSim, nrow=nrow(X))
        #                            for (i in 1:nSim) {
        #                                B <- samples[i, -drop2beta]
        #                                trend <- X %*% B
        #                                rho <- samples[i, irho]
        #                                if (power) {
        #                                    res <- c(as(powerWeights(W, rho=rho,
        #                                    X=trend, order=order, tol=tol), "matrix"))
        #                                } else {
        #                                    res <- c(invIrW(listw, rho) %*% trend)
        #                                }
        #                                outmat[,i] <- res
        #                            }
        #                            pred.se <- apply(outmat, 1, sd)
        #                            attr(res, "pred.se") <- pred.se
        #                        }
        attr(res, "trend") <- c(trendo)
        attr(res, "signal") <- c(signal)
      } else { # Wy
        if (power) {
          W <- as(listw, "CsparseMatrix")
          res <- c(as(powerWeights(W, rho=object$rho,
                                   X=trendo, order=order, tol=tol), "matrix"))
        } else {
          res <- c(invIrW(listw, object$rho) %*% trendo)
        }
        if (legacy) {
          signal <- object$rho * lag.listw(listw, 
                                           res, zero.policy=zero.policy)
          res <- c(trendo + signal)
        } else {
          signal <- res - trendo
        }
        #                        if (pred.se) {
        #                            samples <- attr(lagImpact, "samples")$samples
        #                            irho <- attr(lagImpact, "samples")$irho
        #                            drop2beta <- attr(lagImpact, "samples")$drop2beta
        #                            nSim <- nrow(samples)
        #                            outmat <- matrix(NA, ncol=nSim, nrow=nrow(X))
        #                            for (i in 1:nSim) {
        #                                B <- samples[i, -drop2beta]
        #                                trend <- X %*% B
        #                                rho <- samples[i, irho]
        #                                if (power) {
        #                                    res <- c(as(powerWeights(W, rho=rho,
        #                                    X=trend, order=order, tol=tol), "matrix"))
        #                                } else {
        #                                    res <- c(invIrW(listw, rho) %*% trend)
        #                                }
        #                                outmat[,i] <- res
        #                            }
        #                            pred.se <- apply(outmat, 1, sd)
        #                            attr(res, "pred.se") <- pred.se
        #                        }
        attr(res, "trend") <- c(trendo)
        attr(res, "signal") <- c(signal)
      }
    } else { # new predictors
      if (pred.type %in% c("TS1", "KP4", "KP2", "KP3")) { # need to compute TS1/KP4
        Wos <- .listw.decompose(listw, region.id.data = attr(ys, "names"), region.id.newdata = row.names(newdata), type = "Wos")$Wos
        if (is.null(object$rho)) TS1 <- Xo %*% B
        else TS1 <- Xo %*% B + object$rho * Wos %*% ys
      }
      if (pred.type %in% c("TC", "BP", "BPN")) { # need to compute TC
        #notations of C.Thomas and al (2015)
        if (all.data | pred.type %in% c("BP", "BPN")) {
          # compute s and o units together
          X <- rbind(Xs, Xo)
          trend <- X %*% B
          if (power){
            W <- as(listw, "CsparseMatrix")
            TC <- powerWeights(W, rho = object$rho, X = X, order = order, tol = tol) %*% B
          } else {
            TC <- invIrW(listw, object$rho) %*% trend
          }
        } else { # TCo = TC for out-of-sample spatial units
          listw.d <- .listw.decompose(listw, region.id.data = attr(ys, "names"), region.id.newdata = row.names(newdata), type = c("Wss", "Wos", "Wso", "Woo"))
          Wss <- listw.d$Wss
          Wso <- listw.d$Wso
          Wos <- listw.d$Wos
          Woo <- listw.d$Woo
          rm(listw.d)
          mB <- - object$rho * Wso
          mC <- - object$rho * Wos
          mD <- Diagonal(dim(Woo)[1]) - object$rho * Woo
          if (power){
            mAInvXsB <- powerWeights(Wss, rho = object$rho, X = Xs, order = order, tol = tol) %*% B
            mAInvmB <- powerWeights(Wss, rho = object$rho, X = mB, order = order, tol = tol)
            E <- solve(mD - mC %*% mAInvmB)
            TCo <- - E %*% mC %*% mAInvXsB + E %*% Xo %*% B
          } else {
            #mA <- Diagonal(dim(Wss)[1]) - object$rho * Wss
            #mAInv <- solve(mA)
            mAInv <- invIrW(Wss, object$rho)
            E <- solve(mD - mC %*% mAInv %*% mB)
            TCo <- - E %*% mC %*% mAInv %*% Xs %*% B + E %*% Xo %*% B
          }
        }
      }
      
      if (pred.type == "trend") {
        res <- as.vector(trendo)
      } else if (pred.type %in% c("TS1", "KP4")) {
        res <- as.vector(TS1)
      } else if (pred.type == "TC") {
        if(all.data) {
          res <- as.vector(TC)
        } else {
          res <- as.vector(TCo)
        }
      } else if (pred.type == "BP") {
        is.data <- 1:length(ys)
        is.newdata <- (length(ys)+1):length(TC)
        TCo <- TC[is.newdata]
        TCs <- TC[is.data]
        W <- as(listw, "CsparseMatrix")
        Q <- 1/object$s2 * ( Diagonal(dim(W)[1]) - object$rho * (t(W) + W) + object$rho^2 * (t(W) %*% W) )
        Qoo <- Q[is.newdata, is.newdata]
        Qos <- Q[is.newdata, is.data]
        BPo <- TCo - solve(Qoo) %*% Qos %*% (ys - TCs)
        res <- as.vector(BPo)
      } else if (pred.type == "BPW") {
        if (power){
          W <- as(listw, "CsparseMatrix")
          invW <- powerWeights(W, rho = object$rho, X = Diagonal(dim(W)[1]), order = order, tol = tol)
        } else {
          invW <- invIrW(listw, object$rho)
        }
        X <- rbind(Xs, Xo)
        TC <- invW %*% X %*% B
        is.data <- 1:length(ys)
        is.newdata <- (length(ys)+1):length(TC)
        TCo <- TC[is.newdata]
        TCs <- TC[is.data]
        #Sigma <- object$s2 * solve((Diagonal(dim(W)[1]) - object$rho * t(W)) %*% (Diagonal(dim(W)[1]) - object$rho * W))
        Sigma <- object$s2 * invW %*% t(invW)
        Sos <- Sigma[is.newdata, is.data, drop=F]
        Sss <- Sigma[is.data, is.data, drop=F]
        Wos <- .listw.decompose(listw, region.id.data = attr(ys, "names"), region.id.newdata = row.names(newdata), type = "Wos")$Wos
        BPW <- as.numeric(TCo + Sos %*% t(Wos) %*% solve(Wos %*% Sss %*% t(Wos)) %*% (Wos %*% ys - Wos %*% TCs))
        res <- as.vector(BPW)
      } else if (pred.type == "BPN") {
        is.data <- 1:length(ys)
        is.newdata <- (length(ys)+1):length(TC)
        TCs <- TC[is.data]
        TCo <- TC[is.newdata]
        # compute J = set of all sites in S which are neighbors of at least one site in O
        O <- which(attr(listw,"region.id") %in% row.names(newdata))
        S <- which(attr(listw,"region.id") %in% attr(ys, "names"))
        #bug fix: J = set of S s.t. at least 1 elmt of O is neigbourg with S. And not: set of S s.t. these elmts of S are neigbourgs with at least 1 elmt of S.
        #different iff W is not symetric
        #J.logical <- rep(FALSE, length(listw$neighbours))
        #for (i in S) {
        #  J.logical[i] <- any(O %in% listw$neighbours[[i]])
        #}
        #J <- attr(listw,"region.id")[J.logical]
        J.ref <- NULL
        for (i in O) {
          J.ref <- c(J.ref, listw$neighbours[[i]][listw$neighbours[[i]] %in% S])
        }
        J <- attr(listw,"region.id")[unique(J.ref)]
        J <- attr(listw,"region.id")[attr(listw,"region.id") %in% J] # keep the order of listw
        if (length(J)<1) {
          warning("out-of-sample units have no neighbours")
          BPN <- TCo
        } else {
          W <- as(listw, "CsparseMatrix")
          region.id <- c(J, row.names(newdata))
          W_jo <- W[region.id, region.id]
          rm(W)
          Q_jo <- 1/object$s2 * (Diagonal(length(region.id)) - object$rho * (W_jo + t(W_jo)) + object$rho^2 * (t(W_jo) %*% W_jo))
          is.j <- 1:length(J)
          is.o <- (length(J)+1):length(region.id)
          Qoo <- Q_jo[is.o, is.o]
          Qoj <- Q_jo[is.o, is.j]
          rm(Q_jo)
          yj <- ys[J]
          TCj <- TCs[attr(ys, "names") %in% J]
          BPN <- as.vector(TCo - solve(Qoo) %*% Qoj %*% (yj - TCj))
        }
        res <- as.vector(BPN)
      } else if (pred.type %in% c("TC1", "KP1")) {
        if (nrow(newdata) > 1)
          warning("newdata have more than 1 row and the predictor type is leave-one-out")
        region.id.data <- attr(ys, "names")
        region.id.newdata <- row.names(newdata)
        res <- rep(NA, nrow(newdata))
        W <- as(listw, "CsparseMatrix")
        style <- listw$style
        if (listw$style == "M")
          warning("unknown weight style: a subset of listw is used without re-normalization")
        for (i in 1:nrow(newdata)) {
          region.id.temp <- c(region.id.data, region.id.newdata[i])
          Wi <- W[region.id.temp, region.id.temp]
          listwi <- mat2listw(Wi, row.names = region.id.temp, style = style) # re-normalize
          if (power)
            Wi <- as(listwi, "CsparseMatrix")
          Xi <- rbind(Xs, Xo[i,])
          if (power) {
            res[i] <- c(as(powerWeights(Wi, rho=object$rho, X=Xi, order=order, tol=tol), "matrix") %*% B)[length(region.id.temp)]
          } else {
            trendi <- c(trends, trendo[i])
            res[i] <- (invIrW(listwi, object$rho) %*% trendi)[length(region.id.temp)]
          }
        }
      } else if (pred.type == "BP1") {
        if (nrow(newdata) > 1)
          warning("newdata have more than 1 row and the predictor type is leave-one-out")
        region.id.data <- attr(ys, "names")
        region.id.newdata <- row.names(newdata)
        BP1o <- rep(NA, nrow(newdata))
        W <- as(listw, "CsparseMatrix")
        style <- listw$style
        if (listw$style == "M")
          warning("unknown weight style: a subset of listw is used without re-normalization")
        for (i in 1:nrow(newdata)) {
          region.id.temp <- c(region.id.data, region.id.newdata[i])
          Wi <- W[region.id.temp, region.id.temp]
          listwi <- mat2listw(Wi, row.names = region.id.temp, style = style) # re-normalize
          Wi <- as(listwi, "CsparseMatrix")
          Xi <- rbind(Xs, Xo[i,])
          # compute TC1 for S and o units
          TC1i <- rep(NA, length(region.id.temp))
          is.data <- 1:length(ys)
          is.newdata <- length(region.id.temp)
          if (power) {
            TC1i <- c(as(powerWeights(Wi, rho=object$rho, X=Xi, order=order, tol=tol), "matrix") %*% B)
          } else {
            trendi <- c(trends, trendo[i])
            TC1i <- (invIrW(Wi, object$rho) %*% trendi)
          }
          TC1si <- TC1i[is.data]
          TC1oi <- TC1i[is.newdata]
          
          Qi <- 1/object$s2 * ( Diagonal(dim(Wi)[1]) - object$rho * (t(Wi) + Wi) + object$rho^2 * (t(Wi) %*% Wi) )
          Qooi <- Qi[is.newdata, is.newdata]
          Qosi <- Qi[is.newdata, is.data]
          BP1o[i] <- TC1oi - solve(Qooi) %*% Qosi %*% (ys - TC1si)
        }
        res <- as.vector(BP1o)
      } else if (pred.type == "BPW1") {
        if (nrow(newdata) > 1)
          warning("newdata have more than 1 row and the predictor type is leave-one-out")
        region.id.data <- attr(ys, "names")
        region.id.newdata <- row.names(newdata)
        BPW1o <- rep(NA, nrow(newdata))
        W <- as(listw, "CsparseMatrix")
        style <- listw$style
        if (listw$style == "M")
          warning("unknown weight style: a subset of listw is used without re-normalization")
        for (i in 1:nrow(newdata)) {
          region.id.temp <- c(region.id.data, region.id.newdata[i])
          Wi <- W[region.id.temp, region.id.temp]
          listwi <- mat2listw(Wi, row.names = region.id.temp, style = style) # re-normalize
          Wi <- as(listwi, "CsparseMatrix")
          Xi <- rbind(Xs, Xo[i,])
          is.data <- 1:length(ys)
          is.newdata <- length(region.id.temp)
          if (power){
            invWi <- powerWeights(Wi, rho = object$rho, X = Diagonal(dim(Wi)[1]), order = order, tol = tol)
          } else {
            invWi <- invIrW(Wi, object$rho)
          }
          TC1i <- invWi %*% Xi %*% B
          TC1oi <- TC1i[is.newdata]
          TC1si <- TC1i[is.data]
          #Sigma <- object$s2 * solve((Diagonal(dim(W)[1]) - object$rho * t(W)) %*% (Diagonal(dim(W)[1]) - object$rho * W))
          Sigmai <- object$s2 * invWi %*% t(invWi)
          Sosi <- Sigmai[is.newdata, is.data, drop=F]
          Sssi <- Sigmai[is.data, is.data]
          Wosi <- Wi[region.id.newdata[i], region.id.data, drop=F]
          BPW1o[i] <- as.numeric(TC1oi + Sosi %*% t(Wosi) %*% solve(Wosi %*% Sssi %*% t(Wosi)) %*% (Wosi %*% ys - Wosi %*% TC1si))
        }
        res <- as.vector(BPW1o)
      } else if (pred.type == "BPN1") {
        if (nrow(newdata) > 1)
          warning("newdata have more than 1 row and the predictor type is leave-one-out")
        region.id.data <- attr(ys, "names")
        region.id.newdata <- row.names(newdata)
        BPN1o <- rep(NA, nrow(newdata))
        W <- as(listw, "CsparseMatrix")
        style <- listw$style
        if (listw$style == "M")
          warning("unknown weight style: a subset of listw is used without re-normalization")
        for (i in 1:nrow(newdata)) {
          region.id.temp <- c(region.id.data, region.id.newdata[i])
          Wi <- W[region.id.temp, region.id.temp]
          listwi <- mat2listw(Wi, row.names = region.id.temp, style = style) # re-normalize
          Wi <- as(listwi, "CsparseMatrix")
          Xi <- rbind(Xs, Xo[i,])
          # compute TC1 for S and o units
          TC1i <- rep(NA, length(region.id.temp))
          is.data <- 1:length(ys)
          is.newdata <- length(region.id.temp)
          if (power) {
            TC1i <- c(as(powerWeights(Wi, rho=object$rho, X=Xi, order=order, tol=tol), "matrix") %*% B)
          } else {
            trendi <- c(trends, trendo[i])
            TC1i <- (invIrW(Wi, object$rho) %*% trendi)
          }
          TC1si <- TC1i[is.data]
          TC1oi <- TC1i[is.newdata]
          
          # compute J = set of all sites in S which are neighbors of at least one site in O
          o <- which(attr(listwi,"region.id") %in% row.names(newdata)[i])
          S <- which(attr(listwi,"region.id") %in% attr(ys, "names"))
          J.ref <- NULL
          for (k in o) {
            J.ref <- c(J.ref, listwi$neighbours[[k]][listwi$neighbours[[k]] %in% S])
          }
          J <- attr(listwi,"region.id")[unique(J.ref)]
          J <- attr(listwi,"region.id")[attr(listwi,"region.id") %in% J] # keep the order of listw
          if (length(J)<1) {
            warning("out-of-sample units have no neighbours")
            BPN1o[i] <- TC1oi
          } else {
            region.id <- c(J, row.names(newdata)[i])
            W_jo <- Wi[region.id, region.id, drop=F]
            Q_jo <- 1/object$s2 * (Diagonal(length(region.id)) - object$rho * (W_jo + t(W_jo)) + object$rho^2 * (t(W_jo) %*% W_jo))
            is.j <- 1:length(J)
            is.o <- (length(J)+1):length(region.id)
            Qoo <- Q_jo[is.o, is.o, drop=F]
            Qoj <- Q_jo[is.o, is.j, drop=F]
            yj <- ys[J]
            TC1j <- TC1si[attr(ys, "names") %in% J]
            BPN1o[i] <- as.vector(TC1oi - solve(Qoo) %*% Qoj %*% (yj - TC1j))
          }
        }
        res <- as.vector(BPN1o)
      } else if (pred.type == "KP2") {
        if (nrow(newdata) > 1)
          warning("newdata have more than 1 row and the predictor type is leave-one-out")
        region.id.data <- attr(ys, "names")
        region.id.newdata <- row.names(newdata)
        W <- as(listw, "CsparseMatrix")
        KP2 <- rep(NA, nrow(newdata))
        for (i in 1:nrow(newdata)) {
          region.id.temp <- c(region.id.data, region.id.newdata[i])
          Wi <- W[region.id.temp, region.id.temp]
          Xi <- rbind(Xs, Xo[i,])
          wi <- Wi[length(region.id.temp), ]
          yi <- c(ys, 0)
          if (object$type %in% c("sac", "sacmixed")) { # compute GL, GR, sum.u, sum.y
            if (power){
              GL <- powerWeights(Wi, rho= object$lambda, order= order, tol= tol, X= Diagonal(length(ys)+1))
              GR <- powerWeights(Wi, rho= object$rho, order= order, tol= tol, X= Diagonal(length(ys)+1))
            } else {
              GL <- invIrW(Wi, object$lambda)
              GR <- invIrW(Wi, object$rho)
            }
            sum.u <- GL %*% t(GL)
            sum.y <- GR %*% sum.u %*% t(GR)
          } else if (object$type %in% c("lag", "mixed")) {
            GL <- Diagonal(length(ys)+1)
            if (power){
              GR <- powerWeights(Wi, rho= object$rho, order= order, tol= tol, X= Diagonal(length(ys)+1))
            } else {
              GR <- invIrW(Wi, object$rho)
            }
            sum.u <- Diagonal(length(ys)+1)
            sum.y <- GR %*% t(GR)
          } else if (object$type == "error") {
            GR <- Diagonal(length(ys)+1)
            if (power){
              GL <- powerWeights(Wi, rho= object$lambda, order= order, tol= tol, X= Diagonal(length(ys)+1))
            } else {
              GL <- invIrW(Wi, object$lambda)
            }
            sum.u <- GL %*% t(GL)
            sum.y <- sum.u
          } else stop("unknown model type")
          covar <- as.vector(sum.u[length(region.id.temp), ] %*% t(GR) %*% wi / (wi %*% sum.y %*% wi))
          Ewiy <- as.vector(wi %*% GR %*% Xi %*% B)
          KP2[i] <- TS1[i] + covar %*% (wi %*% yi - Ewiy)
        }
        res <- as.vector(KP2)
      } else if (pred.type == "KP3") {
        if (nrow(newdata) > 1)
          warning("newdata have more than 1 row and the predictor type is leave-one-out")
        region.id.data <- attr(ys, "names")
        region.id.newdata <- row.names(newdata)
        W <- as(listw, "CsparseMatrix")
        KP3 <- rep(NA, nrow(newdata))
        for (i in 1:nrow(newdata)) {
          region.id.temp <- c(region.id.data, region.id.newdata[i])
          Wi <- W[region.id.temp, region.id.temp]
          Xi <- rbind(Xs, Xo[i,])
          wi <- Wi[length(region.id.temp), ]
          yi <- c(ys, 0)
          if (object$type %in% c("sac", "sacmixed")) { # compute GL, GR, sum.u, sum.y
            if (power){
              GL <- powerWeights(Wi, rho= object$lambda, order= order, tol= tol, X= Diagonal(length(ys)+1))
              GR <- powerWeights(Wi, rho= object$rho, order= order, tol= tol, X= Diagonal(length(ys)+1))
            } else {
              GL <- invIrW(Wi, object$lambda)
              GR <- invIrW(Wi, object$rho)
            }
            sum.u <- GL %*% t(GL)
            sum.y <- GR %*% sum.u %*% t(GR)
          } else if (object$type %in% c("lag", "mixed")) {
            GL <- Diagonal(length(ys)+1)
            if (power){
              GR <- powerWeights(Wi, rho= object$rho, order= order, tol= tol, X= Diagonal(length(ys)+1))
            } else {
              GR <- invIrW(Wi, object$rho)
            }
            sum.u <- Diagonal(length(ys)+1)
            sum.y <- GR %*% t(GR)
          } else if (object$type == "error") {
            GR <- Diagonal(length(ys)+1)
            if (power){
              GL <- powerWeights(Wi, rho= object$lambda, order= order, tol= tol, X= Diagonal(length(ys)+1))
            } else {
              GL <- invIrW(Wi, object$lambda)
            }
            sum.u <- GL %*% t(GL)
            sum.y <- sum.u
          } else stop("unknown model type")
          cov <- sum.u[length(region.id.temp), ] %*% t(GR)[, -length(region.id.temp)]
          KP3[i] <- as.vector(TS1[i] + cov %*% solve(sum.y[-length(region.id.temp), -length(region.id.temp)]) %*% (ys - GR[-length(region.id.temp),] %*% Xi %*% B))
        }
        res <- as.vector(KP3)
      } else if (pred.type == "KP5") {
        if (nrow(newdata) > 1)
          warning("newdata have more than 1 row and the predictor type is leave-one-out")
        Wos <- .listw.decompose(listw, region.id.data = attr(ys, "names"), region.id.newdata = row.names(newdata), type = "Wos")$Wos
        res <- as.vector(trendo + object$lambda * Wos %*% (ys - trends))
      } else {
        stop("unknow predictor type")
      }
    }
    # add region.id attribute
    if (length(res) == nrow(newdata)) {
      attr(res, "region.id") <- as.vector(row.names(newdata))
    } else if (length(res) == length(ys)+nrow(newdata)) {
      attr(res, "region.id") <- c(attr(ys, "names"), row.names(newdata))
    } else stop("incorrect final output")
  }
  attr(res, "pred.type") <- pred.type
  attr(res, "call") <- match.call()
  class(res) <- "sarlm.pred"
  res
}

# decompose a listw object into Wss Wso Wos and Woo sparse matrices
.listw.decompose <- function(listw, region.id.data, region.id.newdata, type = c("Wss", "Wos", "Wso", "Woo")) { # TODO: hidden? in this file? zero.policy?
  if (is.null(listw) || !inherits(listw, "listw")) 
    stop ("spatial weights list required")
  region.id <- attr(listw, "region.id")
  if (!all(region.id.data %in% region.id))
    stop("at least one region.id in data is not in listw object")
  if (!all(region.id.newdata %in% region.id))
    stop("at least one region.id in newdata is not in listw object")
  if (!all(type %in% c("Wss", "Wos", "Wso", "Woo")))
    stop("type is incorrect")
  W <- as(listw, "CsparseMatrix")
  s <- list(Wss = NULL, Wos = NULL, Wso = NULL, Woo = NULL)
  if ("Wss" %in% type)
    s$Wss <- W[region.id.data, region.id.data, drop=F]
  if ("Wos" %in% type)
    s$Wos <- W[region.id.newdata, region.id.data, drop=F]
  if ("Wso" %in% type)
    s$Wso <- W[region.id.data, region.id.newdata, drop=F]
  if ("Woo" %in% type)
    s$Woo <- W[region.id.newdata, region.id.newdata, drop=F]
  return(s)
}


print.sarlm.pred <- function(x, ...) {
  res <- as.data.frame(x)
  print(res, ...)
  invisible(res)
}


as.data.frame.sarlm.pred <- function(x, ...) {
  #    res <- data.frame(fit=as.vector(x), trend=attr(x, "trend"), 
  #        signal=attr(x, "signal"))
  #fix bug when no signal or trend attributes
  res <- data.frame(fit=as.vector(x))
  if(!is.null(attr(x, "region.id"))) row.names(res) <- attr(x, "region.id")
  if(!is.null(attr(x, "trend"))) res$trend <- attr(x, "trend")
  if(!is.null(attr(x, "signal"))) res$signal <- attr(x, "signal")
  res
}

