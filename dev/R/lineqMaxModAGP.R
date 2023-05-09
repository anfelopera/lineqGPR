#' @title MaxMod algorithm for \code{"lineqGP"} and \code{"lineqAGP"} Models
#' @description A wrapper function for the maximum Modification algorithm.
#' 
#' @param model an object with class \code{lineqGP} or \code{lineqAGP}
#' @param xtest test data for assessing the modification of the MAP estimate
#' @param tol a number corresponding to the tolerance of algorithm.
#' The algorithm stops if the (MaxMod criterion)/norm(pred_old) < tol
#' @param max_iter an integer corresponding to number of iterations
#' @param reward_new_knot a number corresponding to the reward of adding a new knot in an existing dimension
#' @param reward_new_dim a number corresponding to the reward of adding a new dimension
#' @param print_iter a logical variable to print results at each iteration
#' @param nClusters an integer corresponding to the number of clusters
#' @param save_history a logical variable to save the model at each iteration
# #' @param MCtest a logical variable to approximate MaxMod criterion via MC
#' 
#' @return an object with class \code{lineqGP} or \code{lineqAGP} containing the resulting model
#'
#' @author A. F. Lopez-Lopera
#' 
#' @export
AdditiveMaxMod <- function(model,
                           xtest,
                           tol = 1e-4,
                           max_iter = 10*ncol(model$x),
                           reward_new_knot = tol,
                           reward_new_dim = 1e-9,
                           print_iter = FALSE,
                           nClusters = 1,
                           save_history = FALSE #, MCtest = FALSE,
                           ) {
  
  D <- ncol(model$x) # input dimension
  
  activeDim <- rep(FALSE, D)
  activeDim_IdxSeq <- c()
  MaxMod_values <- c()
  optDecision <- c()
  ulist <- c()
  pred <- c()
  relMaxMod_values <- c()
  relativeMaxMod_criterion <- Inf

  if (save_history)
    models_hist <- list()
  
  MaxMod_iter_values <- matrix(NaN, 3, D)
  rownames(MaxMod_iter_values) <- c("MaxMod criterion", "knot's position", "decision")
  colnames(MaxMod_iter_values) <- paste("dim", 1:D)
  
  if (nClusters > 1) {
    if (!requireNamespace("foreach", quietly = TRUE))
      stop("Package \"foreach\" not found")
    if (!requireNamespace("doParallel", quietly = TRUE))
      stop("Package \"doParallel\" not found")
    
    nClusters_max <- max(parallel::detectCores(), 1)
    if (nClusters > nClusters_max) {
      warning("The maximum number of clusters has been exceeded.
      Further computations will consider only ", nClusters_max, " clusters")
      nClusters <- nClusters_max
    }
    
    cl <- parallel::makeCluster(nClusters)
    doParallel::registerDoParallel(cl)
  }
  
  
  for (k in 1:max_iter) {
    MaxMod_temp_values <- MaxMod_iter_values
    if (nClusters > 1) {
      MaxMod_temp_values[1:2, ] <- try(foreach::"%dopar%"(foreach::foreach(i=1:D, .combine=cbind,
                                                                           .errorhandling='remove'), {
                                                                             MaxModCriterion(model, ulist, 
                                                                                             xtest, activeDim, 
                                                                                             activeDim_IdxSeq, i,
                                                                                             reward_new_dim = reward_new_dim,
                                                                                             reward_new_knot = reward_new_knot, 
                                                                                             iter = k, pred = pred #, MCtest = MCtest,
                                                                                             )
                                                                           }))
      
    } else {
      for (i in 1:D) {
        MaxMod_temp_values[1:2, i] <- MaxModCriterion(model, ulist, 
                                                      xtest, activeDim,
                                                      activeDim_IdxSeq, idx_add = i,
                                                      reward_new_dim = reward_new_dim,
                                                      reward_new_knot = reward_new_knot,
                                                      iter = k, pred = pred #, MCtest = MCtest
                                                      ) # pred can be dropped !
      }
    }
    
    MaxMod_temp_values[3,] <- 0
    idxOpt <- try(which.min(MaxMod_temp_values[1,]))
    if(length(idxOpt) == 0) # a potential error from doParallel/foreach
      stop("A fatal error occured in the parallel loop")
    
    MaxMod_temp_values[3,idxOpt] <- 1L
    MaxMod_temp_values[1,] <- -MaxMod_temp_values[1,]
    
    if(print_iter)
      print(MaxMod_temp_values)
    
    if (isFALSE(activeDim[idxOpt])) {
      activeDim[idxOpt] <- TRUE
      flag_AddDim <- TRUE
      activeDim_IdxSeq <- c(activeDim_IdxSeq, idxOpt)
      idxOptOrdered <- which(activeDim_IdxSeq == idxOpt)
    } else {
      flag_AddDim <- FALSE
    }
    
    model_temp <-  create(class = 'lineqAGP',
                          x = model$x[, activeDim_IdxSeq], y = model$y,
                          constrType = model$constrType[activeDim_IdxSeq])
    model_temp$varnoise <- model$varnoise
    model_temp$nugget <- model$nugget
    model_temp$kernParam <- model$kernParam[activeDim_IdxSeq]
    
    if (isTRUE(flag_AddDim)) {
      if (k == 1) {
        model_temp$ulist[[1]] <- c(0, 1)
      } else {
        model_temp$ulist[-idxOptOrdered] <- ulist
        model_temp$ulist[[idxOptOrdered]] <- c(0, 1)
      }
    } else {
      idxOptOrdered <- which(activeDim_IdxSeq == idxOpt)
      unew <- MaxMod_temp_values[2, idxOpt]
      model_temp$ulist <- ulist
      model_temp$ulist[[idxOptOrdered]] <- sort(c(ulist[[idxOptOrdered]], unew))
    }
    
    model_temp <- lineqGPOptim(model_temp,
                               additive = TRUE,
                               estim.varnoise = TRUE, # to add this info at the MaxMod level
                               bounds.varnoise = c(1e-7, Inf), # to add this info at the MaxMod level
                               lb = rep(1e-2, model_temp$d+1), ub = c(Inf, rep(1, model_temp$d)) # to add this info at the MaxMod level
    )
    pred_temp <- predict(model_temp, xtest[, activeDim_IdxSeq])
    
    idxOptSorted <- which(activeDim_IdxSeq == idxOpt)
    MaxMod_criterion <- MaxMod_temp_values[1, idxOpt]
    
    if (k >= 2) {
      z <- pred$xi.map
      Psi <- GramMatrixPhiAdditive(model_maxmod$ulist)
      Psi_vec <- GramVectorPhiAdditive(model_maxmod$ulist)
      zPsivec <- sapply(1:(length(z)), function(x) t(z[[x]]) %*% Psi_vec[[x]])
      zPsivec2 <- sapply(1:(length(z)), function(x) zPsivec[[x]]^2)
      zPsiz <- sapply(1:(length(z)), function(x) t(z[[x]]) %*% Psi[[x]] %*% z[[x]])
      
      varfOld <- sum(zPsiz) - sum(zPsivec2) + sum(zPsivec)^2
      varfOld <- max(0, varfOld)
      
      relativeMaxMod_criterion <- MaxMod_criterion/varfOld
    }
    
    if (print_iter) {
      if (!is.na(MaxMod_temp_values[2,idxOpt])) {
        message("Iter ", k, " - New knot added: d = ", idxOptSorted,
                " (Relative MaxMod criterion = ", relativeMaxMod_criterion, ")")
      } else {
        message("Iter ", k, " - dimension ", idxOpt, " added as d = ", length(activeDim_IdxSeq),
                " (Relative MaxMod criterion = ", relativeMaxMod_criterion, ")")
      }
    }
    
    if (k > 1 & relativeMaxMod_criterion < tol) {
      message("The sequential algorithm converged")
      break
    } else {
      relMaxMod_values <- c(relMaxMod_values, relativeMaxMod_criterion)
      MaxMod_values <- c(MaxMod_values, MaxMod_criterion)
      model_maxmod <- model_temp
      pred <- pred_temp

      optDecision <- c(optDecision, idxOpt)
      ulist <- model_maxmod$ulist
      pred <- predict(model_maxmod, xtest[, activeDim_IdxSeq])
      
      model_maxmod$kernParam$l <- NULL
      model$kernParam[activeDim_IdxSeq] <- model_maxmod$kernParam
      if (save_history)
        models_hist[[k]] <- model_maxmod
    }
    if (k == max_iter) 
      message("Run out of budget")
  }
  
  if (nClusters > 1)
    parallel::stopCluster(cl)
  
  model_maxmod$MaxMod$optDecision <- optDecision
  model_maxmod$MaxMod$MaxMod_values <- MaxMod_values
  model_maxmod$MaxMod$relMaxMod_values <- relMaxMod_values
  model_maxmod$MaxMod$models_hist <- models_hist
  return(model_maxmod) 
}

#' @title MaxMod criterion for \code{"lineqAGP"} Models
#' @description MaxMod criterion used to add a new knot or a new active dimension (see Bachoc et al., 2020)
#' 
#' @param model an object with class \code{lineqAGP}
#' @param ulist a list with the location of the knots
#' @param xtest test data for assessing the modification of the MAP estimate
#' @param activeDim a sequence containing the active dimensions
#' @param activeDim_IdxSeq a sequence containing the indices of the actived input dimensions 
#' @param idx_add index of the dimension to be activated
#' @param reward_new_dim a number corresponding to the reward of adding a new dimension
#' @param reward_new_knot a number corresponding to the reward of adding a new knot in an existing dimension
#' @param iter an integer corresponding to the iteration of the MaxMod algorithm 
#' @param pred an object with class \code{lineqAGP} containing the predictive model
# #' @param MCtest a logical variable to approximate MaxMod criterion via MC
#' 
#' @return the value of the knot that minimize the MaxMod criterion and the 
#' value of the objective function 
#'
#' @author A. F. Lopez-Lopera
#'
#' @references F. Bachoc, A. F. Lopez-Lopera, and O. Roustant (2022),
#' "Sequential construction and dimension reduction of additive Gaussian processes under inequality constraints".
# #' \emph{ArXiv e-prints}
# #' <arXiv:2009.04188>
#'
#' @importFrom utils tail
#' @export
MaxModCriterion.lineqAGP <- function(model, ulist, xtest,
                                     activeDim, activeDim_IdxSeq, idx_add,
                                     reward_new_dim, reward_new_knot, iter, pred #, MCtest,
                                     ) {
  flag_AddDim <- isFALSE(activeDim[idx_add])
  
  # adding an active dimension
  if (flag_AddDim) {
    activeDim[idx_add] <- TRUE
    activeDim_IdxSeq <- c(activeDim_IdxSeq, idx_add)
    
    # creating a new model with the new input variable
    model_update <-  create(class = 'lineqAGP',
                            x = model$x[, activeDim_IdxSeq],
                            y = model$y,
                            constrType = model$constrType[activeDim_IdxSeq])
    model_update$nugget <- model$nugget
    model_update$varnoise <- model$varnoise
    
    idxAddOrdered <- which(activeDim_IdxSeq == idx_add)
    if (iter == 1) {
      model_update$ulist[[1]] <- c(0, 1)
      model_update$constrType <- "none" # to be defined as an unconstrained problem
    } else {
      model_update$ulist[-idxAddOrdered] <- ulist
      model_update$ulist[[idxAddOrdered]] <- c(0, 1)
      constrIdx <- which(sapply(model_update$ulist, length) != 2)
      if (length(constrIdx) == 0L) { # to be defined as an unconstrained problem
        model_update$constrType <- rep("none", model_update$d)
      } else {
        model_update$constrIdx <- constrIdx
      }
    }
    model_update$localParam$m <- sapply(model_update$ulist, length)
    model_update$varnoise <- model$varnoise
    model_update$kernParam$nugget <- model$kernParam$nugget
    model_update$kernParam <- model$kernParam[activeDim_IdxSeq]
    
    
    # optimizing the constrained model
    model_update <- lineqGPOptim(model_update,
                                 additive = TRUE,
                                 estim.varnoise = TRUE, # to add this info at the MaxMod level
                                 bounds.varnoise = c(1e-7, Inf), # to add this info at the MaxMod level
                                 lb = rep(1e-2, model_update$d+1), ub = c(Inf, rep(10, model_update$d))) # to add this info at the MaxMod level

    pred_update <- predict(model_update, xtest[, activeDim_IdxSeq])
    beta_update <- pred_update$xi.map
    
    if (iter == 1) {
      # z <- beta_update
      # Psi <- GramMatrixPhiAdditive(model_update$ulist)
      # Psi_vec <- GramVectorPhiAdditive(model_update$ulist)
      # zPsivec <- sapply(1:(length(z)), function(x) t(z[[x]]) %*% Psi_vec[[x]])
      # zPsivec2 <- sapply(1:(length(z)), function(x) zPsivec[[x]]^2)
      # zPsiz <- sapply(1:(length(z)), function(x) t(z[[x]]) %*% Psi[[x]] %*% z[[x]])
      # 
      # f <- sum(zPsiz) - sum(zPsivec2) #+ sum(zPsivec)^2
      # f <- max(0, f)
      
      f <- diff(unlist(beta_update))^2/12 + sum(unlist(beta_update))^2/4
    } else {
      beta <- pred$xi.map
      beta[[length(beta_update)]] <- c(0,0)
      
      z <- lapply(1:(length(pred_update$xi.map)-1), function(x) beta[[x]] - beta_update[[x]])
      Psi <- GramMatrixPhiAdditive(model_update$ulist)
      Psi_vec <- GramVectorPhiAdditive(model_update$ulist)
      
      zPsivec <- sapply(1:(length(pred_update$xi.map)-1), function(x) t(z[[x]]) %*% Psi_vec[[x]])
      zPsivec2 <- sapply(1:(length(pred_update$xi.map)-1), function(x) zPsivec[[x]]^2)
      zPsiz <- sapply(1:(length(pred_update$xi.map)-1), function(x) t(z[[x]]) %*% Psi[[x]] %*% z[[x]])
      sumZPsivec <- sum(zPsivec)
      znewdiff <- diff(unlist(tail(beta_update, n=1)))
      znewplus <- sum(unlist(tail(beta_update, n=1)))
      
      f <- sum(zPsiz) - sum(zPsivec2) + znewdiff^2/12 + (sum(zPsivec) - znewplus/2)^2
      f <- max(0, f)
      # # f2 <- 0
      # # for (k in 1:length(z))
      # #   f2 <- f2 + t(z[[k]]) %*% Psi[[k]] %*% z[[k]] - zPsivec[[k]]^2
      # # f2 <- f2 + znewdiff^2/12 + sumZPsivec^2 - znewplus*sumZPsivec + znewplus^2/4
      # # 
      
      # # if (MCtest) {
      #   f_MC <- mean((pred$PhiAll.test %*% unlist(pred$xi.map) - pred_update$PhiAll.test %*% unlist(beta_update))^2)
      #   message("Dimension addition - Analytic value: ", f, " Monte-Carlo: ", f_MC)
      # # }
    }
    f <- - f - reward_new_dim
    g <- NaN
  } else {
    model_update <-  create(class = 'lineqAGP',
                            x = model$x[, activeDim_IdxSeq], y = model$y,
                            constrType = model$constrType[activeDim_IdxSeq])
    model_update$varnoise <- model$varnoise
    model_update$nugget <- model$nugget
    model_update$kernParam <- model$kernParam[activeDim_IdxSeq]
    model_update$ulist <- ulist
    model_update$localParam$m <- sapply(ulist, length)
    
    optim_knot_pos <- optimize(f = distAdditiveModes, interval = c(0, 1),
                               model_update, pred, xtest[, activeDim_IdxSeq],
                               idx_add = which(activeDim_IdxSeq == idx_add),
                               reward_new_knot#, MCtest
                               )
    f <- optim_knot_pos$objective
    g <- optim_knot_pos$minimum
  }
  return(c(f, g))
}


# distModesNum <- function(x, y, modelNew, modelOld, xtest) {
#   PhiNew <- basisCompute.lineqAGP(c(x, y), modelNew$ulist, d = 2)
#   PhiOld <- basisCompute.lineqAGP(c(x, y), modelOld$ulist, d = 2)
#   predNew <- predict(modelNew, xtest = xtest)
#   predOld <- predict(modelOld, xtest = xtest)
#   return((PhiNew %*% predNew$xi.map - PhiOld %*% predOld$xi.map)^2)
# }
# 
# distModesNum <- function(u, model, xtest, idx_add) {
#   model_update <- model
#   model_update$ulist[[idx_add]] <- sort(c(model_update$ulist[[idx_add]], u))
#   f <- - integral2(modeDist,
#                    xmin = 0, xmax = 1, ymin = 0, ymax = 1,
#                    vectorized = FALSE, reltol = 1e-1, maxlist = 1e2,
#                    modelNew = model_update, modelOld = model, xtest = xtest)$Q
#   return(f)
# }


#' @title Gram matrix of the basis functions (\code{"lineqAGP"})
#' @description Compute the Gram matrix of the basis functions for \code{"lineqAGP"} models
#' 
#' @param u a sequence corresponding to the values of the knots
#' 
#' @return Gram matrix of the basis functions
#'
#' @author A. F. Lopez-Lopera
#'
#' @references F. Bachoc, A. F. Lopez-Lopera, and O. Roustant (2020),
#' "Sequential construction and dimension reduction of Gaussian processes under inequality constraints".
#' \emph{ArXiv e-prints}
#' <arXiv:2009.04188>
#'
#' @export
GramMatrixPhiAdditive <- function(u) {
  d <- length(u)
  Psi_list <- vector("list", d)
  for (kk in 1:d) {
    diffUVec <- diff(u[[kk]])
    
    diagPsi <- c(diffUVec[1], diffUVec[-1] + rev(rev(diffUVec)[-1]), rev(diffUVec)[1])/3
    Psi_temp <- diag(diagPsi)
    Psi_temp[abs(row(Psi_temp) - col(Psi_temp)) == 1] <- diffUVec/6
    Psi_list[[kk]] <- Psi_temp
  }
  return(Psi_list)
}

GramVectorPhiAdditive <- function(u) {
  d <- length(u)
  Psi_vec_list <- vector("list", d)
  for (kk in 1:d) {
    diffUVec <- diff(u[[kk]])
    
    Psi_vec <- c(diffUVec[1], diffUVec[-1] + rev(rev(diffUVec)[-1]), rev(diffUVec)[1])/2
    Psi_vec_list[[kk]] <- Psi_vec
  }
  return(Psi_vec_list)
}

#' @title Modification of the MAP estimate (\code{"lineqAGP"})
#' @description Compute the modification of the MAP estimate according to the MaxMod criterion
#' 
#' @param u an extended sequence corresponding to the values of the knots
#' @param model an object with class \code{lineqAGP}
#' @param pred an predictive model with class \code{lineqAGP}
#' @param xtest test input data
#' @param idx_add index of the dimension that will be activated
#' @param reward_new_knot a number corresponding to the reward of adding a new knot in an existing dimension
# #' @param MCtest a logical variable to approximate MaxMod criterion via MC
#' 
#' @return the modification of the MAP estimate after adding a new knot
#'
#' @author A. F. Lopez-Lopera
#'
#' @references F. Bachoc, A. F. Lopez-Lopera, and O. Roustant (2020),
#' "Sequential construction and dimension reduction of Gaussian processes under inequality constraints".
#' \emph{ArXiv e-prints}
#' <arXiv:2009.04188>
#'
#' @export
distAdditiveModes <- function(u, 
                              model,
                              pred, xtest,
                              idx_add,
                              reward_new_knot = 1e-6#, MCtest
                              ) {
  # pred <- predict(model, xtest = xtest)
  
  u_update <- model$ulist[[idx_add]] <- sort(c(model$ulist[[idx_add]], u))
  m <- sapply(model$ulist, length)
  
  model <- lineqGPOptim(model,
                        additive = TRUE,
                        estim.varnoise = TRUE, # to add this info at the MaxMod level
                        bounds.varnoise = c(1e-7, Inf), # to add this info at the MaxMod level
                        lb = rep(1e-2, model$d+1), ub = c(Inf, rep(1, model$d))) # to add this info at the MaxMod level
  pred_update <- predict(model, xtest = xtest)
  
  idx <- rev(which(u_update == u))
  betaIdxDimStr <- paste("beta_mod[", idx, "]", sep = "")
  betaIdxDimStrNeg <- paste("beta_mod[", -idx, "]", sep = "")
  betaIdxDimStrMinus <- paste("beta_mod[", idx-1, "]", sep = "")
  betaIdxDimStrPlus <- paste("beta_mod[", idx+1, "]", sep = "")
  
  beta_mod <- rep(0, length(pred_update$xi.map[[idx_add]]))
  eval(parse(text = paste(betaIdxDimStrNeg, "<- pred$xi.map[[idx_add]]")))
  
  diffU <- u_update[idx+1] - u_update[idx-1]
  hatCte1 <- (u_update[idx+1]-u_update[idx])/diffU
  hatCte2 <- (u_update[idx]-u_update[idx-1])/diffU
  eval(parse(text = paste(betaIdxDimStr, "<-",
                          hatCte1,"*",betaIdxDimStrMinus, "+",
                          hatCte2,"*",betaIdxDimStrPlus)))
  beta <- pred$xi.map
  beta[[idx_add]] <- beta_mod
  beta_update <- pred_update$xi.map
  
  z <- lapply(1:(length(pred_update$xi.map)), function(x) beta[[x]] - beta_update[[x]])
  Psi <- GramMatrixPhiAdditive(model$ulist)
  Psi_vec <- GramVectorPhiAdditive(model$ulist)
  zPsivec <- sapply(1:(length(pred_update$xi.map)), function(x) t(z[[x]]) %*% Psi_vec[[x]])
  zPsivec2 <- sapply(1:(length(pred_update$xi.map)), function(x) zPsivec[[x]]^2)
  zPsiz <- sapply(1:(length(pred_update$xi.map)), function(x) t(z[[x]]) %*% Psi[[x]] %*% z[[x]])
  
  f <- sum(zPsiz) - sum(zPsivec2) + sum(zPsivec)^2
  f <- max(0, f)
  
  # # if (MCtest) {
  #   f_MC <- mean((pred$PhiAll.test %*% unlist(pred$xi.map) - pred_update$PhiAll.test %*% unlist(beta_update))^2)
  #   message("Knot insertion - Analytic value: ", f, " Monte-Carlo: ", f_MC)
  # # }
  
  dist_nn <- min(abs(u_update[idx] - u_update[idx-1]), 
                 abs(u_update[idx+1] - u_update[idx]))
  f <- - f - reward_new_knot*dist_nn
  
  return(f)
}
