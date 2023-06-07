#' @title Creation Method for the \code{"lineqBAGP"} S3 Class
#' @description Creation method for the \code{"lineqBAGP"} S3 class.
#' 
#' @param partition a list of list making a partion of set $\{1,cdots,D\}$
#' @param subdivision a list of list containing sequence of subdivision of segment $[0,1]$
#' @param x a vector or matrix with the input data. The dimensions should be indexed by columns
#' @param y a vector with the output data
#' @param constrType a character string corresponding to the type of the inequality constraint
#' @param m a scalar or vector corresponding to the number of knots per dimension.
#' Options: "boundedness", "monotonicity", "convexity", "linear"
#' Multiple constraints can be also defined, e.g. \code{constrType = c("boundedness", "monotonicity")}
#' 
#' @return A list with the following elements.
#' \item{x,y,constrType}{see \bold{Arguments}}
#' \item{d}{a number corresponding to the input dimension}
#' \item{constrIdx}{for d > 1, a integer vector with the indices of active constrained dimensions}
#' \item{constrParam}{constraint inequalities for each dimension}
#' \item{varnoise}{a scalar with noise variance}
#' \item{localParam}{a list with specific parameters required for \code{"lineqBAGP"} models:
#' \code{m} (number of basis functions), \code{sampler}, and \code{samplingParams}.
#' See \code{\link{simulate.lineqBAGP}}}
#' \item{kernParam}{a list with the kernel parameters: \code{par} (kernel parameters), \code{type}, \code{nugget}.
#' See \code{\link{kernCompute}}}
#' \item{bounds}{the limit values if \code{constrType = "boundedness"}.}
#' \item{(Lambda,lb,ub)}{the linear system of inequalities if \code{constrType = "linear"}}
#'
#' @seealso \code{\link{augment.lineqBAGP}}, \code{\link{predict.lineqBAGP}}, \code{\link{simulate.lineqBAGP}}
#' 
#' @author A. F. Lopez-Lopera
#' 
#' @references A. F. Lopez-Lopera (2019),
#' "Gaussian process modelling under inequality constraints".
#' \emph{PhD thesis, Mines Saint-Etienne}
#' <https://tel.archives-ouvertes.fr/tel-02863891>
#'
#' @examples
#' # creating the model
#' d <- 2
#' fun1 <- function(x) return(4*(x-0.5)^2)
#' fun2 <- function(x) return(2*x)
#' targetFun <- function(x) return(fun1(x[, 1]) + fun1(x[, 2])) 
#' xgrid <- expand.grid(seq(0, 1, 0.01), seq(0, 1, 0.01))
#' ygrid <- targetFun(xgrid)
#' xdesign <- rbind(c(0.5, 0), c(0.5, 0.5), c(0.5, 1), c(0, 0.5), c(1, 0.5))
#' ydesign <- targetFun(xdesign)
#' model <- create(class = "lineqBAGP", x = xdesign, y = ydesign,
#'                 constrType = c("convexity", "monotonicity"))
#' str(model)
#' 
#' @method create lineqBAGP
#' @export
 
create.lineqBAGP <- function(x, y, constrType,
                             partition = as.list(seq(ncol(x))),
                             m = NULL) {
  # changing the data as matrices
  if (!is.matrix(x))
    x <- as.matrix(x)
  if (!is.matrix(y) || ncol(y) != 1)
    y <- matrix(y)
  
  d <- ncol(x) # dimension of the input space
  nblock <- length(partition) # nb of blocks
  dblock <- sapply(partition, function(x) length(x))
  
  #biject is a function that take an element 'i' and gives the couple [j,k] 
  #corresponding such that partition[j,k]=i
  bijection <- function(partition, i){
    for (j in 1:length(partition)){
      if (i %in% partition[[j]]){
        return (c(j,which(partition[[j]]==i)))
      }
    }
  }
  
  #Bij <-lapply(c(1:3), function(x) biject(partition,x)) useful?
  
  #Creation mlist if m.type=list[seq]
  mlist <- lapply(1:nblock, function(x) rep(10, dblock[x]))
  if (length(m)==d){
    for (k in 1:length(m)){
      pos<-bijection(partition, k)
      mlist[[pos[1]]][pos[2]] <- m[k]
    }
  }
  #Creation mlist if m.type=float
  else if (length(m) == 1) {
    mlist <- lapply(1:nblock, function(x) rep(m, dblock[x]))
  }
  # else {
  #   stop("Can not deal with this type of 'm', try: m=NULL, m=c(a_1,....,a_d),
  #        m=a")
  # } 
  # "To see later "
  #### to verify if the partition is a disjoint one and to check that the union is a sequence 1:D ###
  
  subdivision <- vector("list", nblock)
  for (j in 1:nblock) {
    subdivision[[j]] <- vector("list", dblock[j])
    for (k in 1:dblock[j])
      subdivision[[j]][[k]] <- matrix(seq(0, 1, by = 1/(mlist[[j]][k]-1)), ncol = 1) # discretization vector
  }
  
  
  # creating some lists for the model
  localParam <- list(mlist = mlist, 
                     sampler = "ExpT",
                     nblock = nblock, partition = partition, dblock = dblock, 
                     samplingParams = c(thinning = 1, burn.in = 1, scale = 0.1))
  
  
  names(localParam$partition) <- paste('block', 1:localParam$nblock, sep = "")
  
  ##Constraints 
  
  constrFlags <- rep(1, nblock) # to be checked later 
  
  kernParam <- vector("list", localParam$nblock)
  constrParam <- vector("list", localParam$nblock)
  for (k in 1:localParam$nblock) { # to be checked later! We will focus on the monotonicity constraint
    kernParam[[k]] <- list(par = c(sigma2 = 1^2, theta = rep(0.1, dblock[k])), type = "matern52")#, nugget = 1e-7*sd(y))
    switch (constrType[k],
            boundedness = {
              constrParam[[k]]$bounds <- c(lower = min(y) - 0.05*abs(max(y) - min(y)),
                                           upper = max(y) + 0.05*abs(max(y) - max(y)))
            }, monotonicity = {
              constrParam[[k]]$bounds <- c(0, Inf)
            }, decreasing = {
              constrParam[[k]]$bounds <- c(0, Inf)
            }, convexity = {
              constrParam[[k]]$bounds <- c(0, Inf)
            }, linear = {
              # imposing positiveness constraints by default
              constrParam[[k]]$Lambda <- diag(prod(model$localParam$m))
              constrParam[[k]]$lb <- rep(0, nrow(model$Lambda))
              constrParam[[k]]$ub <- rep(Inf, nrow(model$Lambda))
            }, none = {
              constrParam[[k]]$bounds <- c(-Inf, Inf)
              constrFlags[k] <- 0
            }
    )
  }
  constrIdx <- which(constrFlags == 1)
  # creating the full list for the model
  model <- list(x = x, y = y, constrType = constrType,  subdivision = subdivision,
                d = d,  nugget = 1e-9, 
                constrIdx = constrIdx, constrParam = constrParam,
                varnoise = 0,  localParam = localParam, kernParam = kernParam)
  return(model)
}

#' @title Augmenting Method for the \code{"lineqBAGP"} S3 Class
#' @description Augmenting method for the \code{"lineqBAGP"} S3 class.
#' 
#' @param x an object with class \code{lineqGP}
#' @param ... further arguments passed to or from other methods
#' 
#' @return An expanded \code{"lineqGP"} object with the following additional elements
#' \item{Phi}{a matrix corresponding to the hat basis functions.
#' The basis functions are indexed by rows}
#' \item{Gamma}{the covariance matrix of the Gassian vector \eqn{\boldsymbol{\xi}}{\xi}.}
#' \item{(Lambda,lb,ub)}{the linear system of inequalities.}
#' \item{...}{further parameters passed to or from other methods.}
#'
#' @details Some paramaters of the finite-dimensional GP with linear inequality
#' constraints are computed. Here, \eqn{\boldsymbol{\xi}}{\xi} is a centred Gaussian
#' vector with covariance \eqn{\boldsymbol{\Gamma}}{\Gamma}, s.t.
#' \eqn{\boldsymbol{\Phi} \boldsymbol{\xi} = \boldsymbol{y}}{\Phi \xi = y}
#' (interpolation constraints) and
#' \eqn{\boldsymbol{l} \leq \boldsymbol{\Lambda} \boldsymbol{\xi} \leq \boldsymbol{u}}{lb \le \Lambda \xi \le ub}
#' (inequality constraints).
#' 
#' @seealso \code{\link{create.lineqBAGP}}, \code{\link{predict.lineqBAGP}},
#'          \code{\link{simulate.lineqBAGP}}
#'          
#' @author A. F. Lopez-Lopera
#'
#' @references A. F. Lopez-Lopera (2019),
#' "Gaussian process modelling under inequality constraints".
#' \emph{PhD thesis, Mines Saint-Etienne}
#' <https://tel.archives-ouvertes.fr/tel-02863891>
#' 
#' @examples
#' # creating the model
#' d <- 2
#' fun1 <- function(x) return(4*(x-0.5)^2)
#' fun2 <- function(x) return(2*x)
#' targetFun <- function(x) return(fun1(x[, 1]) + fun1(x[, 2])) 
#' xgrid <- expand.grid(seq(0, 1, 0.01), seq(0, 1, 0.01))
#' ygrid <- targetFun(xgrid)
#' xdesign <- rbind(c(0.5, 0), c(0.5, 0.5), c(0.5, 1), c(0, 0.5), c(1, 0.5))
#' ydesign <- targetFun(xdesign)
#' model <- create(class = "lineqBAGP", x = xdesign, y = ydesign,
#'                 constrType = c("convexity", "monotonicity"), m = 50)
#' 
#' # updating and expanding the model
#' model$kernParam[[1]]$par <- c(1, 0.2)
#' model$kernParam[[2]]$par <- c(1, 0.2)
#' model$nugget <- 1e-9
#' model$varnoise <- 1e-5
#' model <- augment(model)
#' str(model)
#'
#' @importFrom broom augment
#' @export

augment.lineqBAGP<- function(x, ...) {
  model <- x
  if (!("nugget" %in% names(model)))
    model$nugget <- 0
  if ("bounds" %in% names(model)) {
    bounds <- model$bounds
  } else {
    bounds <- c(0, Inf)
  }
  # passing some terms from the model
  x <- model$x
  nblock <- model$localParam$nblock
  dblock <- model$localParam$dblock
  
  # computing the kernel matrix for the prior
  subdivision <- model$subdivision
  mlist <- model$localParam$mlist
  mtotal <- model$localParam$mtotal <- sum(sapply(mlist, prod))
  
  Gamma <- Phi <- vector("list", nblock)
  for (j in 1:nblock) {
    Gamma[[j]] <- Phi[[j]] <- vector("list", dblock[j])
    for (k in 1:dblock[j]) {
      Gamma[[j]][[k]] <- kernCompute(subdivision[[j]][[k]], subdivision[[j]][[k]], model$kernParam[[j]]$type,
                                     model$kernParam[[j]]$par[c(1, k+1)])
      
      Phi[[j]][[k]] <- basisCompute.lineqGP(x[, model$localParam$partition[[j]][k]],
                                            subdivision[[j]][[k]])
    }
  }
  
  model$subdivision <- subdivision
  model$Gamma <- Gamma
  model$Phi <- Phi
  
  # # precomputing the linear system for the QP solver and MCMC samplers
  # M <- g <- vector("list", nblock)
  # mvec <- vector("list", nblock)
  # for (k in 1:nblock) {
  #   mvec[[k]] <- c()
  #   if (model$constrType[k] == "linear") {
  #     if (!("Lambda" %in% names(model)))
  #       stop('matrix Lambda is not defined')
  #     Lambda <- model$constrParam[[k]]$Lambda
  #     lb <- model$constrParam[[k]]$lb
  #     ub <- model$constrParam[[k]]$ub
  #     lsys <- lineqGPSys(nrow(Lambda), model$constrType[k], lb, ub,
  #                        Lambda, lineqSysType = "oneside")
  #     lsys2 <- lineqGPSys(nrow(Lambda), model$constrType[k], lb, ub,
  #                         Lambda, rmInf = FALSE)
  #   } else {
  #     bounds <- model$constrParam[[k]]$bounds
  #     lsys <- lineqGPSys(m[k], model$constrType[k], bounds[1], bounds[2],
  #                        constrIdx = model$constrIdx, lineqSysType = "oneside")
  #     lsys2 <- lineqGPSys(m[k], model$constrType[k], bounds[1], bounds[2],
  #                         constrIdx = model$constrIdx, rmInf = FALSE)
  #   }
  #   # oneside linear structure for QP.solver: M = [Lambda,-Lambda] and g = [-lb,ub]
  #   M[[k]] <- lsys$M
  #   g[[k]] <- -matrix(lsys$g)
  #   # twosides linear structure (Lambda, lb, ub) for MCMC samplers
  #   model$constrParam[[k]]$Lambda <- lsys2$A
  #   model$constrParam[[k]]$lb <- lsys2$l
  #   model$constrParam[[k]]$ub <- lsys2$u
  #   # extra term required for HMC sampler
  #   mvec[[k]] <- nrow(lsys2$A)
  # }
  # # adding the parameters to the model structure
  # model$lineqSys$M <- M # for QP solve
  # model$lineqSys$g <- g # for QP solve
  # model$localParam$mvec <- mvec # for HMC sampler
  # 
  # model$lineqSys$MAll <- eval(parse(text = paste("bdiag(",
  #                                                paste("model$lineqSys$M[[", 1:nblock, "]]",
  #                                                      sep = "", collapse = ","),
  #                                                ")", sep = "")))
  # model$lineqSys$MAll <- matrix(model$lineqSys$MAll, ncol = mtotal)
  # model$lineqSys$gAll <- matrix(unlist(model$lineqSys$g), ncol = 1)
  # model$lineqSys$LambdaAll <- eval(parse(text = paste("bdiag(",
  #                                                     paste("model$constrParam[[", 1:nblock, "]]$Lambda",
  #                                                           sep = "", collapse = ","),
  #                                                     ")", sep = "")))
  # model$lineqSys$LambdaAll <- matrix(model$lineqSys$LambdaAll, ncol = mtotal)  
  # model$lineqSys$lb <- eval(parse(text = paste("c(",
  #                                              paste("model$constrParam[[", 1:nblock, "]]$lb",
  #                                                    sep = "", collapse = ","),
  #                                              ")", sep = "")))
  # model$lineqSys$ub <- eval(parse(text = paste("c(",
  #                                              paste("model$constrParam[[", 1:nblock, "]]$ub",
  #                                                    sep = "", collapse = ","),
  #                                              ")", sep = "")))
  
  
  return(model)
}

#' @title Prediction Method for the \code{"lineqBAGP"} S3 Class
#' @description Prediction method for the \code{"lineqBAGP"} S3 class.
#' 
#' @param object an object with class \code{"lineqBAGP"}.
#' @param xtest a vector (or matrix) with the test input design
#' @param return_model If \code{TRUE}, the augmented model is returned (see \code{\link{augment.lineqBAGP}}).
#' @param ... further arguments passed to or from other methods
#' 
#' @return A \code{"lineqBAGP"} object with the following elements.
#' \item{Lambda}{a matrix corresponding to the linear set of inequality constraints}
#' \item{lb}{the lower bound vector of the inequalities constraints}
#' \item{ub}{the upper bound vector of the inequalities constraints}
#' \item{Phi.test}{a matrix corresponding to the hat basis functions evaluated
#' at \code{xtest}. The basis functions are indexed by rows}
#' \item{mu}{the unconstrained GP mean predictor}
#' \item{Sigma}{the unconstrained GP prediction conditional covariance matrix}
#' \item{xi.map}{the GP maximum a posteriori (MAP) predictor given the inequality constraints}
#'
#' @details The posterior paramaters of the finite-dimensional GP with linear inequality
#' constraints are computed. Here, \eqn{\boldsymbol{\xi}}{\xi} is a centred Gaussian
#' vector with covariance \eqn{\boldsymbol{\Gamma}}{\Gamma}, s.t.
#' \eqn{\boldsymbol{\Phi} \boldsymbol{\xi} = \boldsymbol{y}}{\Phi \xi = y}
#' (interpolation constraints) and
#' \eqn{\boldsymbol{l} \leq \boldsymbol{\Lambda} \boldsymbol{\xi} \leq \boldsymbol{u}}{lb \le \Lambda \xi \le ub}
#' (inequality constraints).
#' 
#' @seealso \code{\link{create.lineqBAGP}}, \code{\link{augment.lineqBAGP}},
#'          \code{\link{simulate.lineqBAGP}}
#'          
#' @author A. F. Lopez-Lopera
#'
#' @references A. F. Lopez-Lopera (2019),
#' "Gaussian process modelling under inequality constraints".
#' \emph{PhD thesis, Mines Saint-Etienne}
#' <https://tel.archives-ouvertes.fr/tel-02863891>
#'
#' @examples
#' library(plot3D)
#' # creating the model
#' d <- 2
#' fun1 <- function(x) return(4*(x-0.5)^2)
#' fun2 <- function(x) return(2*x)
#' targetFun <- function(x) return(fun1(x[, 1]) + fun2(x[, 2])) 
#' xgrid <- expand.grid(seq(0, 1, 0.01), seq(0, 1, 0.01))
#' ygrid <- targetFun(xgrid)
#' xdesign <- rbind(c(0.5, 0), c(0.5, 0.5), c(0.5, 1), c(0, 0.5), c(1, 0.5))
#' ydesign <- targetFun(xdesign)
#' model <- create(class = "lineqBAGP", x = xdesign, y = ydesign,
#'                 constrType = c("convexity", "monotonicity"), 10)
#' 
#' # updating and expanding the model
#' model$kernParam[[1]]$type <- "matern52"
#' model$kernParam[[2]]$type <- "matern52"
#' model$kernParam[[1]]$par <- c(1, 0.2)
#' model$kernParam[[2]]$par <- c(1, 0.3)
#' model$nugget <- 1e-9
#' model$varnoise <- 1e-5
#'
#' # predictions from the model
#' ntest <- 25
#' xtest  <- cbind(seq(0, 1, length = ntest), seq(0, 1, length = ntest))
#' ytest <- targetFun(xtest)
#' pred <- predict(model, xtest)
#' persp3D(x = unique(xtest[, 1]), y = unique(xtest[, 2]),
#'         z = outer(c(pred$Phi.test[[1]] %*% pred$xi.map[[1]]),
#'                   c(pred$Phi.test[[2]] %*% pred$xi.map[[2]]), "+"),
#'         xlab = "x1", ylab = "x2", zlab = "mode(x1,x2)", zlim = c(0, 3),
#'         phi = 20, theta = -30, alpha = 1, colkey = FALSE)
#' points3D(x = xdesign[,1], y = xdesign[,2], z = ydesign, col = "black", pch = 19, add = TRUE)
#'
#' @importFrom quadprog solve.QP
#' @importFrom Matrix bdiag
#' @import plot3D
#' @export
predict.lineqBAGP <- function(object, xtest, return_model = FALSE, ...) {
  model <- augment(object)
  if (!is.matrix(xtest))
    xtest <- matrix(xtest, ncol = model$d)
  
  nblock <- model$localParam$nblock
  dblock <- model$localParam$dblock
  
  nt <- length(model$y) 
  mt <- model$localParam$mtotal
  
  # passing some terms from the model
  pred <- list()
  class(pred) <- class(model)
  pred$constrParam <- model$constrParam
  
  # precomputing some terms
  
  Phi.test <- vector("list", nblock)
  for (j in 1:nblock) {
    Phi.test[[j]] <- vector("list", dblock[j])
    for (k in 1:dblock[j]) {
      Phi.test[[j]][[k]] <- basisCompute.lineqGP(xtest[, model$localParam$partition[[j]][k]],
                                                 model$subdivision[[j]][[k]])
    }
  }
  pred$Phi.test <- Phi.test
  pred$PhiAll.test <- list("vector", nblock)
  for (j in 1:nblock) {
    pred$PhiAll.test[[j]] <- matrix(0, nt, prod(model$localParam$mlist[[j]]))
    for (iterObs in 1:nt) {
      hfun <- eval(parse(text = paste("Phi.test[[j]][[", 1:dblock[j], "]][", iterObs, ", ]",
                                      sep = "", collapse = "%x%")))
      pred$PhiAll.test[[j]][iterObs, ] <- eval(hfun)
    }
  }
  
  # # computing the conditional mean vector and conditional covariance matrix
  # # given the interpolation points
  GammaBlocks <- phiBlocks <- list("vector", nblock)
  for (j in 1:nblock){
    phiBlocks[[j]] <- matrix(0, nt, prod(model$localParam$mlist[[j]]))
    for (iterObs in 1:nt) {
      hfun <- eval(parse(text = paste("model$Phi[[j]][[", 1:dblock[j], "]][", iterObs, ", ]",
                                      sep = "", collapse = "%x%")))
      phiBlocks[[j]][iterObs, ] <- eval(hfun)
    }
    
    hfun <- eval(parse(text = paste("model$Gamma[[j]][[", 1:dblock[j], "]]",
                                    sep = "", collapse = "%x%")))
    GammaBlocks[[j]] <- eval(hfun)
  }
  
  hfunBigPhi <- parse(text = paste("cbind(",
                                   paste("phiBlocks[[", 1:nblock, "]]", sep = "", collapse = ","),
                                   ")", sep = ""))
  bigPhi <- eval(hfunBigPhi)
  
  hfunBigGamma <- parse(text = paste("bdiag(",
                                     paste("GammaBlocks[[", 1:nblock, "]]", sep = "", collapse = ","),
                                     ")", sep = ""))
  bigGamma <- eval(hfunBigGamma)
  
  ##Function to have stable inverse of symetric matrices
  stable_inv <-function(Sym){
    cholSym <- lapply(Sym, function(x) chol(x))
    invSym <- lapply(cholSym, function(x) chol2inv(x))
  }
  cholGammaBlocks <- lapply(GammaBlocks, function(x) chol(x))
  invGammaBlocks <- lapply(cholGammaBlocks, function(x) chol2inv(x))
  hfunInvBigGamma <- parse(text = paste("bdiag(",
                                        paste("invGammaBlocks[[", 1:nblock, "]]", sep = "", collapse = ","),
                                        ")", sep = ""))
  invBigGamma <- eval(hfunInvBigGamma)
  invSigmaAll <- as.matrix(invBigGamma) + t(bigPhi)%*%bigPhi/model$varnoise
  
  #τ^{-2}[In − ΨL(τ^2 I_m + (ΨL)^t ΨL)^{−1}(ΨL)^t ]
  ### to be adapted !! ###
  if (mt < nt) {
    hfunBigCholGamma <- parse(text = paste("bdiag(",
                                           paste("cholGammaBlocks[[", 1:nblock, "]]", sep = "", collapse = ","),
                                           ")", sep = ""))
    bigCholGamma <- as.matrix(eval(hfunBigCholGamma))
    PhibigCholGamma <- bigPhi %*% bigCholGamma
    ILtPhitPhiL <- model$varnoise*diag(mt) + t(PhibigCholGamma) %*% PhibigCholGamma
    # ILtPhitPhiL <- ILtPhitPhiL + model$nugget*nrow(ILtPhitPhiL)
    cholILtPhitPhiL <- t(chol(ILtPhitPhiL))
    Lschur <- forwardsolve(cholILtPhitPhiL, t(PhibigCholGamma))
    invPhiGammaPhitFull <- (diag(nt) - t(Lschur)%*% Lschur)/model$varnoise
  } else {
    PhiGammaPhit <- as.matrix(bigPhi %*% bigGamma %*% t(bigPhi))
    PhiGammaPhit <- PhiGammaPhit + model$nugget*diag(nrow(PhiGammaPhit))
    invPhiGammaPhitFull <- chol2inv(chol(PhiGammaPhit + model$varnoise * diag(nt))) # instability issues here
  }
  

  temp <- bigGamma%*%t(bigPhi) %*% invPhiGammaPhitFull
  muAll <- as.vector(temp %*% model$y)
  
  plot(ydesign, bigPhi%*%muAll)
  
  pred$mu <- mu
  pred$Sigma <- Sigma
  pred$muAll <- muAll #matrix(c(mu), ncol = 1)
  pred$SigmaAll <- SigmaAll
  
  
  xiAll.map <- solve.QP(invSigmaAll, t(pred$muAll) %*% invSigmaAll,
                        t(model$lineqSys$MAll), model$lineqSys$gAll)$solution
  ### to be adapted !! ###
  
  
  pred$xiAll.map <- xiAll.map
  
  pred$xi.map <- vector("list", nblock)
  for (k in 1:nblock) {
    pred$xi.map[[k]] <- xiAll.map[(mlist[k]+1):mlist[k+1]]
  }
  
  if (return_model)
    pred$model <- model
  return(pred)
}

#' @title Simulation Method for the \code{"lineqBAGP"} S3 Class
#' @description Simulation method for the \code{"lineqBAGP"} S3 class.
#' 
#' @param object an object with class \code{"lineqBAGP"}
#' @param nsim	the number of simulations
#' @param seed see \code{\link{simulate}}
#' @param xtest a vector (or matrix) with the test input design
#' @param ... further arguments passed to or from other methods
#' 
#' @return A \code{"lineqBAGP"} object with the following elements
#' \item{x}{a vector (or matrix) with the training input design}
#' \item{y}{the training output vector at \code{x}}
#' \item{xtest}{a vector (or matrix) with the test input design}
#' \item{Phi.test}{a matrix corresponding to the hat basis functions evaluated
#' at \code{xtest}. The basis functions are indexed by rows.}
#' \item{xi.sim}{the posterior sample-path of the finite-dimensional Gaussian vector}
#' \item{ysim}{the posterior sample-path of the observed GP
#' Note: \code{ysim = Phi.test \%*\% xi.sim}}
#'
#' @details The posterior sample-path of the finite-dimensional GP with linear inequality
#' constraints are computed. Here, \eqn{\boldsymbol{\xi}}{\xi} is a centred Gaussian
#' vector with covariance \eqn{\boldsymbol{\Gamma}}{\Gamma}, s.t.
#' \eqn{\boldsymbol{\Phi} \boldsymbol{\xi} = \boldsymbol{y}}{\Phi \xi = y}
#' (interpolation constraints) and
#' \eqn{\boldsymbol{l} \leq \boldsymbol{\Lambda} \boldsymbol{\xi} \leq \boldsymbol{u}}{lb \le \Lambda \xi \le ub}
#' (inequality constraints).
#' 
#' @seealso \code{\link{create.lineqBAGP}}, \code{\link{augment.lineqBAGP}},
#'          \code{\link{predict.lineqBAGP}}
#'          
#' @author A. F. Lopez-Lopera
#'
#' @references A. F. Lopez-Lopera (2019),
#' "Gaussian process modelling under inequality constraints".
#' \emph{PhD thesis, Mines Saint-Etienne}
#' <https://tel.archives-ouvertes.fr/tel-02863891>
#'
#' @examples
#' library(plot3D)
#' # creating the model
#' d <- 2
#' fun1 <- function(x) return(4*(x-0.5)^2)
#' fun2 <- function(x) return(2*x)
#' targetFun <- function(x) return(fun1(x[, 1]) + fun2(x[, 2])) 
#' xgrid <- expand.grid(seq(0, 1, 0.01), seq(0, 1, 0.01))
#' ygrid <- targetFun(xgrid)
#' xdesign <- rbind(c(0.5, 0), c(0.5, 0.5), c(0.5, 1), c(0, 0.5), c(1, 0.5))
#' ydesign <- targetFun(xdesign)
#' model <- create(class = "lineqBAGP", x = xdesign, y = ydesign,
#'                 constrType = c("convexity", "monotonicity"), m = 10)
#' 
#' # updating and expanding the model
#' model$kernParam[[1]]$type <- "matern52"
#' model$kernParam[[2]]$type <- "matern52"
#' model$kernParam[[1]]$par <- c(1, 0.2)
#' model$kernParam[[2]]$par <- c(1, 0.3)
#' model$nugget <- 1e-9
#' model$varnoise <- 1e-5
#'
#' # sampling from the model
#' ntest <- 25
#' xtest  <- cbind(seq(0, 1, length = ntest), seq(0, 1, length = ntest))
#' ytest <- targetFun(xtest)
#' sim.model <- simulate(model, nsim = 1e3, seed = 1, xtest = xtest)
#' PhiAll.test <- cbind(sim.model$Phi.test[[1]][rep(1:ntest, times = ntest), ],
#'                      sim.model$Phi.test[[2]][rep(1:ntest, each = ntest), ])
#' persp3D(x = unique(xtest[, 1]), y = unique(xtest[, 2]),
#'         z = matrix(rowMeans(PhiAll.test %*% sim.model$xiAll.sim), ntest, ntest),
#'         xlab = "x1", ylab = "x2", zlab = "mode(x1,x2)", zlim = c(0, 3),
#'         phi = 20, theta = -30, alpha = 1, colkey = FALSE)
#' points3D(x = xdesign[,1], y = xdesign[,2], z = ydesign, col = "black", pch = 19, add = TRUE)
#'
#' @importFrom stats simulate
#' @import plot3D
#' @export
simulate.lineqBAGP <- function(object, nsim = 1, seed = NULL, xtest, ...) {
  model <- object
  predtime <- proc.time()
  pred <- predict(model, xtest, return_model = TRUE)
  predtime <- proc.time() - predtime
  model <- pred$model
  pred$model <- NULL
  
  nblock <- model$localParam$nblock
  
  # # computing the transformed conditional mode and covariance matrix given
  # # the interpolation points and the inequality constraints
  # ysim <- xi.sim <- vector("list", nblock)
  # simtime <- proc.time()
  # for (k in 1:nblock) {
  #   Lambda <- pred$constrParam[[k]]$Lambda
  #   eta.map <- as.vector(Lambda %*% pred$xi.map[[k]])
  #   Sigma.eta <- Lambda %*% pred$Sigma[[k]] %*% t(Lambda)
  #   if (min(eigen(Sigma.eta, symmetric = TRUE)$values) < 0)
  #     Sigma.eta <- Sigma.eta + model$nugget*diag(nrow(Sigma.eta))
  #   
  #   # listing control terms
  #   control <- as.list(unlist(model$localParam$samplingParam))
  #   control$mvec <- model$localParam$mvec[[k]] # for HMC
  #   control$constrType <- model$constrType[k] # for HMC
  #   
  #   # sampling from the truncated multinormal
  #   tmvPar <- list(mu = eta.map, Sigma = Sigma.eta,
  #                  lb = pred$constrParam[[k]]$lb,
  #                  ub = pred$constrParam[[k]]$ub)
  #   class(tmvPar) <- model$localParam$sampler
  #   set.seed(seed)
  #   eta <- tmvrnorm(tmvPar, nsim, control)
  #   xi.sim[[k]] <- qr.solve(Lambda, eta)
  #   ysim[[k]] <- pred$Phi.test[[k]] %*% xi.sim[[k]]
  # }
  # simtime <- proc.time() - simtime
  
  etaAll.map <- as.vector(model$lineqSys$LambdaAll %*% pred$xiAll.map)
  Sigma.etaAll <- model$lineqSys$LambdaAll %*% pred$SigmaAll %*% t(model$lineqSys$LambdaAll)
  Sigma.etaAll <- Sigma.etaAll + model$nugget*diag(nrow(Sigma.etaAll))
  
  # listing control terms
  control <- as.list(unlist(model$localParam$samplingParam))
  control$mvec <- model$localParam$mtotal # for HMC
  control$constrType <- model$constrType # for HMC
  
  tmvPar <- list(mu = etaAll.map, Sigma = Sigma.etaAll,
                 lb = model$lineqSys$lb,
                 ub = model$lineqSys$ub)
  class(tmvPar) <- model$localParam$sampler
  set.seed(seed)
  simtime <- proc.time()
  etaAll <- tmvrnorm(tmvPar, nsim, control)
  simtime <- proc.time() - simtime
  xiAll.sim <- qr.solve(model$lineqSys$LambdaAll, etaAll)
  
  # passing some terms to the simulated model
  simModel <- list()
  simModel$x <- model$x
  simModel$y <- model$y
  simModel$xtest <- xtest
  simModel$Phi.test <- pred$Phi.test
  # simModel$xi.map <- pred$xi.map
  # simModel$xi.sim <- xi.sim
  # simModel$ysim <- ysim
  simModel$predtime <- predtime
  simModel$simtime <- simtime
  
  simModel$muAll <- pred$muAll
  simModel$xiAll.map <- pred$xiAll.map
  simModel$xiAll.sim <- xiAll.sim
  simModel$PhiAll.test <- pred$PhiAll.test
  
  class(simModel) <- class(model)
  return(simModel)
}

