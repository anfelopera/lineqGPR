#' @title Creation Method for the \code{"lineqBAGP"} S3 Class
#' @description Creation method for the \code{"lineqBAGP"} S3 class.
#' 
#' @param x a vector or matrix with the input data. The dimensions should be indexed by columns
#' @param y a vector with the output data
#' @param constrType a character string corresponding to the type of the inequality constraint
#' @param m a scalar or vector corresponding to the number of knots per dimension.
#' @param partition a list of list making a partion of set $\{1,cdots,D\}$
#' Options: "boundedness", "monotonicity", "convexity", "linear"
#' Multiple constraints can be also defined, e.g. \code{constrType = c("boundedness", "monotonicity")}
#' 
#' @return A list with the following elements.
#' \item{x,y,constrType}{see \bold{Arguments}}
#' \item{d}{a number corresponding to the input dimension}
#' \item{constrIdx}{for d > 1, a integer vector with the indices of active constrained dimensions}
#' \item{constrParam}{constraint inequalities for each dimension}
#' \item{varnoise}{a scalar with noise variance}
#' \item{subdivision}{a list of list containing sequence of subdivision of segment $[0,1]$}
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
#' targetFun <- function(x) {
#'   return(x[, 1]*x[, 3] + x[, 2])
#' }
#' xdesign <- matrix(runif(12), 4, 3)
#' ydesign <- targetFun(xdesign)
#' d <- 3 # number of dimensions
#' nblocks <- 2
#' model <- create(class = "lineqBAGP", x = xdesign, y = ydesign,
#'                 constrType = rep("monotonicity", nblocks),
#'                 partition = list(c(1,3), 2),
#'                 m = c(5, 3, 4))
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
  
  # bijection is a function that take an element 'i' and gives the couple [j,k] 
  # corresponding such that partition[j,k]=i
  bijection <- function(partition, i) {
    for (j in 1:length(partition)) {
      if (i %in% partition[[j]]) {
        return (c(j, which(partition[[j]] == i)))
      }
    }
  }
  
  #Bij <-lapply(c(1:3), function(x) biject(partition,x)) useful?
  
  #Creation mlist if m.type=list[seq]
  mlist <- lapply(1:nblock, function(x) rep(10, dblock[x]))
  if (length(m) == d){
    for (k in 1:length(m)){
      pos <- bijection(partition, k)
      mlist[[pos[1]]][pos[2]] <- m[k]
    }
  } else if (length(m) == 1) {
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
    subdivision[[j]] <- lapply(1:dblock[j],
                               function(k) matrix(seq(0, 1, by = 1/(mlist[[j]][k]-1)), ncol = 1))
    names(subdivision[[j]]) <- names(mlist[[j]]) <- paste("x", partition[[j]], sep = "")
  }
    
  # creating some lists for the model
  names(mlist) <- names(subdivision) <- names(partition) <- paste("block", 1:nblock, sep = "")
  localParam <- list(mlist = mlist, 
                     sampler = "ExpT",
                     nblock = nblock, partition = partition, dblock = dblock, 
                     samplingParams = c(thinning = 1, burn.in = 1, scale = 0.1))
  
  ##Constraints 
  
  constrFlags <- rep(1, nblock) # to be checked later 
  
  kernParam <- constrParam <- vector("list", nblock)
  names(kernParam) <- names(constrParam) <- paste("block", 1:nblock, sep = "")
  for (j in 1:nblock) { # to be checked later! We will focus on the monotonicity constraint
    kernParam[[j]] <- list(par = c(sigma2 = 1^2, theta = rep(0.1, dblock[j])), type = "matern52")#, nugget = 1e-7*sd(y))
    switch (constrType[j],
            boundedness = {
              constrParam[[j]]$bounds <- c(lower = min(y) - 0.05*abs(max(y) - min(y)),
                                           upper = max(y) + 0.05*abs(max(y) - max(y)))
            }, monotonicity = {
              constrParam[[j]]$bounds <- c(0, Inf)
            }, decreasing = {
              constrParam[[j]]$bounds <- c(0, Inf)
            }, convexity = {
              constrParam[[j]]$bounds <- c(0, Inf)
            }, linear = { # to be checked!
              # imposing positiveness constraints by default
              constrParam[[j]]$Lambda <- diag(prod(model$localParam$m))
              constrParam[[j]]$lb <- rep(0, nrow(model$Lambda))
              constrParam[[j]]$ub <- rep(Inf, nrow(model$Lambda))
            }, none = {
              constrParam[[j]]$bounds <- c(-Inf, Inf)
              constrFlags[j] <- 0
            }
    )
  }
  constrIdx <- which(constrFlags == 1) # to be checked !!
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
  partition <- model$localParam$partition
  subdivision <- model$subdivision
  mlist <- model$localParam$mlist
  mtotal <- model$localParam$mtotal <- sum(sapply(mlist, prod))
  
  Gamma.perVar <- Phi.perVar <- vector("list", nblock)
  names(Gamma.perVar) <- names(Phi.perVar) <- paste("block", 1:nblock, sep = "")
  for (j in 1:nblock) {
    Gamma.perVar[[j]] <- lapply(1:dblock[j],
                                function(k) kernCompute(subdivision[[j]][[k]],
                                                        subdivision[[j]][[k]],
                                                        model$kernParam[[j]]$type,
                                                        model$kernParam[[j]]$par[c(1, k+1)]))
    Phi.perVar[[j]] <- lapply(1:dblock[j],
                              function(k) basisCompute.lineqGP(x[, partition[[j]][k]],
                                                               subdivision[[j]][[k]]))
    
    names(Gamma.perVar[[j]]) <- names(Phi.perVar[[j]]) <- paste("x", partition[[j]], sep = "")
  }
  
  model$subdivision <- subdivision
  model$Gamma.perVar <- Gamma.perVar
  model$Phi.perVar <- Phi.perVar
  
  # precomputing the linear system for the QP solver and MCMC samplers
  M.perBlock <- g.perBlock <- m.perBlock <- vector("list", nblock)
  names(M.perBlock) <- names(g.perBlock) <- names(m.perBlock) <- paste("x", 1:nblock, sep = "")
  for (j in 1:nblock) {
    if (model$constrType[j] == "linear") { #to be checked!
      if (!("Lambda" %in% names(model)))
        stop('matrix Lambda is not defined')
      Lambda <- model$constrParam[[j]]$Lambda
      lb <- model$constrParam[[j]]$lb
      ub <- model$constrParam[[j]]$ub
      lsys <- lineqGPSys(nrow(Lambda), model$constrType[j], lb, ub,
                         Lambda, lineqSysType = "oneside")
      lsys2 <- lineqGPSys(nrow(Lambda), model$constrType[j], lb, ub,
                          Lambda, rmInf = FALSE)
    } else {
      bounds <- model$constrParam[[j]]$bounds
      
      
      lsys <- lineqGPSys(mlist[[j]], model$constrType[[j]], bounds[1], bounds[2],
                         constrIdx = model$constrIdx, lineqSysType = "oneside")
      lsys2 <- lineqGPSys(mlist[[j]], model$constrType[[j]], bounds[1], bounds[2],
                          constrIdx = model$constrIdx, rmInf = FALSE)
      
      model$constrParam[[j]]$Lambda <- lsys2$A
      model$constrParam[[j]]$lb <- lsys2$l
      model$constrParam[[j]]$ub <- lsys2$u
    }
    
    # oneside linear structure for QP.solver: M = [Lambda,-Lambda] and g = [-lb,ub]
    M.perBlock[[j]] <- lsys$M
    g.perBlock[[j]] <- -matrix(lsys$g)
    # twosides linear structure (Lambda, lb, ub) for MCMC samplers
    
    # extra term required for HMC sampler
    m.perBlock[[j]] <- nrow(lsys2$A)
  }
  # adding the parameters to the model structure
  model$lineqSys$M.perBlock <- M.perBlock # for QP solve
  model$lineqSys$g.perBlock <- g.perBlock # for QP solve
  model$localParam$m.perblock <- m.perBlock # for HMC sampler

  model$lineqSys$M <- eval(parse(text = paste("bdiag(",
                                              paste("M.perBlock[[", 1:nblock, "]]", sep = "", collapse = ","),
                                              ")", sep = "")))
  model$lineqSys$M <- matrix(model$lineqSys$M, ncol = mtotal)
  model$lineqSys$g <- matrix(unlist(model$lineqSys$g.perBlock), ncol = 1)
  
  model$lineqSys$Lambda <- eval(parse(text = paste("bdiag(",
                                                   paste("model$constrParam[[", 1:nblock, "]]$Lambda",
                                                         sep = "", collapse = ","),
                                                      ")", sep = "")))
  model$lineqSys$Lambda <- matrix(model$lineqSys$Lambda, ncol = mtotal)
  model$lineqSys$lb <- eval(parse(text = paste("c(",
                                               paste("model$constrParam[[", 1:nblock, "]]$lb",
                                                     sep = "", collapse = ","),
                                               ")", sep = "")))
  model$lineqSys$ub <- eval(parse(text = paste("c(",
                                               paste("model$constrParam[[", 1:nblock, "]]$ub",
                                                     sep = "", collapse = ","),
                                               ")", sep = "")))

  
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
  mlist <- model$localParam$mlist
  partition <- model$localParam$partition
    
  nt <- length(model$y) # nb of training points 
  mt <- model$localParam$mtotal # total nb of knots
  
  # passing some terms from the model
  pred <- list()
  class(pred) <- class(model)
  pred$constrParam <- model$constrParam
  
  # precomputing some terms
  
  Gamma.perBlock <- list("vector", nblock)
  Phi.perBlock <- list("vector", nblock)
  
  Phi.test.perVar <- vector("list", nblock)
  Phi.test.perBlock <- vector("list", nblock)
  
  for (j in 1:nblock) {
    Gamma.perBlock[[j]] <- eval(parse(text = paste("model$Gamma.perVar[[j]][[", 1:dblock[j], "]]",
                                                   sep = "", collapse = "%x%")))
    
    Phi.perBlock[[j]] <- matrix(0, nt, prod(mlist[[j]]))
    for (n in 1:nt)
      Phi.perBlock[[j]][n, ] <- eval(parse(text = paste("model$Phi.perVar[[j]][[", 1:dblock[j], "]][", n, ", ]",
                                                        sep = "", collapse = "%x%")))
    
    Phi.test.perVar[[j]] <- lapply(1:dblock[j],
                                   function(k) basisCompute.lineqGP(xtest[, partition[[j]][k]],
                                                                    model$subdivision[[j]][[k]]))
    Phi.test.perBlock[[j]] <- matrix(0, nrow(xtest), prod(mlist[[j]]))
    for (n in 1:nrow(xtest))
      Phi.test.perBlock[[j]][n, ] <- eval(parse(text = paste("Phi.test.perVar[[j]][[", 1:dblock[j], "]][", n, ", ]",
                                                            sep = "", collapse = "%x%")))
    
  }
    
  

  # # computing the conditional mean vector and conditional covariance matrix
  # # given the interpolation points
  bigPhi <- eval(parse(text = paste("cbind(",
                                    paste("Phi.perBlock[[", 1:nblock, "]]", sep = "", collapse = ","),
                                    ")", sep = "")))

  bigGamma <- eval(parse(text = paste("bdiag(",
                                      paste("Gamma.perBlock[[", 1:nblock, "]]", sep = "", collapse = ","),
                                      ")", sep = "")))
  
  # ##Function to have stable inverse of symetric matrices
  # stable_inv <-function(Sym){
  #   cholSym <- lapply(Sym, function(x) chol(x))
  #   invSym <- lapply(cholSym, function(x) chol2inv(x))
  # }
  
  cholGamma.perBlock <- lapply(Gamma.perBlock, function(x) chol(x))
  invGamma.perBlock <- lapply(cholGamma.perBlock, function(x) chol2inv(x))
  invBigGamma <- eval(parse(text = paste("bdiag(",
                                         paste("invGamma.perBlock[[", 1:nblock, "]]", sep = "", collapse = ","),
                                         ")", sep = "")))
  invSigma <- as.matrix(invBigGamma) + t(bigPhi)%*%bigPhi/model$varnoise
  
  ### to be adapted !! ###
  ## to fix the matrix inversion lemma in the computation of the mean
  # if (mt < nt) {
  #   PhicholGamma.perBlock <- lapply(1:nblock, function(j) Phi.perBlock[[j]]%*%cholGamma.perBlock[[j]])
  #   PhicholGamma <- eval(parse(text = paste("cbind(", paste("PhicholGamma.perBlock[[", 1:nblock, "]]",
  #                                                           sep = "", collapse = ","), ")", sep = "")))
  #   ILtPhitPhiL <- model$varnoise*diag(mt) + t(PhicholGamma) %*% PhicholGamma
  #   
  #   # ILtPhitPhiL <- ILtPhitPhiL + model$nugget*nrow(ILtPhitPhiL)
  #   cholILtPhitPhiL <- t(chol(ILtPhitPhiL))
  #   Lschur <- forwardsolve(cholILtPhitPhiL, t(PhicholGamma))
  #   invPhiGammaPhitFull <- (diag(nt) - t(Lschur)%*% Lschur)/model$varnoise
  # } else {
    PhiGammaPhit.perBlock <- lapply(1:nblock, function(j) Phi.perBlock[[j]]%*%Gamma.perBlock[[j]]%*%t(Phi.perBlock[[j]]))
    PhiGammaPhit <- eval(parse(text = paste("PhiGammaPhit.perBlock[[", 1:nblock, "]]", sep = "", collapse = "+")))
    # PhiGammaPhit <- as.matrix(bigPhi %*% bigGamma %*% t(bigPhi))
    PhiGammaPhit <- PhiGammaPhit + model$nugget*diag(nrow(PhiGammaPhit))
    invPhiGammaPhitFull <- chol2inv(chol(PhiGammaPhit + model$varnoise * diag(nt))) # instability issues here
  # }

  temp <- bigGamma%*%t(bigPhi) %*% invPhiGammaPhitFull
  predMean <- as.vector(temp %*% model$y)
  pred$mean <- predMean #matrix(c(mu), ncol = 1)
  pred$Sigma <- chol2inv(chol(invSigma))


  pred$mode <- solve.QP(invSigma, t(predMean) %*% invSigma,
                        t(model$lineqSys$M), model$lineqSys$g)$solution

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
  
  etaAll.map <- as.vector(model$lineqSys$Lambda %*% pred$xiAll.map)
  Sigma.etaAll <- model$lineqSys$Lambda %*% pred$Sigma %*% t(model$lineqSys$Lambda)
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
  xiAll.sim <- qr.solve(model$lineqSys$Lambda, etaAll)
  
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

