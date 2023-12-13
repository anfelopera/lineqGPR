#' @title Creation Method for the \code{"lineqBAGP"} S3 Class
#' @description Creation method for the \code{"lineqBAGP"} S3 class.
#' 
#' @param x A vector or matrix with the input data. The dimensions should be indexed by columns.
#' @param y A vector with the output data.
#' @param constrType A character string corresponding to the type of the inequality constraint.
#' @param partition A list containing the partion of set \eqn{\{1, \cdots, D\}}.
#' @param subdivision A list containing sequences of subdivision of segment $[0,1]$.
#' @param subdivision_size A list containing the sizes of the sequences of subdivision of segment $[0,1]$.
# #' @param subdivision A list[list[seq]] The subdivisions for each axe subdivision[[j]][[k]]
# #' @param subdivision_taille A list[seq] subdivision_size[[j]][k] =length(subdivision[[j]][[k]])
# #' @param dim_block :seq, the size of the blocks, dim_block[j]=length(partition[[j]])
# #'  will give the subdivision of $[0,1]$ for the variable partition[[j]][[k]]    
# #' @param X A matrix such that X[,i] is an observation $x^{(i)}$
# #' @param x A float supposed to be a real number in $[0,1]$
# #' @param Y The vector of observation Y[i]=f(X[,i]) 
# #' @param A list["vector"] Matrices per block of type 
# #' @param B list["vector"] Matrices per block of type 
# #' @param d a number corresponding to the dimension of the input space.
# #' 
# #' Options: "boundedness", "monotonicity", "convexity", "linear"
# #' Multiple constraints can be also defined, e.g. \code{constrType = c("boundedness", "monotonicity")}
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
#' @author M. Deronzier and A. F. Lopez-Lopera
#' 
#' @references A. F. Lopez-Lopera (2019),
#' "Gaussian process modelling under inequality constraints".
#' \emph{PhD thesis, Mines Saint-Etienne}
#' <https://tel.archives-ouvertes.fr/tel-02863891>
#'
#' @examples
#' 
#' # synthetic data
#' targetFun <- function(x)  return(x[, 1]*x[, 3] + x[, 2])
#' 
#' d <- 3 # nb of dimensions
#' nblocks <- 2 # nb of blocks
#' 
#' nbdesign <- 4
#' xdesign <- matrix(runif(nbdesign*d), nbdesign, d)
#' ydesign <- targetFun(xdesign)
#' 
#' # creating the model
#' model <- create(class = "lineqBAGP", x = xdesign, y = ydesign,
#'                 constrType = rep("monotonicity", nblocks),
#'                 partition = list(c(1,3), 2),
#'                 subdivision = list(list(seq(0, 1, length = 3), seq(0, 1, length = 5)), 
#'                                    list(c(0, 1))))
#' str(model)
#' 
#' model <- create(class = "lineqBAGP", x = xdesign, y = ydesign,
#'                 constrType = rep("monotonicity", nblocks),
#'                 partition = list(c(1, 3), 2),
#'                 subdivision_size = list(c(3, 5), 2))
#' str(model)
#' 
#' @method create lineqBAGP
#' @export

create.lineqBAGP <- function(x, y, constrType,
                             partition = NULL,
                             subdivision = NULL,
                             subdivision_size = NULL
                             ) {
  # changing the data as matrices
  if (is.null(partition)&&is.null(subdivision)){
    model <- list(x = x, y = y, constrType = constrType,  subdivision = subdivision,
                  partition = partition)
    return(model)
  }
  if (!is.matrix(x)) 
    x <- as.matrix(x)
  if (!is.matrix(y) || ncol(y) != 1) 
    y <- matrix(y)
  
  d <- ncol(x) # dimension of the input space
  nblocks <- length(partition) # nb of blocks
  dim_block <- sapply(partition, function(x) length(x))
  
  if (!is.null(subdivision)) {
    subdivision_size <- lapply(subdivision, function(j) sapply(j, function(k) length(k)))
  } else {
    if (!is.null(subdivision_size)) {
      subdivision <- lapply(subdivision_size, function(j) lapply(j, function(k) seq(0, 1, length = k)))
    } else {
      warning("Either 'subdivision' or 'subdivision_size' must be non null.")
    }
  }
  names(subdivision) <- names(subdivision_size) <- names(partition) <- paste("block", 1:nblocks, sep = "")
  
  #Creation of the subdivision_size from m, type: list[seq]
  
  nknots <- sum(sapply(subdivision_size, prod))
  
  for (j in 1:nblocks) {
    names(partition[[j]]) <- names(subdivision[[j]]) <-
      names(subdivision_size[[j]]) <- paste("var", partition[[j]], sep = "")
  }
  
  # creating some lists for the model
  
  localParam <- list(subdivision_size = subdivision_size,
                     sampler = "ExpT", nblocks = nblocks, dim_block = dim_block, 
                     samplingParams = c(thinning = 1, burn.in = 1, scale = 0.1))
  
  ##  to be adapted to deal with componentwise constraints per blocks 
  # constrFlags <- lapply(subdivision_size, function(sub_temp) rep(1, length(sub_temp)))
  
  kernParam <- list()
  kernParam$type <- "matern52" # to consider different kernels per blocks
  kernParam$par <- constrParam <- vector("list", nblocks)
  names(constrParam) <- paste("block", 1:nblocks, sep = "")
  
  # to initialize the list objects with the correct names
  constrTypePartition <- constrFlag <- partition
  
  
  for (j in 1:nblocks) { # to be checked later! We will focus on the monotonicity constraint
    kernParam$par[[j]] <- c(sigma2 = 1^2, theta = rep(0.1, dim_block[j]))
    
    nvar <- length(partition[[j]])
    constrParam[[j]] <- vector("list", nvar) 
    names(constrParam[[j]]) <- names(partition[[j]])
    for (k in 1:nvar) {
      constrFlag[[j]][k] <- 1
      idxVar <- partition[[j]][k]
      constrTypePartition[[j]][k] <- constrType[idxVar]
  
      switch (constrType[idxVar],
              boundedness = {
                constrParam[[j]][[k]]$bounds <- c(lower = min(y) - 0.05*abs(max(y) - min(y)),
                                             upper = max(y) + 0.05*abs(max(y) - max(y)))
              }, monotonicity = {
                constrParam[[j]][[k]]$bounds <- c(0, Inf)
              }, decreasing = {
                constrParam[[j]][[k]]$bounds <- c(0, Inf)
              }, convexity = {
                constrParam[[j]][[k]]$bounds <- c(0, Inf)
              }, linear = { # to be checked!
                # imposing positiveness constraints by default
                constrParam[[j]][[k]]$Lambda <- diag(prod(model$localParam$m))
                constrParam[[j]][[k]]$lb <- rep(0, nrow(model$Lambda))
                constrParam[[j]][[k]]$ub <- rep(Inf, nrow(model$Lambda))
              }, none = {
                constrParam[[j]][[k]]$bounds <- c(-Inf, Inf)
                constrFlag[[j]][k] <- 0
              }
      )
    }
    
    
  }
  # constrIdx <- which(constrFlags == 1) # to be checked !!
  # creating the full list for the model
  constrIdx <- lapply(constrFlag,  function(x) which(x == 1))
  
  model <- list(x = x, y = y, constrType = constrTypePartition,  subdivision = subdivision,
                partition = partition, d = d,  nugget = 1e-6, nblock = nblocks,
                constrIdx = constrIdx, 
                constrParam = constrParam, nknots = nknots,
                varnoise = 0.05*sd(y)^2,  localParam = localParam, kernParam = kernParam)
  return(model)
}

#' @title Augmenting Method for the \code{"lineqBAGP"} S3 Class
#' @description Augmenting method for the \code{"lineqBAGP"} S3 class.
#' 
#' @param x An object with class \code{lineqBAGP}.
#' @param ... further arguments passed to or from other methods.
#' 
#' @return An expanded \code{"lineqBAGP"} An object with the following additional elements.
#' 
#' \item{Phi}{a matrix corresponding to the hat basis functions.
#' The basis functions are indexed by rows}
#' \item{kernel$Param}{an object indicating for each variable the associated kernel} #will evolve
#' \item{Gamma.var}{a list of list of matrix where Gamma.var[[j]][[k]][i] is 
#' the covariance matrix of the 1D GP with kernel kernel$Param[[j]][[k]] for the points 
#' \eqn{x^{(1)}_{\sigma(i)},..., x^{n}_{\sigma(i)}}}.
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
#' @author M. Deronzier and A. F. Lopez-Lopera
#'
#' @references A. F. Lopez-Lopera (2019),
#' "Gaussian process modelling under inequality constraints".
#' \emph{PhD thesis, Mines Saint-Etienne}
#' <https://tel.archives-ouvertes.fr/tel-02863891>
#' 
#' @examples
#' 
#' # synthetic data
#' targetFun <- function(x)  return(x[, 1]*x[, 3] + x[, 2])
#' 
#' d <- 3 # nb of dimensions
#' nblocks <- 2 # nb of blocks
#' 
#' nbdesign <- 4
#' xdesign <- matrix(runif(nbdesign*d), nbdesign, d)
#' ydesign <- targetFun(xdesign)
#' 
#' # creating the model
#' model <- create(class = "lineqBAGP", x = xdesign, y = ydesign,
#'                 constrType = rep("monotonicity", nblocks),
#'                 partition = list(c(1, 3), 2),
#'                 subdivision_size = list(c(3, 5), 2))
#' str(model)
#'                 
#' # computing some matrices used in the "lineqBAGP" model
#' model <- augment(model)
#' str(model)
#'
#' @importFrom broom augment
#' @export

augment.lineqBAGP <- function(x, ...) {
  model <- x
  if (!("nugget" %in% names(model)))
    model$nugget <- 0
  # if ("bounds" %in% names(model)) {
  #   bounds <- model$bounds
  # } else {
  #   bounds <- c(0, Inf)
  # }
  
  # passing some terms from the model
  x <- model$x
  nblocks <- model$localParam$nblocks
  dim_block <- model$localParam$dim_block
  partition <- model$partition
  subdivision <- model$subdivision
  subdivision_size <- model$localParam$subdivision_size
  nknots <- model$localParam$nknots <- model$nknots 
  
  # Beginning of the construction of the two matrices Gamma and Phi which will 
  # help us to express the law of the Gaussian vector from the conditional law.
  
  Gamma.var <- Gamma_var(subdivision, model$kernParam$par, model$kernParam$type)
  Phi.var <- Phi_per_var(subdivision, partition, x)
  
  model$Gamma.var <- Gamma.var
  model$Phi.var <- Phi.var
  
  # precomputing the linear system for the QP solver and MCMC samplers
  M.block <- g.block <- vector("list", nblocks)
  m.block <- rep(0, nblocks)
  names(M.block) <- names(g.block) <- names(m.block) <- paste("block", 1:nblocks, sep = "")
  for (j in 1:nblocks) {
    if (length(model$constrType[[j]]) == 1 && model$constrType[[j]] == "linear") { # to be checked!
      if (!("Lambda" %in% names(constrParam[[j]])))
        stop('matrix Lambda is not defined')
      Lambda <- model$constrParam[[j]]$Lambda
      lb <- model$constrParam[[j]]$lb
      ub <- model$constrParam[[j]]$ub
      lsys <- lineqGPSys(nrow(Lambda), model$constrType[[j]], lb, ub,
                         Lambda, lineqSysType = "oneside")
      lsys2 <- lineqGPSys(nrow(Lambda), model$constrType[[j]], lb, ub,
                          Lambda, rmInf = FALSE)
    } else {
      nvar <- length(partition[[j]])
      
      bounds <- t(sapply(1:nvar, function(k) model$constrParam[[j]][[k]]$bounds))
      
      # The constraints are imposed over all the variables in a block
      # lsys <- lineqGPSys(subdivision_size[[j]], model$constrType[[j]], bounds[,1], bounds[,2],
      #                    constrIdx = model$constrIdx[[k]], lineqSysType = "oneside")
      
      temp <- sapply(1:nvar,
                     function(k) lineqGPSys(subdivision_size[[j]][k], model$constrType[[j]][k], bounds[k,1], bounds[k,2],
                                            constrIdx = model$constrIdx[[j]][k], lineqSysType = "twosides", rmInf = FALSE))
      
      Abase <- temp[1, ]
      lbase <- temp[2, ]
      ubase <- temp[3, ]
      Idiag <- lapply(as.list(model$localParam$subdivision_size[[j]]), diag)
      ones <- lapply(Idiag, diag)
      
      Anames <- matrix(paste("Idiag[[", seq(nvar),"]]", sep = ""),
                       nvar, nvar, byrow = TRUE)
      diag(Anames) <- paste("Abase[[",seq(nvar),"]]", sep = "")
      lbnames <- matrix(paste("ones[[",seq(nvar),"]]", sep = ""),
                        nvar, nvar, byrow = TRUE)
      diag(lbnames) <- paste("lbase[[",seq(nvar),"]]", sep = "")
      ubnames <- matrix(paste("ones[[",seq(nvar),"]]", sep = ""),
                        nvar, nvar, byrow = TRUE)
      diag(ubnames) <- paste("ubase[[",seq(nvar),"]]", sep = "")
      A <- l <- u <- c()
      for (k in 1:length(partition[[j]])) {
        Atemp <- eval(parse(text = paste(Anames[k, ], collapse = " %x% ")))
        A <- rbind(A, Atemp)
        lbtemp <- eval(parse(text = paste(lbnames[k, ], collapse = " %x% ")))
        l <- c(l, lbtemp)
        ubtemp <- eval(parse(text = paste(ubnames[k, ], collapse = " %x% ")))
        u <- c(u, ubtemp)
      }
      lsys <- bounds2lineqSys(nrow(A), l, u, A, lineqSysType = "oneside")
      lsys2 <- bounds2lineqSys(nrow(A), l, u, A, rmInf = FALSE)
        
      
      # lsys <- lapply(1:nvar, function(k)
      #                  lineqGPSys(subdivision_size[[j]][k], model$constrType[[j]][k], bounds[k,1], bounds[k,2],
      #                             constrIdx = model$constrIdx[[j]][k], lineqSysType = "oneside"))
      # lsys2 <-  lapply(1:nvar, 
      #                  function(k) lineqGPSys(subdivision_size[[j]][k], model$constrType[[j]][k], bounds[k,1], bounds[k,2],
      #                                                 constrIdx = model$constrIdx[[j]][k], rmInf = FALSE))
      
      
      model$constrParam[[j]]$Lambda <- lsys2$A
      model$constrParam[[j]]$lb <- lsys2$l
      model$constrParam[[j]]$ub <- lsys2$u
    }
    
    # oneside linear structure for QP.solver: M = [Lambda,-Lambda] and g = [-lb,ub]
    M.block[[j]] <- lsys$M
    g.block[[j]] <- -matrix(lsys$g)
    # twosides linear structure (Lambda, lb, ub) for MCMC samplers
    
    # extra term required for HMC sampler
    m.block[j] <- nrow(lsys2$A)
  }
  # adding the parameters to the model structure
  model$lineqSys$M.block <- M.block # for QP solve
  model$lineqSys$g.block <- g.block # for QP solve
  model$localParam$m.block <- m.block # for HMC sampler

  model$lineqSys$M <- eval(parse(text = paste("bdiag(",
                                              paste("M.block[[", 1:nblocks, "]]", sep = "", collapse = ","),
                                              ")", sep = "")))
  model$lineqSys$M <- matrix(model$lineqSys$M, ncol = nknots)
  model$lineqSys$g <- matrix(unlist(model$lineqSys$g.block), ncol = 1)
  
  model$lineqSys$Lambda <- eval(parse(text = paste("bdiag(",
                                                   paste("model$constrParam[[", 1:nblocks, "]]$Lambda",
                                                         sep = "", collapse = ","),
                                                      ")", sep = "")))
  model$lineqSys$Lambda <- matrix(model$lineqSys$Lambda, ncol = nknots)
  model$lineqSys$lb <- eval(parse(text = paste("c(",
                                               paste("model$constrParam[[", 1:nblocks, "]]$lb",
                                                     sep = "", collapse = ","),
                                               ")", sep = "")))
  model$lineqSys$ub <- eval(parse(text = paste("c(",
                                               paste("model$constrParam[[", 1:nblocks, "]]$ub",
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
#' \item{Gamma}{ Matrix representing \eqn{k_{S}(\xi,\xi)}}
#' \item{Gamma.block}{ Gamma represented per block as a list of matrices}
#' \item{Phi}{ Matrix representing \eqn{\Phi}}
#' \item{Phi.block}{ Phi represented per block as a list of matrices}
#' \item{t_phi}{ The transpose of Phi}
#' \item{t_Phi.block}{ t_Phi represented per block as a list of matrices}
#' \item{invSigma}{ The covariance matrix of the gaussian vector }
#' \item{y.mean}{Prediction of x knowing that }
#' \item{mid.term}{ The mid term is computed}
#' \item{xi.mode}{Mode of the truncated Gaussian vector}
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
#' @author M. Deronzier and A. F. Lopez-Lopera
#'
#' @references A. F. Lopez-Lopera (2019),
#' "Gaussian process modelling under inequality constraints".
#' \emph{PhD thesis, Mines Saint-Etienne}
#' <https://tel.archives-ouvertes.fr/tel-02863891>
#'
#' @examples
#' 
#' # synthetic data
#' d <- 3 # number of active input variables
#' partition <- list(c(1,3), 2) # partition of the block structure
#' nblocks <- length(partition) # nb of blocks
#' 
#' targetFun <- function(x, partition)
#'   return(x[, partition[[1]][1]]*x[, partition[[1]][2]] + x[, partition[[2]][1]])
#' 
#' # building a random design 
#' nbdesign <- 6*d
#' xdesign <- matrix(runif(nbdesign*d), nbdesign, d)
#' ydesign <- targetFun(xdesign, partition)
#' 
#' # defining the 3D grid for predictions 
#' n1D <- 20
#' xbase <- seq(0, 1, length = n1D)
#' xtest <- as.matrix(expand.grid(xbase, xbase, xbase))
#' ytest <- targetFun(xtest, partition)
#'
#' # creating the model
#' model <- create(class = "lineqBAGP", x = xdesign, y = ydesign,
#'                 constrType = rep("monotonicity", nblocks), 
#'                 partition = partition,
#'                 subdivision_size = list(c(3, 5), 2))
#'
#' # modifying the covariance parameters of each block
#' for (k in 1:nblocks)
#'   model$kernParam$par[[k]] <- c(1, rep(0.1, model$localParam$dim_block[k]))
#'
#' # computing the unconstrained GP mean and the constrained GP mode
#' pred <- predict(model, xtest)
#'
#' # Q2 criterion
#' var_design <- mean((ytest- mean(ytest))^2)
#' message("Unconstrained GP mean: ", 1 - mean((ytest- pred$y.mean)^2)/var_design)
#' message("Constrained GP mode: ", 1 - mean((ytest- pred$y.mode)^2)/var_design)
#' 
#' @importFrom quadprog solve.QP
#' @importFrom Matrix bdiag
#' @import plot3D
#' @export
predict.lineqBAGP <- function(object, xtest, return_model = FALSE, ...) {
  model <- augment(object)
  if (!is.matrix(xtest))
    xtest <- matrix(xtest, ncol = model$d)
  
  nblocks <- model$localParam$nblocks
  dim_block <- model$localParam$dim_block
  subdivision <- model$subdivision
  subdivision_size <- model$localParam$subdivision_size
  block_tensor_size <- sapply(subdivision_size, prod)
  partition <- model$partition
  inv_tau <- 1/model$varnoise
    
  nobs <- length(model$y) # nb of training points 
  nknots <- model$localParam$nknots # total nb of knots
  
  # passing some terms from the model
  pred <- list()
  class(pred) <- class(model)
  pred$constrParam <- model$constrParam
  
  # The real block matrices
  
  Gamma.block <- Gamma_var_to_tensor(model$Gamma.var)
  Gamma <- block_to_matrix(Gamma.block, "bdiag")
  #nugget.block <- lapply(1:nblocks, function(x) 1e-9*diag(block_tensor_size[x]))
  #Gamma.block <- block_compute(Gamma.block, "sum", nugget.block)
  
  Phi.block <- Phi_var_to_tensor(model$Phi.var)
  Phi <- block_to_matrix(Phi.block, "cbind")
  t_Phi.block <- block_compute(Phi.block, "transpose")
  t_Phi <- block_to_matrix(t_Phi.block, "rbind")
    
  # One Block matrix for testing our results
  
  #Phi.test.var <- Phi_per_var(subdivision,xtest)
  #Phi.test.block <- Phi_var_to_tensor(Phi.test.var, subdivision)
  pred$Phi.test <- Phi.test <- block_to_matrix(Block_Phi(subdivision, partition, xtest), "cbind")
  
  
  #The cholesky decomposition to stabilize the inverse operation
  
  cholGamma.block <- block_compute(Gamma.block, "chol")
  invGamma.block <- block_compute(cholGamma.block, "chol2inv")
  
  #Computation of the conditional covariance matrix of the vector 
  t_PhiPhi <- t_Phi%*%Phi
  invSigma <- block_to_matrix(invGamma.block, "bdiag") + inv_tau*t_PhiPhi
  
  #Computation of the mean vector for the conditional vector law, according to the 
  #complexity study we have two cases to optimize the computation
  
  In <- diag(nobs)  
  Gammat_Phi.block <- block_compute(Gamma.block, "prod",  t_Phi.block)
  Gammat_Phi <- block_to_matrix(Gammat_Phi.block, "rbind")
  
  #Computation of the inverse of the mid term in the most efficient way
  if (nobs<=2*nknots) {
    mid.term <- chol2inv(chol(block_to_matrix(block_compute(Phi.block, "prod", Gammat_Phi.block), "sum") + 
                                model$varnoise*In))
  } else {
    mid.term <- inv_tau*(In-inv_tau*Phi %*%
                           chol2inv(chol(block_to_matrix(invGamma.block) + inv_tau*t_PhiPhi)) %*% t_Phi)
    #message("computation of mu using Woodbury formula")
  }
  Gammat_Phimid.term <- Gammat_Phi %*% mid.term
  xi.mean <- Gammat_Phimid.term %*% model$y
  pred$xi.mean <- xi.mean
  pred$Sigma <- Gamma - Gammat_Phimid.term %*% t(Gammat_Phi)
  
  # pred$Sigma <- chol2inv(chol(invSigma))

  pred$xi.mode <- as.matrix(solve.QP(invSigma, t(pred$xi.mean) %*% invSigma,
                           t(model$lineqSys$M), model$lineqSys$g)$solution)
  
  
  pred$y.mean <- pred$Phi.test%*%pred$xi.mean
  pred$y.mode <- pred$Phi.test%*%pred$xi.mode

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
#' @author M. Deronzier and A. F. Lopez-Lopera
#'
#' @references A. F. Lopez-Lopera (2019),
#' "Gaussian process modelling under inequality constraints".
#' \emph{PhD thesis, Mines Saint-Etienne}
#' <https://tel.archives-ouvertes.fr/tel-02863891>
#'
#' @examples
#' 
#' # synthetic data
#' d <- 3 # number of active input variables
#' partition <- list(c(1,3), 2) # partition of the block structure
#' nblocks <- length(partition) # nb of blocks
#' 
#' targetFun <- function(x, partition)
#'   return(x[, partition[[1]][1]]*x[, partition[[1]][2]] + x[, partition[[2]][1]])
#' 
#' # building a random design 
#' nbdesign <- 6*d
#' xdesign <- matrix(runif(nbdesign*d), nbdesign, d)
#' ydesign <- targetFun(xdesign, partition)
#' 
#' # defining the 3D grid for predictions 
#' n1D <- 20
#' xbase <- seq(0, 1, length = n1D)
#' xtest <- as.matrix(expand.grid(xbase, xbase, xbase))
#' ytest <- targetFun(xtest, partition)
#'
#' # creating the model
#' model <- create(class = "lineqBAGP", x = xdesign, y = ydesign,
#'                 constrType = rep("monotonicity", nblocks), 
#'                 partition = partition,
#'                 subdivision_size = list(c(3, 5), 2))
#'
#' # modifying the covariance parameters of each block
#' for (k in 1:nblocks)
#'   model$kernParam$par[[k]] <- c(1, rep(0.1, model$localParam$dim_block[k]))
#'   
#' # simulating constrained MCMC samples
#' model.sim <- simulate(model, 1e2, seed = 1, xtest)
#'
#' # Q2 criterion
#' var_design <- mean((ytest- mean(ytest))^2)
#' message("Unconstrained GP mean: ", 1 - mean((ytest- model.sim$y.mean)^2)/var_design)
#' message("Constrained GP mode: ", 1 - mean((ytest- model.sim$y.mode)^2)/var_design)
#' message("Constrained GP mean via MCMC: ", 1 - mean((ytest- rowMeans(model.sim$y.sim))^2)/var_design)
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
  
  nblocks <- model$localParam$nblocks
  
  eta.mode <- as.vector(model$lineqSys$Lambda %*% pred$xi.mode)
  Sigma.eta <- model$lineqSys$Lambda %*% pred$Sigma %*% t(model$lineqSys$Lambda)
  Sigma.eta <- Sigma.eta + model$nugget*diag(nrow(Sigma.eta))
  
  # listing control terms
  control <- as.list(unlist(model$localParam$samplingParam))
  control$mvec <- sum(model$localParam$m.block) # for HMC
  control$constrType <- model$constrType # for HMC
  
  tmvPar <- list(mu = eta.mode, Sigma = Sigma.eta,
                 lb = model$lineqSys$lb,
                 ub = model$lineqSys$ub)
  class(tmvPar) <- model$localParam$sampler
  set.seed(seed)
  simtime <- proc.time()
  eta.sim <- tmvrnorm(tmvPar, nsim, control)
  simtime <- proc.time() - simtime
  xi.sim <- qr.solve(model$lineqSys$Lambda, eta.sim)
  
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
  
  simModel$xi.mean <- pred$xi.mean
  simModel$xi.mode <- pred$xi.mode
  simModel$y.mean <- pred$y.mean
  simModel$y.mode <- pred$y.mode
  simModel$xi.sim <- xi.sim
  simModel$y.sim <- pred$Phi.test %*% xi.sim
  # simModel$PhiAll.test <- pred$PhiAll.test
  
  class(simModel) <- class(model)
  return(simModel)
}

