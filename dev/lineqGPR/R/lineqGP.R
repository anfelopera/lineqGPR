#' @title Hat Basis Functions for \code{"lineqGP"} Models
#' @description Evaluate the hat basis functions for \code{"lineqGP"} models.
#' 
#' @param x a vector (or matrix) with the input data
#' @param u a vector (or matrix) with the locations of the knots
#' @param d a number corresponding to the dimension of the input space
#' 
#' @return A matrix with the hat basis functions. The basis functions are indexed by rows
#'
#' @section Comments:
#' This function was tested mainly for 1D or 2D input spaces. It could change
#' in future versions for higher dimensions.
#' 
#' @author A. F. Lopez-Lopera
#'
#' @references A. F. Lopez-Lopera, F. Bachoc, N. Durrande and O. Roustant (2017),
#' "Finite-dimensional Gaussian approximation with linear inequality constraints".
#' \emph{SIAM/ASA Journal on Uncertainty Quantification}, 6(3): 1224–1255.
#' <doi:10.1137/17M1153157>
#'
#' @examples
#' x <- seq(0, 1, 1e-3)
#' m <- 5
#' u <- seq(0, 1, 1/(m-1))
#' Phi <- basisCompute.lineqGP(x, u, d = 1)
#' matplot(Phi, type = "l", lty = 2, main = "Hat basis functions with m = 5")
#'
#' @export
basisCompute.lineqGP <- function(x, u, d = 1) {
  if (!is.matrix(x) || ncol(x) != d)
    x <- matrix(x, ncol = d)
  if (is.list(u) && length(u) == d) {
    m <- as.integer(lapply(u, length))
  } else if (is.double(u)) {
    u <- matrix(u, ncol = d)
    m <- nrow(u)
  } else {
    m <- nrow(u)
  }
  
  # precomputing some terms
  n <- nrow(x)
  # delta <- diff(u)
  
  # computing the hat basis functions
  if (d == 1){
    Phi <- matrix(0, n, m)
    for (i in 1:length(x)) {
      for (j in 1:length(u)) {
        if (x[i] > u[j] + 1e-12) { # there is a numerical stability issue
          distAbs <- abs((x[i] - u[j])/(u[j+1] - u[j]))
        } else if (x[i] < u[j]) {
          distAbs <- abs((x[i] - u[j])/(u[j] - u[j-1]))
        } else {
          distAbs <- 0
        }
        if (distAbs <= 1)
          Phi[i, j] <- 1 - distAbs
      }
    }
  } else if (d >= 2){
    PhiList <- list()
    for (k in seq(d)) {
      PhiList[[k]] <- basisCompute.lineqGP(x[, k], u[[k]], d = 1)
    }
    Phi <- matrix(0, n, prod(m))
    for (i in seq(n)) {
      Phi_temp2 <- PhiList[[1]][i, ]
      for (k in 2:d) {
        Phi_temp2 <- Phi_temp2 %x% PhiList[[k]][i, ]
      }
      Phi[i, ] <- Phi_temp2
    }
  }
  return(Phi)
}

#' @title Linear Systems of Inequalities for \code{"lineqGP"} Models
#' @description Build the linear system of inequalities for \code{"lineqGP"} models.
#' 
#' @param m the number of linear inequality constraints
#' @param constrType a character string corresponding to the type of the inequality constraint
#' Options: "boundedness", "monotonicity", "convexity", "linear"
#' @param l the value (or vector) with the lower bound
#' @param u the value (or vector) with the upper bound
#' @param A a matrix containing the structure of the linear equations
#' @param d the value with the input dimension
#' @param lineqSysType a character string corresponding to the type of the
#' linear system. Options: \code{twosides}, \code{oneside} (see \code{\link{bounds2lineqSys}} for more details)
#' @param constrIdx for d > 1, a logical vector with the indices of active constrained dimensions
#' @param rmInf If \code{TRUE}, inactive constraints are removed
#' @param knots_pos Position of the knots (required for some constraints such convexity)
#' (e.g. \eqn{-\infty \leq x \leq \infty}{-Inf \le x \le Inf}).
#' 
#' @return  A list with the linear system of inequalities: \code{list(A,l,u)} (\code{twosides}) or \code{list(M,g)} (\code{oneside}).
#'
#' @section Comments:
#' This function could change in future versions for more types of inequality
#' constraints in higher dimensions.
#' 
#' @seealso \code{\link{bounds2lineqSys}}
#' 
#' @author A. F. Lopez-Lopera
#'
#' @examples
#' linSys1 <- lineqGPSys(m = 5, constrType = "boundedness", l = 0, u = 1, lineqSysType = "twosides")
#' linSys1
#' linSys2 <- lineqGPSys(m = 5, constrType = "boundedness", l = 0, u = 1, lineqSysType = "oneside")
#' linSys2
#'
#' @export
lineqGPSys <- function(m = nrow(A),
                       constrType = c("boundedness", "monotonicity",
                                      "decreasing", "convexity", "linear", "none"),
                       l = -Inf, u = Inf,  A = diag(m), d = length(m),
                       lineqSysType = c("twosides", "oneside"),
                       constrIdx = seq(length(m)),
                       rmInf = TRUE,
                       knots_pos = NULL
                       ) {
  constrType <- match.arg(constrType)
  lineqSysType <- match.arg(lineqSysType)
  
  if (!is.null(knots_pos) & is.list(knots_pos) & length(knots_pos) != d)
    stop("The length of 'knots_pos' has to be equal to the input dimension")
  
  if (constrType != "linear") {
    if (d == 1) {
      switch(constrType,
             boundedness = {
               linSys <- bounds2lineqSys(m, l, u, lineqSysType = lineqSysType, rmInf = rmInf)
             }, monotonicity = {
               if (length(l) == 1) l <- c(-Inf, rep(l, m-1))
               if (length(u) == 1) u <- c(Inf, rep(u, m-1))
               A <- diag(m)
               if (m == 2) { # a bug with only 2 knots
                 A[2, 1] <- -1
               } else {
                 diag(A[-1, -ncol(A)]) <- -1
               }
               linSys <- bounds2lineqSys(nrow(A), l, u, A, lineqSysType, rmInf)
             }, decreasing = {
               if (length(l) == 1) l <- c(-Inf, rep(l, m-1))
               if (length(u) == 1) u <- c(Inf, rep(u, m-1))
               A <- -diag(m)
               if (m == 2) { # a bug with only 2 knots
                 A[2, 1] <- 1
               } else {
                 diag(A[-1, -ncol(A)]) <- 1
               }
               A[1,1] <- 1
               linSys <- bounds2lineqSys(nrow(A), l, u, A, lineqSysType, rmInf)
             }, convexity = {
               if (length(l) == 1) l <- c(-rep(Inf, 2), rep(l, m-2))
               if (length(u) == 1) u <- c(rep(Inf, 2), rep(u, m-2))
               A <- diag(m)
               # if (m == 3) { # a bug with only 3 knots
               #   A[3, 1:2] <- c(1, -2)
               # } else {
               #   diag(A[-seq(2), -c(ncol(A)-1,ncol(A))]) <- 1
               #   diag(A[-seq(2), -c(1,ncol(A))]) <- -2
               # }
               
               distKnots <- diff(knots_pos)
               for (i in 3:(m)) {
                 A[i, i-2] <- 1 / distKnots[i-2]
                 A[i, i] <-  1 / distKnots[i-1]
                 A[i, i-1] <- - A[i, i-2] - A[i, i]
               }
               linSys <- bounds2lineqSys(nrow(A), l, u, A, lineqSysType, rmInf)
             }, none = {
               linSys <- bounds2lineqSys(m, l, u, lineqSysType = lineqSysType, rmInf = rmInf)
             }
      )
    } else {
      # For the of case convexity for d >= 2, constrType == "convexity" is a
      # weak version of convexity
      temp <- sapply(1:length(m), function(i) lineqGPSys(m[i], constrType, l, u,
                     lineqSysType = "twosides", rmInf = FALSE, knots_pos = knots_pos[[i]]) )
      Abase <- temp[1, ]
      lbase <- temp[2, ]
      ubase <- temp[3, ]
      Idiag <- lapply(as.list(m), diag)
      ones <- lapply(Idiag, diag)
      
      Anames <- matrix(paste("Idiag[[",seq(d),"]]", sep = ""),
                       d, d, byrow = TRUE)
      diag(Anames) <- paste("Abase[[",seq(d),"]]", sep = "")
      lbnames <- matrix(paste("ones[[",seq(d),"]]", sep = ""),
                        d, d, byrow = TRUE)
      diag(lbnames) <- paste("lbase[[",seq(d),"]]", sep = "")
      ubnames <- matrix(paste("ones[[",seq(d),"]]", sep = ""),
                        d, d, byrow = TRUE)
      diag(ubnames) <- paste("ubase[[",seq(d),"]]", sep = "")
      A <- l <- u <- c()
      
      # if (length(constrIdx) < d) {
        for (k in constrIdx) {
          Atemp <- eval(parse(text = paste(Anames[k, ], collapse = " %x% ")))
          A <- rbind(A, Atemp)
          lbtemp <- eval(parse(text = paste(lbnames[k, ], collapse = " %x% ")))
          l <- c(l, lbtemp)
          ubtemp <- eval(parse(text = paste(ubnames[k, ], collapse = " %x% ")))
          u <- c(u, ubtemp)
        }
      # }
      linSys <- bounds2lineqSys(nrow(A), l, u, A, lineqSysType, rmInf)
    }
  } else {
    linSys <- bounds2lineqSys(nrow(A), l, u, A, lineqSysType, rmInf)
  }
  return(linSys)
}

#' @title Creation Method for the \code{"lineqGP"} S3 Class
#' @description Creation method for the \code{"lineqGP"} S3 class.
#' 
#' @param x a vector or matrix with the input data. The dimensions should be indexed by columns
#' @param y a vector with the output data
#' @param constrType a character string corresponding to the type of the inequality constraint.
#' @param m a scalar or vector corresponding to the number of knots per dimension.
#' Options: "boundedness", "monotonicity", "convexity", "linear";
#' Multiple constraints can be also defined, e.g. \code{constrType = c("boundedness", "monotonicity")}
#' 
#' @return A list with the following elements
#' \item{x,y,constrType}{see \bold{Arguments}}
#' \item{d}{a number corresponding to the input dimension}
#' \item{constrIdx}{for d > 1, a logical vector with the indices of active constrained dimensions.}
#' \item{localParam}{a list with specific parameters required for \code{"lineqGP"} models:
#' \code{m} (number of basis functions), \code{sampler}, and \code{samplingParams}.
#' See \code{\link{simulate.lineqGP}}.}
#' \item{kernParam}{a list with the kernel parameters: \code{par} (kernel parameters), \code{type}, \code{nugget}.
#' See \code{\link{kernCompute}}}
#' \item{bounds}{the limit values if \code{constrType = "boundedness"}}
#' \item{(Lambda,lb,ub)}{the linear system of inequalities if \code{constrType = "linear"}}
#'
#' @seealso \code{\link{augment.lineqGP}}, \code{\link{predict.lineqGP}}, \code{\link{simulate.lineqGP}}
#' 
#' @author A. F. Lopez-Lopera
#' 
#' @references A. F. Lopez-Lopera, F. Bachoc, N. Durrande and O. Roustant (2017),
#' "Finite-dimensional Gaussian approximation with linear inequality constraints".
#' \emph{SIAM/ASA Journal on Uncertainty Quantification}, 6(3): 1224–1255.
#' <doi:10.1137/17M1153157>
#'
#' @examples
#' # creating the model
#' sigfun <- function(x) return(1/(1+exp(-7*(x-0.5))))
#' x <- seq(0, 1, length = 5)
#' y <- sigfun(x)
#' model <- create(class = "lineqGP", x, y, constrType = "monotonicity")
#' model
#'
#' @method create lineqGP
#' @export
create.lineqGP <- function(x, y, constrType, m = NULL) {
  # changing the data as matrices
  if (!is.matrix(x))
    x <- as.matrix(x)
  if (!is.matrix(y) || ncol(y) != 1)
    y <- matrix(y)
  
  d <- ncol(x) # dimension of the input space
  
  if (is.null(m)) {
    m <- rep(10*length(y)/d, d)
  } else if (length(m) == 1) {
    m <- rep(m, d)
  } else if (length(m) != d) {
    stop("The length of 'm' has to be equal to the input dimension")
  }
  
  u <- vector("list", d)
  for (k in 1:d)
    u[[k]] <- matrix(seq(0, 1, by = 1/(m[k]-1)), ncol = 1) # discretization vector
  
  # creating some lists for the model
  kernParam <- list(par = c(sigma2 = 1^2, theta = rep(0.2, d)),
                    type = "gaussian", nugget = 1e-7*sd(y)) # Matern52 if MaxMod
  localParam <- list(m = m, sampler = "ExpT",
                     samplingParams = c(thinning = 1, burn.in = 1, scale = 0.1))
  # creating the full list for the model
  model <- list(x = x, y = y, constrType = constrType,  ulist = u,
                d = d, constrIdx = seq(d), varnoise = 0,
                localParam = localParam, kernParam = kernParam, varnoise = 0)
  model$bounds <- c()
  
  if (any(constrType == "none")) {
    model$bounds <- c(-Inf, Inf)
    # constrFlags[k] <- 0
  } else {
    if (any(constrType == "boundedness"))
      model$bounds <- rbind(model$bounds,
                            c(lower = min(y) - 0.05*abs(max(y) - min(y)),
                              upper = max(y) + 0.05*abs(max(y) - max(y))))
    if (any(constrType == "monotonicity"))
      model$bounds <- rbind(model$bounds, c(0, Inf))
    if (any(constrType == "decreasing"))
      model$bounds <- rbind(model$bounds, c(0, Inf))
    if (any(constrType == "convexity"))
      model$bounds <- rbind(model$bounds, c(0, Inf))
  }
  
  if (any(constrType == "linear")) {
    # imposing positiveness constraints by default
    model$Lambda <- diag(prod(model$localParam$m))
    model$lb <- rep(0, nrow(model$Lambda))
    model$ub <- rep(Inf, nrow(model$Lambda))
  }
  return(model)
}

#' @title Augmenting Method for the \code{"lineqGP"} S3 Class
#' @description Augmenting method for the \code{"lineqGP"} S3 class.
#' 
#' @param x an object with class \code{lineqGP}
#' @param ... further arguments passed to or from other methods
#' 
#' @return An expanded \code{"lineqGP"} object with the following additional elements
#' \item{Phi}{a matrix corresponding to the hat basis functions
#' The basis functions are indexed by rows}
#' \item{Gamma}{the covariance matrix of the Gassian vector \eqn{\boldsymbol{\xi}}{\xi}.}
#' \item{(Lambda,lb,ub)}{the linear system of inequalities}
#' \item{...}{further parameters passed to or from other methods}
#'
#' @details Some paramaters of the finite-dimensional GP with linear inequality
#' constraints are computed. Here, \eqn{\boldsymbol{\xi}}{\xi} is a centred Gaussian
#' vector with covariance \eqn{\boldsymbol{\Gamma}}{\Gamma}, s.t.
#' \eqn{\boldsymbol{\Phi} \boldsymbol{\xi} = \boldsymbol{y}}{\Phi \xi = y}
#' (interpolation constraints) and
#' \eqn{\boldsymbol{l} \leq \boldsymbol{\Lambda} \boldsymbol{\xi} \leq \boldsymbol{u}}{lb \le \Lambda \xi \le ub}
#' (inequality constraints).
#' 
#' @seealso \code{\link{create.lineqGP}}, \code{\link{predict.lineqGP}},
#'          \code{\link{simulate.lineqGP}}
#'          
#' @author A. F. Lopez-Lopera
#'
#' @references A. F. Lopez-Lopera, F. Bachoc, N. Durrande and O. Roustant (2017),
#' "Finite-dimensional Gaussian approximation with linear inequality constraints".
#' \emph{SIAM/ASA Journal on Uncertainty Quantification}, 6(3): 1224–1255.
#' <doi:10.1137/17M1153157>
#'
#' @examples
#' # creating the model
#' sigfun <- function(x) return(1/(1+exp(-7*(x-0.5))))
#' x <- seq(0, 1, length = 5)
#' y <- sigfun(x)
#' model <- create(class = "lineqGP", x, y, constrType = "monotonicity")
#'
#' # updating and expanding the model
#' model$localParam$m <- 30
#' model$kernParam$par <- c(1, 0.2)
#' model2 <- augment(model)
#' image(model2$Gamma, main = "covariance matrix")
#'
#' @importFrom broom augment
#' @export
augment.lineqGP <- function(x, ...) {
  model <- x
  if (!("nugget" %in% names(model$kernParam)))
    model$kernParam$nugget <- 0
  if ("bounds" %in% names(model)) {
    bounds <- model$bounds
  } else {
    bounds <- c(0, Inf)
  }
  # passing some terms from the model
  x <- model$x
  
  # computing the kernel matrix for the prior
  if (model$d == 1) {
    u <- model$ulist # discretization vector
    if (is.list(model$ulist))
      u <- model$ulist[[1]]
    m <- model$localParam$m <- length(u)
    Gamma <- kernCompute(u, u, model$kernParam$type, model$kernParam$par)
  } else if (model$d >= 2) {
    m <- model$localParam$m <- sapply(model$ulist, length)
    u <- model$ulist
    AllGammas <- vector("list", model$d) # list with the d covariance matrices
    if (length(model$kernParam$par) == 2 && length(unique(m)) == 1) {
      uBase <- u[[1]]
      GammaBase <- kernCompute(uBase, uBase,
                               model$kernParam$type,
                               model$kernParam$par)
      for (j in 1:model$d) {
        AllGammas[[j]] <- GammaBase
      }
    } else {
      if (length(model$kernParam$par[-1]) != model$d) {
        model$kernParam$par <- c(model$kernParam$par[1],
                                 rep(model$kernParam$par[2], model$d))
        names(model$kernParam$par) <- c(names(model$kernParam$par)[1],
                                        paste(names(model$kernParam$par)[2],
                                              seq(model$d), sep = ""))
      }
      for (j in 1:model$d) {
        AllGammas[[j]] <- kernCompute(u[[j]], u[[j]],
                                      model$kernParam$type,
                                      model$kernParam$par[c(1, j+1)])
      }
    }
    expr <- paste("AllGammas[[", seq(model$d), "]]", collapse = " %x% ")
    Gamma <- eval(parse(text = expr))
  }
  model$Gamma <- Gamma# + model$kernParam$nugget*diag(nrow(Gamma))
  model$Phi <- basisCompute.lineqGP(x, u, ncol(x))
  
  # precomputing the linear system for the QP solver and MCMC samplers
  mt <- nrow(model$Gamma)
  M <- g <- numeric()
  Lambda <- lb <- ub <- numeric()
  mvec <- rep(0, length(model$constrType))
  for (i in seq(length(model$constrType))) {
    if (model$constrType[i] == "linear") {
      if (!("Lambda" %in% names(model)))
        stop('matrix Lambda is not defined')
      lsys <- lineqGPSys(nrow(model$Lambda), model$constrType[i], model$lb, model$ub,
                         model$Lambda, lineqSysType = "oneside")
      lsys2 <- lineqGPSys(nrow(model$Lambda), model$constrType[i], model$lb, model$ub,
                          model$Lambda, rmInf = FALSE)
    } else {
      bounds <- matrix(bounds, ncol = 2)
      lsys <- lineqGPSys(m, model$constrType[i], bounds[i,1], bounds[i,2],
                         d = model$d, constrIdx = model$constrIdx,
                         lineqSysType = "oneside", knots_pos = u)
      lsys2 <- lineqGPSys(m, model$constrType[i], bounds[i,1], bounds[i,2],
                          d = model$d, constrIdx = model$constrIdx,
                          rmInf = FALSE, knots_pos = u)
    }
    # oneside linear structure for QP.solver: M = [Lambda,-Lambda] and g = [-lb,ub]
    M <- rbind(M, lsys$M)
    g <- rbind(g, -matrix(lsys$g))
    # twosides linear structure (Lambda, lb, ub) for MCMC samplers
    Lambda <- rbind(Lambda, lsys2$A)
    lb <- c(lb, lsys2$l)
    ub <- c(ub, lsys2$u)
    # extra term required for HMC sampler
    mvec[i] <- nrow(lsys2$A)
  }
  # adding the parameters to the model structure
  model$Lambda <- Lambda
  model$lb <- lb
  model$ub <- ub
  model$lineqSys$M <- M # for QP solve
  model$lineqSys$g <- g # for QP solve
  model$localParam$mvec <- mvec # for HMC sampler
  model$localParam$u <- u
  return(model)
}

#' @title Prediction Method for the \code{"lineqGP"} S3 Class
#' @description Prediction method for the \code{"lineqGP"} S3 class.
#' 
#' @param object an object with class \code{"lineqGP"}
#' @param xtest a vector (or matrix) with the test input design
#' @param return_model If \code{TRUE}, the augmented model is returned (see \code{\link{augment.lineqGP}}).
#' @param ... further arguments passed to or from other methods
#' @return A \code{"lineqGP"} object with the following elements
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
#' @seealso \code{\link{create.lineqGP}}, \code{\link{augment.lineqGP}},
#'          \code{\link{simulate.lineqGP}}
#'          
#' @author A. F. Lopez-Lopera
#'
#' @references A. F. Lopez-Lopera, F. Bachoc, N. Durrande and O. Roustant (2017),
#' "Finite-dimensional Gaussian approximation with linear inequality constraints".
#' \emph{SIAM/ASA Journal on Uncertainty Quantification}, 6(3): 1224–1255.
#' <doi:10.1137/17M1153157>
#'
#' @examples
#' # creating the model
#' sigfun <- function(x) return(1/(1+exp(-7*(x-0.5))))
#' x <- seq(0, 1, length = 5)
#' y <- sigfun(x)
#' model <- create(class = "lineqGP", x, y, constrType = "monotonicity", m = 30)
#'
#' # modifying the covariance parameters
#' model$kernParam$par <- c(1, 0.2)
#'
#' # predictions from the model
#' xtest <- seq(0, 1, length = 100)
#' ytest <- sigfun(xtest)
#' pred <- predict(model, xtest)
#' plot(xtest, ytest, type = "l", lty = 2, main = "Kriging predictions")
#' lines(xtest, pred$Phi.test %*% pred$mu, type = "l", col = "blue")
#' lines(xtest, pred$Phi.test %*% pred$xi.map, type = "l", col = "red")
#' legend("right", c("ytest", "mean", "mode"), lty = c(2,1,1),
#'        col = c("black", "blue", "red"))
#'
#' @importFrom quadprog solve.QP
#' @export
predict.lineqGP <- function(object, xtest, return_model = FALSE, ...) {
  model <- augment(object)
  if (!is.matrix(xtest))
    xtest <- matrix(xtest, ncol = model$d)

  # passing some terms from the model
  pred <- list()
  class(pred) <- class(model)
  pred$Lambda <- model$Lambda
  pred$lb <- model$lb
  pred$ub <- model$ub

  # precomputing some terms
  pred$Phi.test <- basisCompute.lineqGP(xtest, model$localParam$u, d = model$d)

  # # computing the linear system for the QP solver
  # A <- rbind(model$Phi, model$lineqSys$M)
  # b <- rbind(matrix(model$y, ncol = 1), model$lineqSys$g)

  # computing the conditional mean vector and conditional covariance matrix
  # given the interpolation points
  GammaPhit <- model$Gamma %*% t(model$Phi)
  PhiGammaPhit <- model$Phi %*% GammaPhit 
  PhiGammaPhit <- PhiGammaPhit + (model$varnoise + model$kernParam$nugget)*diag(nrow(model$Phi))
  invPhiGammaPhit <- chol2inv(chol(PhiGammaPhit))
  pred$mu <- GammaPhit %*% invPhiGammaPhit %*% model$y
  pred$Sigma <- model$Gamma - GammaPhit %*% invPhiGammaPhit %*% t(GammaPhit)

  # computing the conditional mode vector given the interpolation points
  # and the inequality constraints
  # invGamma <- chol2inv(chol(model$Gamma))
  # d <- matrix(0, nrow = nrow(model$Gamma))
  # pred$xi.map <- solve.QP(invGamma, d, t(A), b, nrow(model$x))$solution
  if (min(eigen(pred$Sigma, symmetric = TRUE)$values) <= 0) # numerical stability
    pred$Sigma <- pred$Sigma + model$kernParam$nugget*diag(nrow(pred$Sigma))
  invSigma <- chol2inv(chol(pred$Sigma))
  pred$xi.map <- solve.QP(invSigma, t(pred$mu) %*% invSigma,
                          t(model$lineqSys$M), model$lineqSys$g)$solution
  pred$ymap <- pred$Phi.test %*% pred$xi.map
  if (return_model)
    pred$model <- model
  return(pred)
}

#' @title Simulation Method for the \code{"lineqGP"} S3 Class
#' @description Simulation method for the \code{"lineqGP"} S3 class.
#' 
#' @param object an object with class \code{"lineqGP"}
#' @param xtest a vector (or matrix) with the test input design
#' @param nsim	the number of simulations
#' @param seed see \code{\link{simulate}}
#' @param ... further arguments passed to or from other methods
#' 
#' @return A \code{"lineqGP"} object with the following elements
#' \item{x}{a vector (or matrix) with the training input design}
#' \item{y}{the training output vector at \code{x}}
#' \item{xtest}{a vector (or matrix) with the test input design}
#' \item{Phi.test}{a matrix corresponding to the hat basis functions evaluated
#' at \code{xtest}. The basis functions are indexed by rows}
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
#' @seealso \code{\link{create.lineqGP}}, \code{\link{augment.lineqGP}},
#'          \code{\link{predict.lineqGP}}
#'          
#' @author A. F. Lopez-Lopera
#'
#' @references A. F. Lopez-Lopera, F. Bachoc, N. Durrande and O. Roustant (2017),
#' "Finite-dimensional Gaussian approximation with linear inequality constraints".
#' \emph{SIAM/ASA Journal on Uncertainty Quantification}, 6(3): 1224–1255.
#' <doi:10.1137/17M1153157>
#'
#' @examples
#' # creating the model
#' sigfun <- function(x) return(1/(1+exp(-7*(x-0.5))))
#' x <- seq(0, 1, length = 5)
#' y <- sigfun(x)
#' model <- create(class = "lineqGP", x, y, constrType = "monotonicity", m = 30)
#'
#' # modifying the covariance parameters
#' model$kernParam$par <- c(1, 0.2)
#'
#' # sampling from the model
#' xtest <- seq(0, 1, length = 100)
#' ytest <- sigfun(xtest)
#' sim.model <- simulate(model, xtest, nsim = 50, seed = 1)
#' mu <- apply(sim.model$ysim, 1, mean)
#' qtls <- apply(sim.model$ysim, 1, quantile, probs = c(0.05, 0.95))
#' matplot(xtest, t(qtls), type = "l", lty = 1, col = "gray90",
#'         main = "Constrained Kriging model")
#' polygon(c(xtest, rev(xtest)), c(qtls[2,], rev(qtls[1,])), col = "gray90", border = NA)
#' lines(xtest, ytest, lty = 2)
#' lines(xtest, mu, type = "l", col = "darkgreen")
#' points(x, y, pch = 20)
#' legend("right", c("ytrain", "ytest", "mean", "confidence"), lty = c(NaN,2,1,NaN),
#'        pch = c(20,NaN,NaN,15), col = c("black", "black", "darkgreen", "gray80"))
#'
#' @importFrom stats simulate
#' @export
simulate.lineqGP <- function(object, nsim = 1, seed = NULL, xtest, ...) {
  model <- object
  predtime <- proc.time()
  pred <- predict(model, xtest, return_model = TRUE)
  predtime <- proc.time() - predtime
  model <- pred$model
  pred$model <- NULL
  
  # computing the transformed conditional mode and covariance matrix given
  # the interpolation points and the inequality constraints
  eta.mu <- as.vector(pred$Lambda %*% pred$mu)
  eta.map <- as.vector(pred$Lambda %*% pred$xi.map)
  Sigma.eta <- pred$Lambda %*% pred$Sigma %*% t(pred$Lambda)
  Sigma.eta <- Sigma.eta + model$kernParam$nugget*diag(nrow(Sigma.eta))

  # listing control terms
  control <- as.list(unlist(model$localParam$samplingParam))
  control$mvec <- sum(model$localParam$mvec) # for HMC
  control$constrType <- model$constrType # for HMC

  # sampling from the truncated multinormal
  tmvPar <- list(mu = eta.mu, map = eta.map, Sigma = Sigma.eta, lb = pred$lb, ub = pred$ub)
  class(tmvPar) <- model$localParam$sampler
  set.seed(seed)
  simtime <- proc.time()
  eta <- tmvrnorm(tmvPar, nsim, control)
  simtime <- proc.time() - simtime
  
  # passing some terms to the simulated model
  simModel <- list()
  simModel$x <- model$x
  simModel$y <- model$y
  simModel$xtest <- xtest
  simModel$Phi.test <- pred$Phi.test
  simModel$predtime <- predtime
  simModel$simtime <- simtime
  
  simModel$xi.map <- pred$xi.map
  simModel$xi.sim <- qr.solve(pred$Lambda, eta)
  simModel$ymap <- pred$ymap
  simModel$ysim <- pred$Phi.test %*% simModel$xi.sim
  class(simModel) <- class(model)
  return(simModel)
}
