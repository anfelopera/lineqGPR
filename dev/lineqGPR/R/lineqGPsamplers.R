#' @title  \code{"tmvrnorm"} Sampler for \code{"RSM"} (Rejection Sampling from the Mode) S3 Class
#' @description Sampler for truncated multivariate normal distributions
#' via RSM according to (Maatouk and Bay, 2017).
#' 
#' @param object an object with \code{"RSM"} S3 class containing:
#'        \code{mu} (mean vector), \code{Sigma} (covariance matrix)
#'        \code{lb} (lower bound vector), \code{ub} (upper bound vector)
#' @param nsim an integer corresponding to the number of simulations
#' @param control extra parameters required for the MC/MCMC sampler
#' @param ... further arguments passed to or from other methods
#' 
#' @return A matrix with the simulated samples. Samples are indexed by columns
#'
#' @seealso \code{\link{tmvrnorm.HMC}}, \code{\link{tmvrnorm.ExpT}}
#' 
#' @author A. F. Lopez-Lopera
#'
#' @references H. Maatouk and X. Bay (2017),
#' "Gaussian process emulators for computer experiments with inequality constraints".
#' \emph{Mathematical Geosciences},
#' 49(5):557-582.
#'
#' @examples
#' n <- 100
#' x <- seq(0, 1, length = n)
#' Sigma <- kernCompute(x1 = x, type = "gaussian", par = c(1,0.2))
#' tmgPar <- list(mu = rep(0,n), Sigma = Sigma + 1e-9*diag(n), lb = rep(-1,n), ub = rep(1,n))
#' class(tmgPar) <- "RSM"
#' y <- tmvrnorm(tmgPar, nsim = 10)
#' matplot(x, y, type = 'l', ylim = c(-1,1),
#'         main = "Constrained samples using RSM")
#' abline(h = c(-1,1), lty = 2)
#'
#' @importFrom MASS mvrnorm
#' @export
tmvrnorm.RSM <- function(object, nsim, control = NULL, ...) {
  # precomputing some terms according to Maatouk et al. [2016]
  tmvPar <- object
  if (!("map" %in% names(tmvPar)))
    tmvPar$map <- tmvPar$mu
  
  invSigma <- chol2inv(chol(tmvPar$Sigma))
  invSigmamu_star <-  invSigma %*% tmvPar$map

  # generating the samples
  xi <- matrix(0, nrow = length(tmvPar$map), ncol = nsim)
  for (i in seq(nsim)) {
    condition <- FALSE
    while (condition == FALSE) {
      xi_temp <- mvrnorm(n = 1, tmvPar$map, Sigma = tmvPar$Sigma)
      while(!all(xi_temp >= tmvPar$lb & xi_temp <= tmvPar$ub)) {
        xi_temp <- mvrnorm(n = 1, tmvPar$map, Sigma = tmvPar$Sigma)
      }
      t <- exp( t(tmvPar$map - xi_temp) %*% invSigmamu_star)
      u <- runif(1)
      if (u <= t) condition <- TRUE
    }
    xi[, i] <- xi_temp
  }
  return(xi)
}

#' @title  \code{"tmvrnorm"} Sampler for \code{"HMC"} (Hamiltonian Monte Carlo) S3 Class
#' @description Sampler for truncated multivariate normal distributions
#' via Hamiltonian Monte Carlo using the package \code{tmg}  (Pakman and Paninski, 2014).
#' 
#' @param object an object with \code{"HMC"} S3 class containing:
#'        \code{mu} (mean vector), \code{Sigma} (covariance matrix)
#'        \code{lb} (lower bound vector), \code{ub} (upper bound vector)
#' @param nsim an integer corresponding to the number of simulations
#' @param control extra parameters required for the MC/MCMC sampler
#' @param ... further arguments passed to or from other methods
#' 
#' @return A matrix with the simulated samples. Samples are indexed by columns
#'
#' @seealso \code{\link{tmvrnorm.RSM}}, \code{\link{tmvrnorm.ExpT}}
#' 
#' @author A. F. Lopez-Lopera
#'
#' @references A. Pakman and L. Paninski (2014),
#' "Exact Hamiltonian Monte Carlo for truncated multivariate Gaussians".
#' \emph{Journal of Computational and Graphical Statistics},
#' 23(2):518-542.
#'
#' @examples
#' n <- 100
#' x <- seq(0, 1, length = n)
#' Sigma <- kernCompute(x1 = x, type = "gaussian", par = c(1,0.2))
#' tmgPar <- list(mu = rep(0,n), Sigma = Sigma + 1e-9*diag(n), lb = rep(-1,n), ub = rep(1,n))
#' class(tmgPar) <- "HMC"
#' y <- tmvrnorm(tmgPar, nsim = 10)
#' matplot(x, y, type = 'l', ylim = c(-1,1),
#'         main = "Constrained samples using Hamiltonian MC")
#' abline(h = c(-1,1), lty = 2)
#'
# #' @importFrom tmg rtmg
# #' @importFrom Matrix bdiag
#' @export
tmvrnorm.HMC <- function(object, nsim, control = list(burn.in = 1e2), ...) {
  if (!("burn.in" %in% names(control)))
    control$burn.in <- 1e2
  if (!("mvec" %in% names(control)))
    control$mvec <- length(object$mu)
  if (!("constrType" %in% names(control)))
    control$constrType <- "boundedness"
  tmvPar <- object
  if (!("map" %in% names(tmvPar)))
    tmvPar$map <- tmvPar$mu

  # precomputing some terms
  M <- init_vec <- g <- numeric()
  # idxInit <- 1
  # idxEnd <- 0
  # epsilon <- 1e-9
  # for (i in seq(length(control$constrType))) {
  #   idxInit <- idxEnd + 1
  #   idxEnd <- idxEnd + control$mvec[i]
  #   lsys <- bounds2lineqSys(control$mvec[i], tmvPar$lb[idxInit:idxEnd],
  #                           tmvPar$ub[idxInit:idxEnd], A = diag(control$mvec[i]),
  #                           lineqSysType = 'oneside')
  #   M <- Matrix::bdiag(M, lsys$M)
  #   g <- c(g, lsys$g)
  #   init_vec <- c(init_vec, pmin(pmax(tmvPar$map[idxInit:idxEnd],
  #                                     tmvPar$lb[idxInit:idxEnd]+epsilon),
  #                                tmvPar$ub[idxInit:idxEnd]-epsilon))
  # }
  epsilon <- 1e-9
  lsys <- bounds2lineqSys(control$mvec, tmvPar$lb,
                          tmvPar$ub, A = diag(control$mvec),
                          lineqSysType = 'oneside')
  M <- Matrix::bdiag(M, lsys$M)
  g <- c(g, lsys$g)
  init_vec <- c(init_vec, pmin(pmax(tmvPar$map,
                                    tmvPar$lb+epsilon),
                               tmvPar$ub-epsilon))
  M <- M[,-1]

  # simulating samples using the package "tmg"
  H <- solve(tmvPar$Sigma)
  xi <- tmg::rtmg(nsim, M = H, r = as.vector(t(tmvPar$mu) %*% H), initial = init_vec,
                  f = as.matrix(M), g = as.vector(g), burn.in = control$burn.in)
  return(t(xi))
}

#' @title  \code{"tmvrnorm"} Sampler for \code{"ExpT"} (Exponential Tilting) S3 Class
#' @description Sampler for truncated multivariate normal distributions
#' via exponential tilting using the package \code{TruncatedNormal} (Botev, 2017).
#' 
#' @param object an object with \code{"ExpT"} S3 class containing:
#'        \code{mu} (mean vector), \code{Sigma} (covariance matrix)
#'        \code{lb} (lower bound vector), \code{ub} (upper bound vector)
#' @param nsim an integer corresponding to the number of simulations
#' @param control extra parameters required for the MC/MCMC sampler
#' @param ... further arguments passed to or from other methods
#' 
#' @return A matrix with the simulated samples. Samples are indexed by columns
#'
#' @seealso \code{\link{tmvrnorm.RSM}}, \code{\link{tmvrnorm.HMC}}
#' 
#' @author A. F. Lopez-Lopera
#'
#' @references Z. I. Botev (2017),
#' "The normal law under linear restrictions: simulation and estimation via minimax tilting".
#' \emph{Journal of the Royal Statistical Society: Series B (Statistical Methodology)},
#' 79(1):125-148.
#'
#' @examples
#' n <- 100
#' x <- seq(0, 1, length = n)
#' Sigma <- kernCompute(x1 = x, type = "gaussian", par = c(1,0.2))
#' tmgPar <- list(mu = rep(0,n), Sigma = Sigma + 1e-9*diag(n), lb = rep(-1,n), ub = rep(1,n))
#' class(tmgPar) <- "ExpT"
#' y <- tmvrnorm(tmgPar, nsim = 10)
#' matplot(x, y, type = 'l', ylim = c(-1,1),
#'         main = "Constrained samples using exponential tilting")
#' abline(h = c(-1,1), lty = 2)
#'
#' @export
tmvrnorm.ExpT <- function(object, nsim, control = NULL, ...) {
  tmvPar <- object
  # simulating samples using the package "TruncatedNormal" (Truncated Multivariate Normal)
  xi <- matrix(tmvPar$mu, nrow = nrow(tmvPar$Sigma), ncol = nsim) +
    TruncatedNormal::mvrandn(tmvPar$lb-tmvPar$mu, tmvPar$ub-tmvPar$mu, tmvPar$Sigma, nsim)
  return(xi)
}
