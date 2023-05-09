#' @title Gaussian Process Model Optimizations
#' @description Function for optimizations of \code{"lineqGP"} S3 class objects.
#' 
#' @param model a list with the structure of the constrained Kriging model.
#' @param x0 the initial values for the parameters to be optimized over.
#' @param eval_f a function to be minimized, with first argument the vector of parameters
#' over which minimization is to take place. It should return a scalar result.
#' @param lb a vector with lower bounds of the params. The params are forced to be positive.
#' See \code{\link{nloptr}}.
#' @param ub a vector with upper bounds of the params. See \code{\link{nloptr}}.
#' @param opts see \code{\link{nl.opts}}. Parameter \code{parfixed} indices of
#' fixed parameters to do not be optimised. If \code{estim.varnoise} is true, the
#' noise variance is estimated.
#' @param seed an optional number. Set a seed to replicate results.
#' @param estim.varnoise an optional logical. If \code{TRUE}, a noise variance is estimated.
#' @param bounds.varnoise a vector with bounds of noise variance.
#' @param add.constr an optional logical. If \code{TRUE}, the inequality constraints are taken
#' into account in the optimization.
#' @param additive an optional logical. If \code{TRUE}, the likelihood of an additive GP model
#' is computed in the optimization.
#' @param mcmc.opts if \code{add.constr}, mcmc options passed to methods.
#' @param max.trials the value of the maximum number of trials when errors are produced by instabilities.
#' @param ... further arguments passed to or from other methods.
#' 
#' @return An optimized \code{lineqGP} model.
#'
#' @section Comments:
#' This function has to be improved in the future for more stable procedures.
#' Cros-validation (CV) methods could be implemented in future versions.
#' 
#' @seealso \code{\link{nloptr}}
#' 
#' @author A. F. Lopez-Lopera
#'
#' @import nloptr
#' @importFrom  purrr map
#' @export
lineqGPOptim <- function(model,
                         x0 = model$kernParam$par,
                         eval_f = "logLik",
                         lb = rep(0.01, length(x0)),
                         ub = rep(Inf, length(x0)),
                         opts = list(algorithm = "NLOPT_LD_MMA",
                                     print_level = 0,
                                     ftol_abs = 1e-3,
                                     maxeval = 50,
                                     check_derivatives = FALSE,
                                     parfixed = rep(FALSE, length(x0))),
                         seed = 1,
                         estim.varnoise = FALSE,
                         bounds.varnoise = c(0, Inf),
                         add.constr = FALSE,
                         additive = FALSE,
                         mcmc.opts = list(probe = "Genz", nb.mcmc = 1e3),
                         max.trials = 10, ...) {
  if (additive) {
    xmodel <- unlist(purrr::map(model$kernParam, "par"))
    if (length(x0) != length(xmodel)) {
      x0 <- xmodel
      lb <- rep(0.01, length(x0))
      ub <- rep(Inf, length(x0))
      opts$parfixed <-  rep(FALSE, length(x0))
    }
  }
  
  model <- augment(model)
  if (!("parfixed" %in% names(opts)))
    opts$parfixed <- rep(FALSE, length(par))
  # if (!("bounds.varnoise" %in% names(opts)))
  #   opts$estim.varnoise = FALSE
  
  if (!("algorithm" %in% names(opts)))
    opts$algorithm <- "NLOPT_LD_MMA"
  
  if (!("print_level" %in% names(opts)))
    opts$print_level <- 0
  
  if (!("ftol_abs" %in% names(opts)))
    opts$ftol_abs <- 1e-3
  
  if (!("maxeval" %in% names(opts)))
    opts$maxeval <- 50
  
  if (!("check_derivatives" %in% names(opts)))
    opts$check_derivatives <- FALSE
  
  # changing the functions according to fn
  if (eval_f == "logLik") {
    if (additive) 
      eval_f = paste(eval_f, "Additive", sep = "")
    fn <- paste(eval_f, "Fun", sep = "")
    gr <- paste(eval_f, "Grad", sep = "")
    if (add.constr) {
      fn <- paste("constr", fn, sep = "")
      gr <- paste("constr", gr, sep = "")
    }
  }   # CV methods will be implemented in future versions

  if (estim.varnoise) {
    if (!("varnoise" %in% names(model)))
      model$varnoise <- 0
    x0 <- c(x0, model$varnoise)
    lb <- c(lb, bounds.varnoise[1])
    ub <- c(ub, bounds.varnoise[2])
  }
  # seed <- 5e1^sum(x0) # to try with different and fixed initial parameters
  trial <- 1
  while(trial <= max.trials) {
    optim <- try(nloptr(x0 = x0,
                        eval_f = get(fn),
                        eval_grad_f = get(gr),
                        lb = lb,
                        ub = ub,
                        opts = opts,
                        model = model,
                        mcmc.opts = mcmc.opts,
                        parfixed = opts$parfixed,
                        estim.varnoise = estim.varnoise, ...))
    if (class(optim) == "try-error") {
      set.seed(seed)
      x0Rand <- runif(length(x0), lb, ub)
      x0[!opts$parfixed] <- x0Rand[!opts$parfixed]
      message("Non feasible solution. Re-initialising the optimizer.")
      seed <- seed + 1
      trial <- trial + 1
    } else {
      break
    }
  }
  if (trial > max.trials) {
    stop('Max number of trials has been exceeded with errors.')
  } else {
    # expanding the model according to the optimized parameters
    if (estim.varnoise) {
      model$varnoise <- optim$solution[length(optim$solution)]
      optim$solution <- optim$solution[-length(optim$solution)]
    }
    if (additive) {
      for (k in 1:model$d)
        model$kernParam[[k]]$par <- optim$solution[(k-1)*length(model$kernParam[[k]]$par) + 1:2]
    } else {
      model$kernParam$par <- optim$solution
    }
    model$kernParam$l <- optim$objective # used for multistart procedures
    model <- augment(model)
    return(model)
  }
}

#' @title Log-Likelihood of a Gaussian Process.
#' @description Compute the negative log-likelihood of a Gaussian Process.
#' 
#' @param par the values of the covariance parameters.
#' @param model an object with \code{"lineqGP"} S3 class.
#' @param parfixed not used.
#' @param mcmc.opts not used.
#' @param estim.varnoise If \code{true}, a noise variance is estimated.
#' 
#' @return The value of the negative log-likelihood.
#' 
#' @seealso \code{\link{logLikGrad}}, \code{\link{constrlogLikFun}},
#'          \code{\link{constrlogLikGrad}}
#'          
#' @author A. F. Lopez-Lopera
#'
#' @references Rasmussen, C. E. and Williams, C. K. I. (2005),
#' "Gaussian Processes for Machine Learning (Adaptive Computation and Machine Learning)".
#' \emph{The MIT Press}.
#'
#' @export
logLikFun <- function(par = model$kernParam$par, model,
                      parfixed = NULL, mcmc.opts = NULL,
                      estim.varnoise = FALSE) {
  switch (class(model),
    lineqGP = {
      m <- model$localParam$m
      u <- vector("list", model$d) # list with d sets of knots
      for (j in 1:model$d)
        u[[j]] <- matrix(seq(0, 1, by = 1/(m[j]-1)))
    }, lineqMaxModGP = {
      u <- model$localParam$u
      m <- c(unlist(model$localParam$m))
    } 
  )
  
  if (estim.varnoise) {
    varnoise <- par[length(par)]
    par <- par[-length(par)]
  } else {
    varnoise <- model$varnoise
  }

  AllGammas <- vector("list", model$d) # list with the d covariance matrices
  if (length(par) == 2 && length(unique(m)) == 1) {
    if (is.list(u)) {
      uBase <- u[[1]]
    } else {
      uBase <- u
    }
    
    GammaBase <- kernCompute(uBase, uBase,
                             model$kernParam$type,
                             par)
    for (j in 1:model$d) {
      AllGammas[[j]] <- GammaBase
    }
  } else {
    if (length(par[-1]) != model$d) {
      par <- c(par[1], rep(par[2], model$d))
      names(par) <- c(names(par)[1],
                      paste(names(par)[2], seq(model$d), sep = ""))
    }
    for (j in 1:model$d) {
      AllGammas[[j]] <- kernCompute(u[[j]], u[[j]],
                                    model$kernParam$type,
                                    par[c(1, j+1)])
    }
  }
  expr <- paste("AllGammas[[", seq(model$d), "]]", collapse = " %x% ")
  Gamma <- eval(parse(text = expr))

  Phi <- model$Phi
  Kyy <- Phi %*% Gamma %*% t(Phi)
  # if (estim.varnoise)
  Kyy <- Kyy + varnoise*diag(nrow(Kyy))
  cholKyy <- chol(Kyy + model$kernParam$nugget*diag(nrow(Kyy)))
  logDetKyy <- 2*sum(diag(log(cholKyy)))
  invKyy <- chol2inv(cholKyy)
  f <- 0.5*(logDetKyy + nrow(Kyy)*log(2*pi) + t(model$y)%*%invKyy%*%model$y)
  return(f)
}

#' @title Gradient of the Log-Likelihood of a Gaussian Process.
#' @description Compute the gradient of the negative log-likelihood of a Gaussian Process.
#' 
#' @param par the values of the covariance parameters.
#' @param model an object with \code{"lineqGP"} S3 class.
#' @param parfixed indices of fixed parameters to do not be optimised.
#' @param mcmc.opts not used.
#' @param estim.varnoise If \code{true}, a noise variance is estimated.
#' 
#' @return the gradient of the negative log-likelihood.
#'
#' @seealso \code{\link{logLikFun}}, \code{\link{constrlogLikFun}},
#'          \code{\link{constrlogLikGrad}}
#'          
#' @author A. F. Lopez-Lopera
#'
#' @references Rasmussen, C. E. and Williams, C. K. I. (2005),
#' "Gaussian Processes for Machine Learning (Adaptive Computation and Machine Learning)".
#' \emph{The MIT Press}.
#'
#' @export
logLikGrad <- function(par = model$kernParam$par, model,
                       parfixed = rep(FALSE, length(par)), mcmc.opts = NULL,
                       estim.varnoise = FALSE) {
  switch (class(model),
          lineqGP = {
            m <- model$localParam$m
            u <- vector("list", model$d) # list with d sets of knots
            for (j in 1:model$d)
              u[[j]] <- matrix(seq(0, 1, by = 1/(m[j]-1)))
          }, lineqMaxModGP = {
            u <- model$localParam$u
            m <- c(unlist(model$localParam$m))
          } 
  )

  if (estim.varnoise) {
    varnoise <- par[length(par)]
    par <- par[-length(par)]
  }
  
  AllGammas <- vector("list", model$d) # list with the d covariance matrices
  if (length(par) == 2 && length(unique(m)) == 1) {
    if (is.list(u)) {
      uBase <- u[[1]]
    } else {
      uBase <- u
    }
    GammaBase <- kernCompute(uBase, uBase,
                             model$kernParam$type,
                             par)
    for (j in 1:model$d) {
      AllGammas[[j]] <- GammaBase
    }
  } else {
    if (length(par[-1]) != model$d) {
      par <- c(par[1], rep(par[2], model$d))
      names(par) <- c(names(par)[1],
                      paste(names(par)[2], seq(model$d), sep = ""))
    }
    for (j in 1:model$d) {
      AllGammas[[j]] <- kernCompute(u[[j]], u[[j]],
                                    model$kernParam$type,
                                    par[c(1, j+1)])
    }
  }
  expr <- paste("AllGammas[[", seq(model$d), "]]", collapse = " %x% ")
  Gamma <- eval(parse(text = expr))

  Phi <- model$Phi
  Kyy <- Phi %*% Gamma %*% t(Phi)
  if (estim.varnoise)
    Kyy <- Kyy + varnoise*diag(nrow(Kyy))
  cholKyy <- chol(Kyy + model$kernParam$nugget*diag(nrow(Kyy)))
  invKyy <- chol2inv(cholKyy)

  gradf <- rep(0, length(par))
  idx_iter <- seq(length(par))
  idx_iter <- idx_iter[parfixed == FALSE]

  expr <- matrix(paste("AllGammas[[", seq(model$d),"]]", sep = ""),
                 model$d, model$d, byrow = TRUE)
  diag(expr) <- paste("attr(AllGammas[[", seq(model$d),
                      "]], 'gradient')[[2]]", sep = "")
  for (i in idx_iter) {
    if (i == 1) {
      gradGammaTemp <- model$d*Gamma/par[1]
    } else {
      expr2 <- paste(expr[i-1, ], collapse = " %x% ")
      gradGammaTemp <- eval(parse(text = expr2))
    }
    gradKyyTemp <- Phi %*% gradGammaTemp %*% t(Phi)
    alpha <- invKyy %*% model$y
    gradf[i] <- 0.5*sum(diag((invKyy - alpha%*%t(alpha)) %*% gradKyyTemp))
  }
  if (estim.varnoise) {
    gradKyynoise <- diag(nrow(Kyy))
    alpha <- invKyy %*% model$y
    gradf <- c(gradf, 
               0.5*sum(diag((invKyy - alpha%*%t(alpha)) %*% gradKyynoise)))
  }
  
  # gradf <- nl.grad(par, logLikFun, heps = 1e-9, model, parfixed, mcmc.opts)
  # gradf[parfixed] <- 0
  return(gradf)
}

#' @title Log-Constrained-Likelihood of a Gaussian Process.
#' @description Compute the negative log-constrained-likelihood of a Gaussian Process
#' conditionally to the inequality constraints (Lopez-Lopera et al., 2019).
#' 
#' @param par the values of the covariance parameters.
#' @param model an object with \code{"lineqGP"} S3 class.
#' @param parfixed not used.
#' @param mcmc.opts mcmc options. \code{mcmc.opts$probe} A character string corresponding
#' to the estimator for the orthant multinormal probabilities.
#' Options: \code{"Genz"} (Genz, 1992), \code{"ExpT"} (Botev, 2017).
#' If \code{probe == "ExpT"}, \code{mcmc.opts$nb.mcmc}  is the number of MCMC
#' samples used for the estimation.
#' @param estim.varnoise If \code{true}, a noise variance is estimated.
#' 
#' @return The value of the negative log-constrained-likelihood.
#' 
#' @details Orthant multinormal probabilities are estimated according to
#' (Genz, 1992; Botev, 2017). See (Lopez-Lopera et al., 2017).
#'
#' @seealso \code{\link{constrlogLikGrad}}, \code{\link{logLikFun}},
#'          \code{\link{logLikGrad}}
#'          
#' @author A. F. Lopez-Lopera
#'
#' @references A. F. Lopez-Lopera, F. Bachoc, N. Durrande and O. Roustant (2017),
#' "Finite-dimensional Gaussian approximation with linear inequality constraints".
#' \emph{SIAM/ASA Journal on Uncertainty Quantification}, 6(3): 1224–1255.
#' <doi:10.1137/17M1153157>
#'
#' @references F. Bachoc, A. Lagnoux and A. F. Lopez-Lopera (2019),
#' "Maximum likelihood estimation for Gaussian processes under inequality constraints".
#' \emph{Electronic Journal of Statistics}, 13 (2): 2921-2969.
#' <doi:10.1214/19-EJS1587>
#'
#' A. Genz (1992),
#' "Numerical computation of multivariate normal probabilities".
#' \emph{Journal of Computational and Graphical Statistics},
#' 1:141-150.
#'
#' Z. I. Botev (2017),
#' "The normal law under linear restrictions: simulation and estimation via minimax tilting".
#' \emph{Journal of the Royal Statistical Society: Series B (Statistical Methodology)},
#' 79(1):125-148.
#'
# #' @importFrom TruncatedNormal mvrandn mvNcdf
# #' @import mvtnorm
#' @export
constrlogLikFun <- function(par = model$kernParam$par, model, parfixed = NULL,
                            mcmc.opts = list(probe = c("Genz"), nb.mcmc = 1e3),
                            estim.varnoise = FALSE) {
  # u <- matrix(seq(0, 1, by = 1/(model$localParam$m-1)), ncol = 1) # discretization vector
  # Gamma <- kernCompute(u, u, model$kernParam$type, par = par)
  # if (model$d == 2)
    # Gamma <- kronecker(Gamma, Gamma)
  m <- model$localParam$m
  u <- vector("list", model$d) # list with d sets of knots
  for (j in 1:model$d)
    u[[j]] <- matrix(seq(0, 1, by = 1/(m[j]-1)))

  if (estim.varnoise) {
    varnoise <- par[length(par)]
    par <- par[-length(par)]
  }
  
  AllGammas <- vector("list", model$d) # list with the d covariance matrices
  if (length(par) == 2 && length(unique(m)) == 1) {
    uBase <- u[[1]]
    GammaBase <- kernCompute(uBase, uBase,
                             model$kernParam$type,
                             par)
    for (j in 1:model$d) {
      AllGammas[[j]] <- GammaBase
    }
  } else {
    if (length(par[-1]) != model$d) {
      par <- c(par[1], rep(par[2], model$d))
      names(par) <- c(names(par)[1],
                      paste(names(par)[2], seq(model$d), sep = ""))
    }
    for (j in 1:model$d) {
      AllGammas[[j]] <- kernCompute(u[[j]], u[[j]],
                                    model$kernParam$type,
                                    par[c(1, j+1)])
    }
  }
  expr <- paste("AllGammas[[", seq(model$d), "]]", collapse = " %x% ")
  Gamma <- eval(parse(text = expr))
  Phi <- model$Phi
  GammaPhit <- Gamma %*% t(Phi)

  # computing the negative of the unconstrained log-likelihood
  Kyy <- Phi %*% GammaPhit
  if (estim.varnoise)
    Kyy <- Kyy + varnoise*diag(nrow(Kyy))
  Kyy <- Kyy + model$kernParam$nugget*diag(nrow(Kyy))
  cholKyy <- chol(Kyy)
  logDetKyy <- 2*sum(diag(log(cholKyy)))
  invKyy <- chol2inv(cholKyy)
  f <- 0.5*(logDetKyy + nrow(Kyy)*log(2*pi) + t(model$y)%*%invKyy%*%model$y)

  # computing the other terms of the constrained log-likelihood
  GammaPhitInvKyy <- GammaPhit %*% invKyy
  mu <- GammaPhitInvKyy %*% model$y
  mu.eta <- model$Lambda %*% mu
  Sigma <- Gamma - GammaPhitInvKyy %*% t(GammaPhit)
  Sigma.eta <- model$Lambda %*% Sigma %*% t(model$Lambda)
  if (min(eigen(Sigma.eta, symmetric = TRUE)$values) <= 0) # numerical stability
    Sigma.eta <- Sigma.eta + 1e-5*diag(nrow(Sigma.eta))
  Gamma.eta <- model$Lambda %*% Gamma %*% t(model$Lambda)
  if (min(eigen(Gamma.eta, symmetric = TRUE)$values) <= 0) # numerical stability
    Gamma.eta <- Gamma.eta + 1e-5*diag(nrow(Gamma.eta))

  set.seed(7)
  switch(mcmc.opts$probe,
         Genz = {
           f <- f - log(mvtnorm::pmvnorm(model$lb, model$ub,
                                         as.vector(mu.eta), sigma = Sigma.eta)[[1]]) +
             log(mvtnorm::pmvnorm(model$lb, model$ub,
                                  rep(0, nrow(Gamma.eta)), sigma = Gamma.eta)[[1]])
         }, ExpT = {
           f <- f - log(TruncatedNormal::mvNcdf(model$lb-mu.eta, model$ub-mu.eta,
                               Sigma.eta, mcmc.opts$nb.mcmc)$prob) +
             log(TruncatedNormal::mvNcdf(model$lb, model$ub, Gamma.eta, mcmc.opts$nb.mcmc)$prob)
         })
  return(f)
}

#' @title  Numerical Gradient of the Log-Constrained-Likelihood of a Gaussian Process.
#' @description Compute the gradient numerically of the negative log-constrained-likelihood of a Gaussian Process
#' conditionally to the inequality constraints (Lopez-Lopera et al., 2019).
#' 
#' @param par the values of the covariance parameters.
#' @param model an object with class \code{lineqGP}.
#' @param parfixed indices of fixed parameters to do not be optimised.
#' @param mcmc.opts mcmc options. \code{mcmc.opts$probe} A character string corresponding
#' to the estimator for the orthant multinormal probabilities.
#' Options: \code{"Genz"} (Genz, 1992), \code{"ExpT"} (Botev, 2017).
#' If \code{probe == "ExpT"}, \code{mcmc.opts$nb.mcmc}  is the number of MCMC
#' samples used for the estimation.
#' @param estim.varnoise If \code{true}, a noise variance is estimated.
#' 
#' @return The gradient of the negative log-constrained-likelihood.
#' 
#' @details Orthant multinormal probabilities are estimated via (Genz, 1992; Botev, 2017).
#'
#' @seealso \code{\link{constrlogLikFun}}, \code{\link{logLikFun}},
#'          \code{\link{logLikGrad}}
#' 
#' @section Comments:
#' As orthant multinormal probabilities don't have explicit expressions,
#' the gradient is implemented numerically based on \code{\link{nl.grad}}.
#' 
#' @author A. F. Lopez-Lopera
#'
#' @references A. F. Lopez-Lopera, F. Bachoc, N. Durrande and O. Roustant (2017),
#' "Finite-dimensional Gaussian approximation with linear inequality constraints".
#' \emph{SIAM/ASA Journal on Uncertainty Quantification}, 6(3): 1224–1255.
#' <doi:10.1137/17M1153157>
#'
#' @references F. Bachoc, A. Lagnoux and A. F. Lopez-Lopera (2019),
#' "Maximum likelihood estimation for Gaussian processes under inequality constraints".
#' \emph{Electronic Journal of Statistics}, 13 (2): 2921-2969.
#' <doi:10.1214/19-EJS1587>
#'
#' A. Genz (1992),
#' "Numerical computation of multivariate normal probabilities".
#' \emph{Journal of Computational and Graphical Statistics},
#' 1:141-150.
#'
#' Z. I. Botev (2017),
#' "The normal law under linear restrictions: simulation and estimation via minimax tilting".
#' \emph{Journal of the Royal Statistical Society: Series B (Statistical Methodology)},
#' 79(1):125-148.
#'
#' @export
constrlogLikGrad <- function(par = model$kernParam$par, model,
                             parfixed = rep(FALSE, length(par)),
                             mcmc.opts = list(probe = "Genz", nb.mcmc = 1e3),
                             estim.varnoise = FALSE) {
  # gradf <- nl.grad(par, constrlogLikFun, heps = 1e-9, model, parfixed, mcmc.opts)
  # gradf[parfixed] <- 0

  if (estim.varnoise)
    parfixed <- c(parfixed, FALSE)
  
  # this is a faster version of "nl.grad" when some parameters are fixed.
  if (!is.numeric(par))
    stop("Argument 'par' must be a numeric value.")
  fun <- match.fun(constrlogLikFun)
  fn <- function(x) fun(x, model, parfixed, mcmc.opts, estim.varnoise)
  if (length(fn(par)) != 1)
    stop("Function 'f' must be a univariate function of 2 variables.")
  n <- length(par)
  hh <- rep(0, n)
  gradf <- rep(0, n)

  # the gradient is computed only w.r.t. the parameters to be optimised
  idx_iter <- seq(n)
  idx_iter <- idx_iter[parfixed == FALSE]
  heps <- 1e-9
  for (i in idx_iter) {
    hh[i] <- heps
    gradf[i] <- (fn(par + hh) - fn(par - hh))/(2 * heps)
    hh[i] <- 0
  }
  return(gradf)
}

#' @title Log-Likelihood of a Additive Gaussian Process.
#' @description Compute the negative log-likelihood of an Additive Gaussian Process.
#' 
#' @param par the values of the covariance parameters.
#' @param model an object with \code{"lineqAGP"} S3 class.
#' @param parfixed not used.
#' @param mcmc.opts not used.
#' @param estim.varnoise If \code{true}, a noise variance is estimated.
#' 
#' @return The value of the negative log-likelihood.
#' 
#' @seealso \code{\link{logLikAdditiveGrad}}
# #'          \code{\link{constrlogLikAdditiveFun}},
# #'          \code{\link{constrlogLikAdditiveGrad}}

#' @author A. F. Lopez-Lopera
#' 
#' @importFrom Matrix bdiag
#' 
#' @export
logLikAdditiveFun <- function(par = unlist(purrr::map(model$kernParam, "par")),
                              model,
                              parfixed = NULL, mcmc.opts = NULL,
                              estim.varnoise = FALSE) {
  m <- model$localParam$m
  mt <- sum(m)
  nt <- length(model$y) 
  ngroups <- model$localParam$ngroups
  
  u <- Gamma <- vector("list", ngroups)
  Phi <- vector("list", ngroups) # to be computed once and called it
  logDetGamma <- rep(0, ngroups)
  invGamma <- vector("list", ngroups)
  
  for (j in 1:ngroups)
    u[[j]] <- matrix(seq(0, 1, by = 1/(m[j]-1)), ncol = 1) # discretization vector
  
  if (estim.varnoise) {
    varnoise <- par[length(par)]
    par <- par[-length(par)]
  } else {
    varnoise <- model$varnoise
  }
  
  # computing the kernel matrix for the prior
  for (k in 1:ngroups) {
    Gamma[[k]] <- kernCompute(u[[k]], u[[k]], model$kernParam[[k]]$type,
                              par[(k-1)*length(model$kernParam[[k]]$par) + 1:2])
    Phi[[k]] <- basisCompute.lineqGP(model$x[, k], u[[k]])
  }
  
  if (mt < nt) {
    hfunBigPhi <- parse(text = paste("cbind(",
                                     paste("Phi[[", 1:ngroups, "]]", sep = "", collapse = ","),
                                     ")", sep = ""))
    bigPhi <- eval(hfunBigPhi)
    cholGamma <- lapply(Gamma, function(x) t(chol(x)))
    hfunBigCholGamma <- parse(text = paste("bdiag(",
                                           paste("cholGamma[[", 1:ngroups, "]]", sep = "", collapse = ","),
                                           ")", sep = ""))
    bigCholGamma <- as.matrix(eval(hfunBigCholGamma))
    PhibigCholGamma <- bigPhi %*% bigCholGamma
    ILtPhitPhiL <- varnoise*diag(mt) + t(PhibigCholGamma) %*% PhibigCholGamma
    cholILtPhitPhiL <- t(chol(ILtPhitPhiL))
    Lschur <- forwardsolve(cholILtPhitPhiL, t(PhibigCholGamma))
    invKyy <- (diag(nt) - t(Lschur)%*% Lschur)/varnoise
    logDetKyy <- (nt-mt)*log(varnoise) + 2*sum(log(diag(cholILtPhitPhiL)))
    
    f <- 0.5*(logDetKyy + nrow(invKyy)*log(2*pi) + t(model$y)%*%invKyy%*%model$y)
  } else {
    # ptm1 <- proc.time()
    Kyy <- matrix(0, nt, nt)
    for (k in 1:ngroups)
      Kyy <- Kyy + Phi[[k]] %*% Gamma[[k]] %*% t(Phi[[k]])
    if (estim.varnoise)
      Kyy <- Kyy + varnoise*diag(nt)
    # ptm1 <- proc.time() - ptm1

    ## Alternative computation using block matrices
    # ptm2 <- proc.time()
    # hfunBigPhi <- parse(text = paste("cbind(",
    #                                  paste("Phi[[", 1:ngroups, "]]", sep = "", collapse = ","),
    #                                  ")", sep = ""))
    # bigPhi <- eval(hfunBigPhi)
    # hfunBigGamma <- parse(text = paste("bdiag(",
    #                                    paste("Gamma[[", 1:ngroups, "]]", sep = "", collapse = ","),
    #                                    ")", sep = ""))
    # bigGamma <- eval(hfunBigGamma)
    # Kyy <- bigPhi %*% bigGamma %*% t(bigPhi) + model$varnoise * diag(nt)
    # Kyy <- matrix(Kyy, nt, nt)
    # ptm2 <- proc.time() - ptm2
    # print(ptm2 - ptm1)

    cholKyy <- t(chol(Kyy + model$nugget*diag(nt)))
    alpha <- forwardsolve(cholKyy, model$y)
    logDetKyy <- 2*sum(log(diag(cholKyy)))
    
    f <- 0.5*(logDetKyy + nrow(cholKyy)*log(2*pi) + t(alpha) %*% alpha)
  }
  return(f)
}

#' @title  Gradient of the Log-Likelihood of a Additive Gaussian Process.
#' @description Compute the gradient of the negative log-likelihood of an Additive Gaussian Process.
#' 
#' @param par the values of the covariance parameters.
#' @param model an object with \code{"lineqAGP"} S3 class.
#' @param parfixed indices of fixed parameters to do not be optimised.
#' @param mcmc.opts not used.
#' @param estim.varnoise If \code{true}, a noise variance is estimated.
#' 
#' @return the gradient of the negative log-likelihood.
#'
#' @seealso \code{\link{logLikAdditiveFun}}
#' 
#' @author A. F. Lopez-Lopera
#'
#' @export
logLikAdditiveGrad <- function(par = unlist(purrr::map(model$kernParam, "par")), 
                               model,  parfixed = rep(FALSE, model$d*length(par)),
                               mcmc.opts = NULL,
                               estim.varnoise = FALSE) {
  m <- model$localParam$m
  mt <- sum(m)
  nt <- length(model$y) 
  ngroups <- model$localParam$ngroups
  
  u <- Gamma <- vector("list", ngroups)
  Phi <- vector("list", ngroups)

  for (j in 1:ngroups)
    u[[j]] <- matrix(seq(0, 1, by = 1/(m[j]-1)), ncol = 1) # discretization vector
  
  if (estim.varnoise) {
    varnoise <- par[length(par)]
    par <- par[-length(par)]
  } else {
    varnoise <- model$varnoise
  }
  
  # computing the kernel matrix for the prior
  for (k in 1:ngroups) {
    Gamma[[k]] <- kernCompute(u[[k]], u[[k]], model$kernParam[[k]]$type,
                              par[(k-1)*length(model$kernParam[[k]]$par) + 1:2])
    Phi[[k]] <- basisCompute.lineqGP(model$x[, k], u[[k]])
  }
  
  if (mt < nt) {
    hfunBigPhi <- parse(text = paste("cbind(",
                                     paste("Phi[[", 1:ngroups, "]]", sep = "", collapse = ","),
                                     ")", sep = ""))
    bigPhi <- eval(hfunBigPhi)
    cholGamma <- lapply(Gamma, function(x) t(chol(x)))
    hfunBigCholGamma <- parse(text = paste("bdiag(",
                                           paste("cholGamma[[", 1:ngroups, "]]", sep = "", collapse = ","),
                                           ")", sep = ""))
    bigCholGamma <- as.matrix(eval(hfunBigCholGamma))
    PhibigCholGamma <- bigPhi %*% bigCholGamma
    ILtPhitPhiL <- varnoise*diag(mt) + t(PhibigCholGamma) %*% PhibigCholGamma
    cholILtPhitPhiL <- t(chol(ILtPhitPhiL))
    Lschur <- forwardsolve(cholILtPhitPhiL, t(PhibigCholGamma))
    invKyy <- (diag(nt) - t(Lschur)%*% Lschur)/varnoise
  } else {
    Kyy <- matrix(0, nt, nt)
    for (k in 1:ngroups)
      Kyy <- Kyy + Phi[[k]] %*% Gamma[[k]] %*% t(Phi[[k]])
    if (estim.varnoise)
      Kyy <- Kyy + varnoise*diag(nt)

    ## Alternative computation using block matrices
    # hfunBigPhi <- parse(text = paste("cbind(",
    #                                  paste("Phi[[", 1:ngroups, "]]", sep = "", collapse = ","),
    #                                  ")", sep = ""))
    # bigPhi <- eval(hfunBigPhi)
    # hfunBigGamma <- parse(text = paste("bdiag(",
    #                                    paste("Gamma[[", 1:ngroups, "]]", sep = "", collapse = ","),
    #                                    ")", sep = ""))
    # bigGamma <- eval(hfunBigGamma)
    # Kyy <- bigPhi %*% bigGamma %*% t(bigPhi) + model$varnoise * diag(nt)
    # Kyy <- matrix(Kyy, nt, nt)

    invKyy <- chol2inv(chol(Kyy + model$nugget*diag(nt)))
 
  }

  gradKyyTemp <- c()
  idx_iter <- seq(length(par))
  idx_iter <- idx_iter[parfixed == FALSE]

  alpha <- invKyy %*% model$y
  cteTermLik <- invKyy - alpha%*%t(alpha)
  gradf <- rep(0, length(par))
  for (k in 1:ngroups) {
    gradGammaSigma2 <- attr(Gamma[[k]], 'gradient')[[1]]
    gradKyySigma2 <- Phi[[k]] %*% gradGammaSigma2 %*% t(Phi[[k]])
    gradf[(k-1)*length(model$kernParam[[k]]$par) + 1] <- 0.5*sum(diag(cteTermLik %*% gradKyySigma2))
    
    gradGammaTheta <- attr(Gamma[[k]], 'gradient')[[2]]
    gradKyyTheta <- Phi[[k]] %*% gradGammaTheta %*% t(Phi[[k]])
    gradf[(k-1)*length(model$kernParam[[k]]$par) + 2] <- 0.5*sum(diag(cteTermLik %*% gradKyyTheta))
  }
  gradf[parfixed == TRUE] <- 0
  
  if (estim.varnoise) {
    gradKyynoise <- diag(nrow(invKyy))
    gradf <- c(gradf, 0.5*sum(diag(cteTermLik %*% gradKyynoise)))
  }
  
  # gradf <- nl.grad(par, logLikFun, heps = 1e-9, model, parfixed, mcmc.opts)
  # gradf[parfixed] <- 0
  return(gradf)
}

# #' @title Log-Constrained-Likelihood of an Additive Gaussian Process.
# #' @description Compute the negative log-constrained-likelihood of an Additive
# #' Gaussian Process conditionally to the inequality constraints
# #' (Lopez-Lopera et al., 2018).
# #' @param par the values of the covariance parameters.
# #' @param model an object with \code{"lineqAGP"} S3 class.
# #' @param parfixed not used.
# #' @param mcmc.opts mcmc options. \code{mcmc.opts$probe} A character string corresponding
# #' to the estimator for the orthant multinormal probabilities.
# #' Options: \code{"Genz"} (Genz, 1992), \code{"ExpT"} (Botev, 2017).
# #' If \code{probe == "ExpT"}, \code{mcmc.opts$nb.mcmc}  is the number of MCMC
# #' samples used for the estimation.
# #' @param estim.varnoise If \code{true}, a noise variance is estimated.
# #' @return The value of the negative log-constrained-likelihood.
# #' @details Orthant multinormal probabilities are estimated according to
# #' (Genz, 1992; Botev, 2017). See (Lopez-Lopera et al., 2018).
# #'
# #' @seealso \code{\link{constrlogLikAdditiveGrad}}, \code{\link{logLikAdditiveFun}},
# #'          \code{\link{logLikAdditiveGrad}}
# #' @author A. F. Lopez-Lopera
# #'
# #' @references Lopez-Lopera, A. F., Bachoc, F., Durrande, N., and Roustant, O. (2018),
# #' "Finite-dimensional Gaussian approximation with linear inequality constraints".
# #' \emph{SIAM/ASA Journal on Uncertainty Quantification}, 6(3): 1224-1255.
# #' <doi:10.1137/17M1153157>
# #'
# #' @references F. Bachoc, A. Lagnoux and A. F. Lopez-Lopera (2019),
# #' "Maximum likelihood estimation for Gaussian processes under inequality constraints".
# #' \emph{Electronic Journal of Statistics}, 13 (2): 2921-2969.
# #' <doi:10.1214/19-EJS1587>
# #'
# #' Genz, A. (1992),
# #' "Numerical computation of multivariate normal probabilities".
# #' \emph{Journal of Computational and Graphical Statistics},
# #' 1:141-150.
# #'
# #' Botev, Z. I. (2017),
# #' "The normal law under linear restrictions: simulation and estimation via minimax tilting".
# #' \emph{Journal of the Royal Statistical Society: Series B (Statistical Methodology)},
# #' 79(1):125-148.
# #'
# #' @import TruncatedNormal mvtnorm
# #' @export
# constrlogLikAdditiveFun <- function(par = unlist(purrr::map(model$kernParam, "par")),
#                                     model,
#                                     parfixed = NULL, mcmc.opts = NULL,
#                                     estim.varnoise = FALSE) {
#   m <- model$localParam$m
#   u <- Gamma <- vector("list", model$d)
#   Phi <- vector("list", model$d)
#   
#   for (j in 1:model$d)
#     u[[j]] <- matrix(seq(0, 1, by = 1/(m[j]-1)), ncol = 1) # discretization vector
#   
#   if (estim.varnoise) {
#     varnoise <- par[length(par)]
#     par <- par[-length(par)]
#   }
#   
#   # computing the kernel matrix for the prior
#   for (k in 1:model$d) {
#     Gamma[[k]] <- kernCompute(u[[k]], u[[k]], model$kernParam[[k]]$type,
#                               par[(k-1)*model$d + 1:2])
#     Phi[[k]] <- basisCompute.lineqGP(model$x[, k], u[[k]])
#   }
#   
#   GammaFull <- Gamma[[1]]
#   GammaPhitList <- vector("list", model$d)
#   GammaPhitList[[1]] <- Gamma[[1]] %*% t(Phi[[1]])
#   GammaPhitFull <- GammaPhitList[[1]]
#   PhiGammaPhitFull <- Phi[[1]] %*% GammaPhitFull
#   for (k in 2:model$d) {
#     GammaPhitList[[k]] <- Gamma[[k]] %*% t(Phi[[k]])
#     GammaPhitFull <- GammaPhitFull + GammaPhitList[[k]]
#     PhiGammaPhitFull <- PhiGammaPhitFull + Phi[[k]] %*% GammaPhitList[[k]]
#   }
#   if (estim.varnoise)
#     PhiGammaPhitFull <- PhiGammaPhitFull + varnoise*diag(nrow(PhiGammaPhitFull))
#   
#   Kyy <- PhiGammaPhitFull
#   cholKyy <- chol(PhiGammaPhitFull + model$nugget*diag(nrow(Kyy)))
#   logDetKyy <- 2*sum(diag(log(cholKyy)))
#   invKyy <- chol2inv(cholKyy)
#   f <- 0.5*(logDetKyy + nrow(Kyy)*log(2*pi) + t(model$y)%*%invKyy%*%model$y)
#   
#   # computing the other terms of the constrained log-likelihood
#   GammaPhitInvKyy <- GammaPhit %*% invKyy
#   mu <- GammaPhitInvKyy %*% model$y
#   mu.eta <- model$Lambda %*% mu
#   Sigma <- Gamma - GammaPhitInvKyy %*% t(GammaPhit)
#   Sigma.eta <- model$Lambda %*% Sigma %*% t(model$Lambda)
#   if (min(eigen(Sigma.eta, symmetric = TRUE)$values) <= 0) # numerical stability
#     Sigma.eta <- Sigma.eta + 1e-5*diag(nrow(Sigma.eta))
#   Gamma.eta <- model$Lambda %*% Gamma %*% t(model$Lambda)
#   if (min(eigen(Gamma.eta, symmetric = TRUE)$values) <= 0) # numerical stability
#     Gamma.eta <- Gamma.eta + 1e-5*diag(nrow(Gamma.eta))
#   
#   set.seed(7)
#   switch(mcmc.opts$probe,
#          Genz = {
#            f <- f - log(pmvnorm(model$lb, model$ub,
#                                 as.vector(mu.eta), sigma = Sigma.eta)[[1]]) +
#              log(pmvnorm(model$lb, model$ub,
#                          rep(0, nrow(Gamma.eta)), sigma = Gamma.eta)[[1]])
#          }, ExpT = {
#            f <- f - log(mvNcdf(model$lb-mu.eta, model$ub-mu.eta,
#                                Sigma.eta, mcmc.opts$nb.mcmc)$prob) +
#              log(mvNcdf(model$lb, model$ub, Gamma.eta, mcmc.opts$nb.mcmc)$prob)
#          })
#   
#   
#   return(f)
# }

# #' @title  Numerical Gradient of the Log-Likelihood of a Additive Gaussian Process.
# #' @description Compute the gradient numerically of the negative log-likelihood of an Additive Gaussian Process.
# #' @param par the values of the covariance parameters.
# #' @param model an object with \code{"lineqGP"} S3 class.
# #' @param parfixed indices of fixed parameters to do not be optimised.
# #' @param mcmc.opts not used.
# #' @param estim.varnoise If \code{true}, a noise variance is estimated.
# #' @return the gradient of the negative log-likelihood.
# #'
# #' @seealso \code{\link{logLikAdditiveFun}}, \code{\link{logLikFun}},
# #'          \code{\link{logLikGrad}}
# #' @author A. F. Lopez-Lopera
# #'
# #' @references Rasmussen, C. E. and Williams, C. K. I. (2005),
# #' "Gaussian Processes for Machine Learning (Adaptive Computation and Machine Learning)".
# #' \emph{The MIT Press}.
# #'
# #' @export
# constrlogLikAdditiveGrad <- function(par = unlist(purrr::map(model$kernParam, "par")),
#                                model,  parfixed = rep(FALSE, model$d*length(par)),
#                                mcmc.opts = list(probe = "Genz", nb.mcmc = 1e3),
#                                estim.varnoise = FALSE) {
#   # gradf <- nl.grad(par, constrlogLikFun, heps = 1e-9, model, parfixed, mcmc.opts)
#   # gradf[parfixed] <- 0
#   if (estim.varnoise)
#     parfixed <- c(parfixed, FALSE)
# 
#   # this is a faster version of "nl.grad" when some parameters are fixed.
#   if (!is.numeric(par))
#     stop("Argument 'par' must be a numeric value.")
#   fun <- match.fun(logLikAdditiveFun)
#   fn <- function(x) fun(x, model, parfixed, mcmc.opts, estim.varnoise)
#   if (length(fn(par)) != 1)
#     stop("Function 'f' must be a univariate function of 2 variables.")
#   n <- length(par)
#   hh <- rep(0, n)
#   gradf <- rep(0, n)
# 
#   # the gradient is computed only w.r.t. the parameters to be optimised
#   idx_iter <- seq(n)
#   idx_iter <- idx_iter[parfixed == FALSE]
#   heps <- 1e-9
#   dimIdx <- 0
#   for (i in idx_iter) {
#     hh[i] <- heps
#     gradf[i] <- (fn(par + hh) - fn(par - hh))/(2 * heps)
#     hh[i] <- 0
#   }
#   return(gradf)
# }
# 
