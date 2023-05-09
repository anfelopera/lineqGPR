#' @title Kernel Matrix for \code{"lineqGP"} Models.
#' @description Compute the kernel matrix for \code{"lineqGP"} models.
#' attr: "gradient".
#' 
#' @param x1 a vector with the first input locations.
#' @param x2 a vector with the second input locations.
#' @param type a character string corresponding to the type of the kernel.
#' Options: "gaussian", "matern32", "matern52", "exponential".
#' @param par the values of the kernel parameters (variance, lengthscale).
#' @param d a number corresponding to the dimension of the input space.
#' 
#' @return Kernel matrix \eqn{K(x_1,x_2)}{K(x1,x2)}
#' (or \eqn{K(x_1,x_1)}{K(x1,x1)} if \eqn{x_2}{x2} is not defined).
#'
#' @author A. F. Lopez-Lopera
#'
#' @examples
#' x <- seq(0, 1, 0.01)
#' K <- kernCompute(x, type = "gaussian", par =  c(1, 0.1))
#' image(K, main = "covariance matrix")
#'
#' @export
kernCompute <- function(x1, x2 = NULL, type, par, d = 1L) {
  kernName <- paste("k", d, type, sep = "")
  kernFun <- try(get(kernName))
  if (inherits(kernFun,  "try-error")) {
    stop('kernel "', type, '" is not supported')
  } else {
    if (!is.null(x2)) {
      kern <- kernFun(x1, x2, par, d)
    } else {
      kern <- kernFun(x1, x1, par, d)
    }
    attr(kern, "gradient") <- attr(kern, "gradient")
    return(kern)
  }
}

#' @title 1D Gaussian Kernel Matrix for \code{"lineqGP"} Models.
#' @description Compute the 1D Gaussian kernel matrix for \code{"lineqGP"} models.
#' attr: "gradient", "derivative".
#' 
#' @param x1 a vector with the first input locations.
#' @param x2 a vector with the second input locations.
#' @param par the values of the kernel parameters (variance, lengthscale).
#' @param d a number corresponding to the dimension of the input space.
#' 
#' @return Kernel matrix \eqn{K(x_1,x_2)}{K(x1,x2)}
#' (or \eqn{K(x_1,x_1)}{K(x1,x1)} if \eqn{x_2}{x2} is not defined).
#'
#' @author A. F. Lopez-Lopera
#'
#' @examples
#' x <- seq(0, 1, 0.01)
#' K <- k1gaussian(x, x, par =  c(1, 0.1))
#' image(K, main = "covariance matrix using a Squared Exponential kernel")
#'
#' @export
k1gaussian <- function(x1, x2, par, d = 1) {
  if (!is.matrix(x1) || ncol(x1) != d) x1 <- matrix(x1, ncol = 1)
  if (!is.matrix(x2) || ncol(x2) != d) x2 <- matrix(x2, ncol = 1)
  sigma2 <- par[1]
  theta <- par[2]
  dist <- outer(x1[, 1]/theta, x2[, 1]/theta, "-")
  dist2 <- dist^2
  kern <- sigma2*exp(-0.5*dist2)

  # gradients
  dsigma2 <- kern/sigma2
  dtheta <- kern*dist2/theta
  attr(kern, "gradient") <- list(sigma2 = dsigma2, theta = dtheta)

  # derivatives w.r.t. the input variable
  d1x1 <- -(dist/theta)*kern
  d2x1 <- (-1 + dist2)*kern/theta^2
  d3x1 <- dist*(3 - dist2)*kern/theta^3
  d4x1 <- (3 - 6*dist2 + dist2^2)*kern/theta^4
  attr(kern, "derivative") <- list(x1 = d1x1, x2 = -d1x1,
                                   x1x1 = d2x1, x1x2 = -d2x1,
                                   x2x2 = d2x1, x2x1 = -d2x1,
                                   x1x1x1 = d3x1, x1x1x2 = -d3x1, x2x2x1 = d3x1,
                                   x1x1x1x1 = d4x1, x1x1x2x2 = d4x1)
  return(kern)
}

#' @title 2D Gaussian Kernel Matrix for \code{"lineqGP"} Models.
#' @description Compute the 2D Gaussian kernel matrix for \code{"lineqGP"} models.
#' attr: "gradient".
#' 
#' @param x1 a matrix with the first couple of input locations.
#' @param x2 a matrix with the second couple of input locations.
#' @param par the values of the kernel parameters (variance, lengthscales).
#' @param d a number corresponding to the dimension of the input space.
#' 
#' @return Kernel matrix \eqn{K(x_1,x_2)}{K(x1,x2)}
#' (or \eqn{K(x_1,x_1)}{K(x1,x1)} if \eqn{x_2}{x2} is not defined).
#'
#' @author A. F. Lopez-Lopera
#'
#' @examples
#' xgrid <- seq(0, 1, 0.1)
#' x <- as.matrix(expand.grid(xgrid, xgrid))
#' K <- k2gaussian(x, x, par =  c(1, 0.1))
#' image(K, main = "covariance matrix using a 2D Gaussian kernel")
#'
#' @export
k2gaussian <- function(x1, x2, par, d = 2) {
  if (!is.matrix(x1) || ncol(x1) != d) x1 <- matrix(x1, ncol = 1)
  if (!is.matrix(x2) || ncol(x2) != d) x2 <- matrix(x2, ncol = 1)
  if (length(par) == 2)
    par <- c(par, par[2])

  sigma2 <- par[1]
  theta_x1 <- par[2]
  theta_x2 <- par[3]

  dist2_x1 <- outer(x1[, 1]/theta_x1, x2[, 1]/theta_x1,'-')^2
  dist2_x2 <- outer(x1[, 2]/theta_x2, x2[, 2]/theta_x2,'-')^2
  kern <- sigma2*exp(-0.5*dist2_x1 - 0.5*dist2_x2)

  # gradients
  attr(kern, "gradient") <- list(sigma2 = kern/sigma2,
                                 theta_x1 = 2*kern*dist2_x1/theta_x1,
                                 theta_x2 = 2*kern*dist2_x2/theta_x2)
  return(kern)
}

#' @title 1D Matern 3/2 Kernel Matrix for \code{"lineqGP"} Models.
#' @description Compute the 1D Matern 3/2 kernel for \code{"lineqGP"} models.
#' attr: "gradient", "derivative".
#' 
#' @param x1 a vector with the first input locations.
#' @param x2 a vector with the second input locations.
#' @param par the values of the kernel parameters (variance, lengthscale).
#' @param d a number corresponding to the dimension of the input space.
#' 
#' @return Kernel matrix \eqn{K(x_1,x_2)}{K(x1,x2)}
#' (or \eqn{K(x_1,x_1)}{K(x1,x1)} if \eqn{x_2}{x2} is not defined).
#'
#' @author A. F. Lopez-Lopera
#'
#' @examples
#' x <- seq(0, 1, 0.01)
#' K <- k1matern32(x, x, par =  c(1, 0.1))
#' image(K, main = "covariance matrix using a Matern 3/2 kernel")
#'
#' @export
k1matern32 <- function(x1, x2, par, d = 1) {
  if (!is.matrix(x1) || ncol(x1) != d) x1 <- matrix(x1, ncol = 1)
  if (!is.matrix(x2) || ncol(x2) != d) x2 <- matrix(x2, ncol = 1)
  sigma2 <- par[1]
  theta <- par[2]
  dist <- outer(x1[, 1]/theta, x2[, 1]/theta, "-")
  distSgn <- abs(dist)
  sqrt3distSgn <- sqrt(3)*distSgn
  expDistSgn <- exp(-sqrt3distSgn)
  kern <- sigma2*(1 + sqrt3distSgn)*expDistSgn

  # gradients
  dsigma2 <- kern/sigma2
  dtheta <- (-sigma2*expDistSgn + kern)*(sqrt3distSgn/theta)
  attr(kern, "gradient") <- list(sigma2 = dsigma2, theta = dtheta)

  # derivatives w.r.t. the input variable
  d1x1 <- -3*(sigma2/theta)*dist*expDistSgn
  d2x1 <- 3*(sigma2/theta^2)*(-1 + sqrt3distSgn)*expDistSgn
  attr(kern, "derivative") <- list(x1 = d1x1, x2 = -d1x1,
                                   x1x1 = d2x1, x1x2 = -d2x1,
                                   x2x2 = d2x1, x2x1 = -d2x1)
  return(kern)
}

#' @title 1D Matern 5/2 Kernel Matrix for \code{"lineqGP"} Models.
#' @description Compute the 1D Matern 5/2 kernel for \code{"lineqGP"} models.
#' attr: "gradient", "derivative".
#' 
#' @param x1 A vector with the first input locations.
#' @param x2 A vector with the second input locations.
#' @param par Values of the kernel parameters (variance, lengthscale).
#' @param d A number corresponding to the dimension of the input space.
#' 
#' @return Kernel matrix \eqn{K(x_1,x_2)}{K(x1,x2)}
#' (or \eqn{K(x_1,x_1)}{K(x1,x1)} if \eqn{x_2}{x2} is not defined).
#'
#' @author A. F. Lopez-Lopera
#'
#' @examples
#' x <- seq(0, 1, 0.01)
#' K <- k1matern52(x, x, par =  c(1, 0.1))
#' image(K, main = "covariance matrix using a Matern 5/2 kernel")
#'
#' @export
k1matern52 <- function(x1, x2, par, d = 1) {
  if (!is.matrix(x1) || ncol(x1) != d) x1 <- matrix(x1, ncol = 1)
  if (!is.matrix(x2) || ncol(x2) != d) x2 <- matrix(x2, ncol = 1)
  sigma2 <- par[1]
  theta <- par[2]
  dist <- outer(x1[, 1]/theta, x2[, 1]/theta, "-")
  dist2 <- dist^2
  distSgn <- abs(dist)
  sqrt5distSgn <- sqrt(5)*distSgn
  expDistSgn <- exp(-sqrt5distSgn)
  kern <- sigma2*(1 + sqrt5distSgn + (5/3)*dist2)*expDistSgn

  # gradients
  dsigma2 <- kern/sigma2
  dtheta <- -(sigma2/theta)*(sqrt5distSgn + (10/3)*dist2)*expDistSgn +
    (sqrt5distSgn/theta)*kern
  attr(kern, "gradient") <- list(sigma2 = dsigma2, theta = dtheta)

  # derivatives w.r.t. the input variable
  sqrt5distSgn2 <- sqrt(5)*sign(dist)*dist2
  const  <- (5/3)*sigma2/theta
  d1x1 <- const*(-dist - sqrt5distSgn2)*expDistSgn
  d2x1 <- (const/theta)*(-1 - sqrt5distSgn + 5*dist2)*expDistSgn
  d3x1 <- (5*const/theta^2)*(3*dist - sqrt5distSgn2)*expDistSgn
  d4x1 <- (5*const/theta^3)*(3 - 5*sqrt5distSgn + 5*dist2)*expDistSgn
  attr(kern, "derivative") <- list(x1 = d1x1, x2 = -d1x1,
                                   x1x1 = d2x1, x1x2 = -d2x1,
                                   x2x2 = d2x1, x2x1 = -d2x1,
                                   x1x1x1 = d3x1, x1x1x2 = -d3x1, x2x2x1 = d3x1,
                                   x1x1x1x1 = d4x1, x1x1x2x2 = d4x1)
  return(kern)
}

#' @title 1D Exponential Kernel Matrix for \code{"lineqGP"} Models.
#' @description Compute the 1D Exponential kernel for \code{"lineqGP"} models.
#' attr: "gradient".
#' 
#' @param x1 a vector with the first input locations.
#' @param x2 a vector with the second input locations.
#' @param par the values of the kernel parameters (variance, lengthscale).
#' @param d a number corresponding to the dimension of the input space.
#' 
#' @return Kernel matrix \eqn{K(x_1,x_2)}{K(x1,x2)}
#' (or \eqn{K(x_1,x_1)}{K(x1,x1)} if \eqn{x_2}{x2} is not defined).
#'
#' @author A. F. Lopez-Lopera
#'
#' @examples
#' x <- seq(0, 1, 0.01)
#' K <- k1exponential(x, x, par =  c(1, 0.1))
#' image(K, main = "covariance matrix using a Exponential kernel")
#'
#' @export
k1exponential <- function(x1, x2, par, d = 1) {
  if (!is.matrix(x1) || ncol(x1) != d) x1 <- matrix(x1, ncol = 1)
  if (!is.matrix(x2) || ncol(x2) != d) x2 <- matrix(x2, ncol = 1)
  sigma2 <- par[1]
  theta <- par[2]
  dist <- outer(x1[, 1]/theta, x2[, 1]/theta, "-")
  distSgn <- abs(dist)
  kern <- sigma2*exp(-distSgn)

  # gradients
  dsigma2 <- kern/sigma2
  dtheta <- kern*distSgn/theta
  attr(kern, "gradient") <- list(sigma2 = dsigma2, theta = dtheta)
  return(kern)
}
