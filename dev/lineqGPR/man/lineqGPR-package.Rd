\name{lineqGPR-package}
\encoding{UTF8}
\alias{lineqGPR-package}
\alias{lineqGPR}
\docType{package}
\title{

  Gaussian Processes with Linear Inequality Constraints

}
\description{

  A package for Gaussian process interpolation, regression and simulation
  under linear inequality constraints based on (López-Lopera et al., 2018).
  Constrained models and constrained additive models are given as objects
  with "lineqGP" and "lineqAGP" S3 class, respectively. Implementations
  according to (Maatouk and Bay, 2017) are also provided as objects with
  "lineqDGP" S3 class.

}

\note{

  This package was initially developed within the frame of the Chair in Applied
  Mathematics OQUAIDO, gathering partners in technological research (BRGM,
  CEA, IFPEN, IRSN, Safran, Storengy) and academia (CNRS, Ecole Centrale
  de Lyon, Mines Saint-Etienne, University of Grenoble, University of Nice,
  University of Toulouse) around advanced methods for Computer Experiments. \cr
  
  \tabular{ll}{
  \tab \bold{Important functions or methods} \cr \cr
  \code{create} \tab Creation function of GP models under linear inequality constraints. \cr
  \code{augment} \tab Augmentation of GP models according to local and covariance parameters. \cr
  \code{lineqGPOptim} \tab Covariance parameter estimation via maximum likelihood. \cr
  \code{predict} \tab Prediction of the objective function at new points using a Kriging model under \cr
  \tab linear inequality constraints. \cr
  \code{simulate} \tab Simulation of kriging models under linear inequality constraints. \cr
  \code{plot} \tab Plot for a constrained Kriging model. \cr
  \code{ggplot} \tab GGPlot for a constrained Kriging model. \cr  
  }
}

\details{
  \packageDESCRIPTION{lineqGPR}
}

\section{Warning}{

  \pkg{lineqGPR} may strongly evolve in the future in order to incorporate
  other packages for Gaussian process regression modelling (see, e.g.,
  \pkg{kergp}, \pkg{DiceKriging}, \pkg{DiceDesign}). It could be also scaled
  to higher dimensions and for a large number of observations.

}

\author{

  Andrés F. López-Lopera (UPHF, Valenciennes)
  with contributions from
  Olivier Roustant (INSA, Toulouse) and
  Yves Deville (Alpestat).

  Maintainer: Andrés F. López-Lopera, <andres.lopezlopera@uphf.fr>

}

\references{

  A. F. López-Lopera, F. Bachoc, N. Durrande and O. Roustant (2018),
  "Finite-dimensional Gaussian approximation with linear inequality constraints".
  \emph{SIAM/ASA Journal on Uncertainty Quantification}, 6(3): 1224–1255.
  <doi:10.1137/17M1153157>

  F. Bachoc, A. Lagnoux and A. F. López-Lopera (2019),
  "Maximum likelihood estimation for Gaussian processes under inequality constraints".
  \emph{Electronic Journal of Statistics}, 13 (2): 2921-2969.
  <doi:10.1214/19-EJS1587>
  
  F. Bachoc, A. F. López-Lopera and O. Roustant (2020),
  "Sequential construction and dimension reduction of Gaussian processes under inequality constraints".
  \emph{SIAM Journal on Mathematics of Data Science}, 4(2): 772-800.
  <doi:10.1137/21M1407513>
  
  A. F. López-Lopera, F. Bachoc and O. Roustant (2022),
  "High-dimensional additive Gaussian processes under monotonicity constraints".
  \emph{NeurIPS}
  <arXiv:2009.04188>
  
  H. Maatouk and X. Bay (2017),
  "Gaussian process emulators for computer experiments with inequality constraints".
  \emph{Mathematical Geosciences}, 49(5): 557-582.

  Roustant, O., Ginsbourger, D., and Deville, Y. (2012),
  "DiceKriging, DiceOptim: Two R Packages for the Analysis of
  Computer Experiments by Kriging-Based Metamodeling and Optimization".
  \emph{Journal of Statistical Software}, 51(1): 1-55.

}
%%\keyword{ package }
%%\seealso{}

\examples{
## ------------------------------------------------------------------
## Gaussian process regression modelling under boundedness constraint
## ------------------------------------------------------------------
library(lineqGPR)

#### generating the synthetic data ####
sigfun <- function(x) return(1/(1+exp(-7*(x-0.5))))
x <- seq(0, 1, 0.001)
y <- sigfun(x)
DoE <- splitDoE(x, y, DoE.idx = c(201, 501, 801))

#### GP with inactive boundedness constraints ####
# creating the "lineqGP" model
model <- create(class = "lineqGP", x = DoE$xdesign, y = DoE$ydesign,
                constrType = c("boundedness"), m = 100)
model$bounds <- c(-10,10)

# sampling from the model
sim.model <- simulate(model, DoE$xtest, nsim = 1e3, seed = 1)
plot(sim.model, xlab = "x", ylab = "y(x)", ylim = range(y),
     main = "Unconstrained GP model")
lines(x, y, lty = 2)
legend("topleft", c("ytrain","ytest","mean","confidence"),
       lty = c(NaN,2,1,NaN), pch = c(20,NaN,NaN,15),
       col = c("black","black","darkgreen","gray80"))

#### GP with active boundedness constraints ####
# creating the "lineqGP" model
model <- create(class = "lineqGP", x = DoE$xdesign, y = DoE$ydesign,
                constrType = c("boundedness"), m = 100)
model$bounds <- c(0,1)

# sampling from the model
sim.model <- simulate(model, DoE$xtest, nsim = 1e2, seed = 1)
plot(sim.model, bounds = model$bounds,
     xlab = "x", ylab = "y(x)", ylim = range(y),
     main = "Constrained GP model under boundedness conditions")
lines(x, y, lty = 2)
legend("topleft", c("ytrain","ytest","mean","confidence"),
       lty = c(NaN,2,1,NaN), pch = c(20,NaN,NaN,15),
       col = c("black","black","darkgreen","gray80"))


## ------------------------------------------------------------------
## Gaussian process regression modelling under multiple constraints
## ------------------------------------------------------------------
library(lineqGPR)

#### generating the synthetic data ####
sigfun <- function(x) return(1/(1+exp(-7*(x-0.5))))
x <- seq(0, 1, 0.001)
y <- sigfun(x)
DoE <- splitDoE(x, y, DoE.idx = c(201, 501, 801))

#### GP with boundedness and monotonicity constraints ####
# creating the "lineqGP" model
model <- create(class = "lineqGP", x = DoE$xdesign, y = DoE$ydesign,
                constrType = c("boundedness","monotonicity"), m = 30)
model$bounds[1, ] <- c(0,1)
model$varnoise <- 1e-5

# sampling from the model
sim.model <- simulate(model, DoE$xtest, nsim = 1e2, seed = 1)
plot(sim.model, bounds = model$bounds,
     xlab = "x", ylab = "y(x)", ylim = range(y),
     main = "Constrained GP model under boundedness & monotonicity conditions")
lines(x, y, lty = 2)
legend("topleft", c("ytrain","ytest","mean","confidence"),
       lty = c(NaN,2,1,NaN), pch = c(20,NaN,NaN,15),
       col = c("black","black","darkgreen","gray80"))


## ------------------------------------------------------------------
## Gaussian process regression modelling under linear constraints
## ------------------------------------------------------------------
library(lineqGPR)
library(Matrix)

#### generating the synthetic data ####
targetFun <- function(x){
  y <- rep(1, length(x))
  y[x <= 0.4] <- 2.5*x[x <= 0.4]
  return(y)
}
x <- seq(0, 1, by = 0.001)
y <- targetFun(x)
DoE <- splitDoE(x, y, DoE.idx = c(101, 301, 501, 701))

#### GP with predefined linear inequality constraints ####
# creating the "lineqGP" model
model <- create(class = "lineqGP", x = DoE$xdesign, y = DoE$ydesign,
                constrType = c("linear"), m = 100)
m <- model$localParam$m
model$localParam$sampler <- "HMC"

# building the predefined linear constraints
bounds1 <- c(0,Inf)
LambdaB1 <- diag(2*m/5)
LambdaM <- diag(2*m/5)
LambdaB2 <- diag(3*m/5)
lsys <- lineqGPSys(m = 2*m/5, constrType = "monotonicity",
                   l = bounds1[1], u = bounds1[2], lineqSysType = "oneside")
LambdaM[-seq(1),] <- lsys$M
model$Lambda <- as.matrix(bdiag(rbind(LambdaM,LambdaB1),LambdaB2))
model$lb <- c(-Inf, rep(0, 2*m/5-1), rep(0, 2*m/5), rep(0.85, 3*m/5))
model$ub <- c(rep(0.1, 2*m/5), rep(1.1, 2*m/5), rep(1.1, 3*m/5))

# sampling from the model
sim.model <- simulate(model, DoE$xtest, nsim = 1e2, seed = 1)
plot(sim.model, bounds = c(0,1.1),
     xlab = "x", ylab = "y(x)", ylim = c(0,1.1),
     main = "Constrained GP model under linear conditions")
lines(x, y, lty = 2)
abline(v = 0.4, lty = 2)
lines(c(0.4, 1), rep(0.85, 2), lty = 2)
legend("bottomright", c("ytrain","ytest","mean","confidence"),
       lty = c(NaN,2,1,NaN), pch = c(20,NaN,NaN,15),
       col = c("black","black","darkgreen","gray80"))


## ------------------------------------------------------------------
## Note:
## 1. More examples are given as demos (run: demo(package="lineqGPR")).
## 2. See also the examples from inner functions of the package
## (run: help("simulate.lineqGP")).
## ------------------------------------------------------------------

}
