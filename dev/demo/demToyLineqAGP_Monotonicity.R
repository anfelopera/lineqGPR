library("lineqGPR")
library("DiceDesign")
library("plot3D")
library("viridis")
library("Matrix")

rm(list=ls())
set.seed(7)

#### Synthetic data ####
d <- 30
modatan <- function(x, a) return(atan(a*x))
targetFun <- function(x, d) {
  y <- 0
  a <- (1-(1:d)/d)*5
  for (k in 1:d)
    y <- y + modatan(x[, k], a[k])
  return(y)
}

xdesign <- lhsDesign(2e2, d, seed = 7)$design
ydesignNoNoise <- targetFun(xdesign, d)
varNoise <- max(range(ydesignNoNoise))*0.01
ydesign <- ydesignNoNoise

#### Constrained model ####
# creating the model
model <- create(class = "lineqAGP", x = xdesign, y = ydesign,
                constrType = rep("monotonicity", d), m = 5)
for (k in 1:d) 
  model$kernParam[[k]]$par <- c(1, 2)
model$nugget <- 1e-3
model$varnoise <- 1e-2 #varNoise

# training the model
# Note: the covariance parameter estimation takes almost 0.5-1 hour
opttime <- proc.time()
model <- lineqGPOptim(model,
                      x0 = unlist(purrr::map(model$kernParam, "par")),
                      eval_f = "logLik",
                      additive = TRUE,
                      opts = list(algorithm = "NLOPT_LD_MMA",
                                  print_level = 3,
                                  ftol_abs = 1e-3,
                                  maxeval = 1e2,
                                  check_derivatives = FALSE),
                      lb = rep(1e-2, 2*d), ub = rep(c(5, 3), d),
                      estim.varnoise = TRUE,
                      bounds.varnoise = c(1e-2, Inf))
opttime <- proc.time() - opttime

# simulating samples from the model
ntest <- 10
xtest  <- matrix(seq(0, 1, length = ntest), nrow = ntest, ncol = d)
ytest <- targetFun(xtest, d)

nsim <- 1e3
sim.model <- simulate(model, nsim = nsim, xtest = xtest)#, seed = 7)

msum <- cumsum(model$localParam$m)
colormap <- rev(viridis(1e2))
for (k in seq(5, d, 5)) {
  par(mfrow = c(2,2), mar=c(1.5,1.5,1.5,1))
  p <- persp3D(x = unique(xtest[, 1]), y = unique(xtest[, 2]),
               z = outer(modatan(xtest[, 1], (1-1/d)*5),
                         modatan(xtest[, k], (1-k/d)*5), "+"),
               xlab = "x1", ylab = paste("x", k, sep = ""), zlab = "y(x1,..., xd)",
               main = "Target function", phi = 20, theta = -30, col = colormap, colkey = FALSE)

  # p <- persp3D(x = unique(xtest[, 1]), y = unique(xtest[, 2]),
  #              z = outer(rowMeans(sim.model$ysim[[1]]),
  #                        rowMeans(sim.model$ysim[[k]]), "+"),
  #              xlab = "x1", ylab = paste("x", k, sep = ""), zlab = "u(x1,..., xd)",
  #              main = "Predictive mean", phi = 20, theta = -30, col = colormap, colkey = FALSE)
  
  PhiAll.test <- cbind(sim.model$Phi.test[[1]][rep(1:ntest, times = ntest), ],
                       sim.model$Phi.test[[k]][rep(1:ntest, each = ntest), ])
  
  p <- persp3D(x = unique(xtest[, 1]), y = unique(xtest[, 2]),
               z = matrix(PhiAll.test %*% sim.model$xiAll.map[c(1:msum[1], (msum[k-1]+1):msum[k])], ntest, ntest),
               xlab = "x1", ylab = paste("x", k, sep = ""), zlab = "y(x1,..., xd)",
               main = "MAP solution", phi = 20, theta = -30, col = colormap, colkey = FALSE)
  
  p <- persp3D(x = unique(xtest[, 1]), y = unique(xtest[, 2]),
               z = matrix(rowMeans(PhiAll.test %*% sim.model$xiAll.sim[c(1:msum[1], (msum[k-1]+1):msum[k]), ]), ntest, ntest),
               xlab = "x1", ylab = paste("x", k, sep = ""), zlab = "y(x1,..., xd)",
               main = "MCMC solution", phi = 20, theta = -30, col = colormap, colkey = FALSE)
  points3D(xdesign[,1], xdesign[,2], ydesignNoNoise, pch = 20, cex = 2,  col = "black", add = TRUE)
  
  
  p <- persp3D(x = unique(xtest[, 1]), y = unique(xtest[, 2]),
               z = matrix(PhiAll.test %*% sim.model$xiAll.sim[c(1:msum[1], (msum[k-1]+1):msum[k]), 2], ntest, ntest),
               xlab = "x1", ylab = paste("x", k, sep = ""), zlab = "y(x1,..., xd)",
               main = "MCMC sample", phi = 20, theta = -30, col = colormap, colkey = FALSE)
  points3D(xdesign[,1], xdesign[,2], ydesignNoNoise, pch = 20, cex = 2,  col = "black", add = TRUE)
  
  
}

