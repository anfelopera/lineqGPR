library("lineqGPR")
library("DiceDesign")
library("plot3D")
library("viridis")
library("Matrix")

rm(list=ls())
set.seed(7)

#### Synthetic data ####
d <- 2
modatan <- function(x, a) return(atan(a*x))
targetFun <- function(x, d) {
  y <- 0
  a <- (1-(1:d)/d)*5
  for (k in 1:d)
    y <- y + modatan(x[, k], a[k])
  return(y)
}

xdesign <- lhsDesign(1e1*d, d, seed = 7)$design
ydesignNoNoise <- targetFun(xdesign, d)
varNoise <- max(range(ydesignNoNoise))*0.01
ydesign <- ydesignNoNoise

#### Constrained model ####
# creating the model
model <- create(class = "lineqAGP", x = xdesign, y = ydesign,
                constrType = rep("monotonicity", d), m = 10)
for (k in 1:d) 
  model$kernParam[[k]]$par <- c(1, 2)
model$localParam$sampler <- "HMC"
model$kernParam
model$nugget <- 1e-7
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
xtestAll <- expand.grid(xtest[,2], xtest[,1])
xtestAll <- xtestAll[d:1]
sim.model <- simulate(model, nsim = nsim, xtest = xtest)
pred <- predict(model, xtest)

colormap <- rev(viridis(1e2))
par(mfrow = c(2,2), mar=c(1.5,1.5,1.5,1))
p <- persp3D(x = unique(xtest[, 1]), y = unique(xtest[, 2]),
             z = outer(modatan(xtest[, 1], (1-1/d)*5),
                       modatan(xtest[, 2], (1-2/d)*5), "+"),
             xlab = "x1", ylab = paste("x", k, sep = ""), zlab = "y(x1,x2)",
             main = "Target function", phi = 20, theta = -30, col = colormap, colkey = FALSE)
points3D(xdesign[,1], xdesign[,2], ydesignNoNoise, pch = 20, cex = 2,  col = "black", add = TRUE)

PhiAll.test <- cbind(pred$Phi.test[[1]][rep(1:ntest, times = ntest), ],
                     pred$Phi.test[[2]][rep(1:ntest, each = ntest), ])
p <- persp3D(x = unique(xtest[, 1]), y = unique(xtest[, 2]),
             z = matrix(PhiAll.test %*% pred$xiAll.map, ntest, ntest),
             xlab = "x1", ylab = paste("x", k, sep = ""), zlab = "y(x1,x2)",
             main = "MAP solution", phi = 20, theta = -30, col = colormap, colkey = FALSE)
points3D(xdesign[,1], xdesign[,2], ydesignNoNoise, pch = 20, cex = 2,  col = "black", add = TRUE)

p <- persp3D(x = unique(xtest[, 1]), y = unique(xtest[, 2]),
             z = matrix(rowMeans(PhiAll.test %*% sim.model$xiAll.sim), ntest, ntest),
             xlab = "x1", ylab = paste("x", k, sep = ""), zlab = "y(x1,x2)",
             main = "MCMC solution", phi = 20, theta = -30, col = colormap, colkey = FALSE)
points3D(xdesign[,1], xdesign[,2], ydesignNoNoise, pch = 20, cex = 2,  col = "black", add = TRUE)


p <- persp3D(x = unique(xtest[, 1]), y = unique(xtest[, 2]),
             z = matrix(PhiAll.test %*% sim.model$xiAll.sim[,1], ntest, ntest),
             xlab = "x1", ylab = paste("x", k, sep = ""), zlab = "y(x1,x2)",
             main = "MCMC sample", phi = 20, theta = -30, col = colormap, colkey = FALSE)
points3D(xdesign[,1], xdesign[,2], ydesignNoNoise, pch = 20, cex = 2,  col = "black", add = TRUE)

