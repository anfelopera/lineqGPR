library("lineqGPR")
library("DiceDesign")
library("plot3D")
library("viridis")
library("Matrix")

rm(list=ls())
set.seed(7)

#### Synthetic data ####
d <- 2
targetFun <- function(x) return(atan(5*x[,2]) + 0.5*x[,1])
# targetFun <- function(x) return(3*x[,1] + 0.5*x[,2])

xdesign <- lhsDesign(2e1*d, d, seed = 7)$design
# xdesign <- factDesign(2, 6)$design
ydesignNoNoise <- targetFun(xdesign)
varNoise <- max(range(ydesignNoNoise))*0.01
ydesign <- ydesignNoNoise

#### Constrained model ####
# creating the model
model <- create(class = "lineqAGP", x = xdesign, y = ydesign,
                constrType = rep("monotonicity", d), m = 5)
for (k in 1:d) 
  model$kernParam[[k]]$par <- c(1, 2)
model$localParam$sampler <- "HMC"
model$nugget <- 1e-7
model$varnoise <- 0.01*sd(ydesign)^2

# ntest <- 1e5
# xtestGrid <- matrix(runif(2*ntest), ntest, 2)
ntest <- 50
xtestGrid  <- as.matrix(expand.grid(seq(0, 1, length = ntest), seq(0, 1, length = ntest)))
ytestGrid <- targetFun(xtestGrid)

# evaluating the model
model <- AdditiveMaxMod(model,
                        xtestGrid,
                        tol = 1e-4, # numerical stability issue
                        max_iter = 10,
                        reward_new_knot = 1e-9,
                        reward_new_dim = 1e-9,
                        print_iter = TRUE,
                        nClusters = 1,
                        save_history = TRUE)
message("\nNumber of active dimensions: ", d)
message("Number of actived dimensions via MaxMod: ", model$d, "\n")
idxAdd <- unique(model$MaxMod$optDecision)

pred <- predict(model, xtestGrid)
model$localParam$sampler <- "HMC"
sim.model <- simulate(model, nsim = 1e3, xtest = xtestGrid)

u <- expand.grid(model$ulist[[1]], model$ulist[[2]])
pred_Knots <- predict(model, as.matrix(u))
PhiAllknots.test <- cbind(pred_Knots$Phi.test[[1]][rep(1:model$localParam$m[1], times = model$localParam$m[2]), ],
                          pred_Knots$Phi.test[[2]][rep(1:model$localParam$m[2], each = model$localParam$m[1]), ])

colormap <- rev(viridis(1e2))
par(mfrow = c(2,2), mar=c(1.5,1.5,1.5,1))
p <- persp3D(x = unique(xtestGrid[, idxAdd[1]]), y = unique(xtestGrid[, idxAdd[2]]),
             z = matrix(targetFun(xtestGrid[, idxAdd]), ntest, ntest),
             xlab = paste("x", idxAdd[1], sep = ""), ylab = paste("x", idxAdd[2], sep = ""), zlab = "y(x1,x2)",
             main = "Target function", phi = 20, theta = -30, col = colormap, colkey = FALSE)
points3D(xdesign[,idxAdd[1]], xdesign[,idxAdd[2]], ydesignNoNoise, pch = 20, cex = 2,  col = "black", add = TRUE)

p <- persp3D(x = unique(xtestGrid[, 1]), y = unique(xtestGrid[, 2]),
             z = matrix(pred$PhiAll.test %*% pred$xiAll.map, ntest, ntest),
             xlab = paste("x", idxAdd[1], sep = ""), ylab = paste("x", idxAdd[2], sep = ""), zlab = "y(x1,x2)",
             main = "MAP solution", phi = 20, theta = -30, col = colormap, colkey = FALSE)
points3D(xdesign[,idxAdd[1]], xdesign[,idxAdd[2]], ydesignNoNoise, pch = 20, cex = 2,  col = "black", add = TRUE)
points(trans3D(x = u[,1], y = u[,2], z = pred_Knots$PhiAll.test %*% pred_Knots$xiAll.map, pmat = p),
       col = 'brown', pch = 4, lwd = 2)

p <- persp3D(x = unique(xtestGrid[, 1]), y = unique(xtestGrid[, 2]),
             z = matrix(rowMeans(sim.model$PhiAll.test %*% sim.model$xiAll.sim), ntest, ntest),
             xlab = paste("x", idxAdd[1], sep = ""), ylab = paste("x", idxAdd[2], sep = ""), zlab = "y(x1,x2)",
             main = "MCMC solution", phi = 20, theta = -30, col = colormap, colkey = FALSE)
points3D(xdesign[,idxAdd[1]], xdesign[,idxAdd[2]], ydesignNoNoise, pch = 20, cex = 2,  col = "black", add = TRUE)
points(trans3D(x = u[,1], y = u[,2], z = pred_Knots$PhiAll.test %*% pred_Knots$xiAll.map, pmat = p),
       col = 'brown', pch = 4, lwd = 2)

p <- persp3D(x = unique(xtestGrid[, 1]), y = unique(xtestGrid[, 2]),
             z = matrix(sim.model$PhiAll.test %*% sim.model$xiAll.sim[,1], ntest, ntest),
             xlab = paste("x", idxAdd[1], sep = ""), ylab = paste("x", idxAdd[2], sep = ""), zlab = "y(x1,x2)",
             main = "MCMC sample", phi = 20, theta = -30, col = colormap, colkey = FALSE)
points3D(xdesign[,idxAdd[1]], xdesign[,idxAdd[2]], ydesignNoNoise, pch = 20, cex = 2,  col = "black", add = TRUE)
points(trans3D(x = u[,1], y = u[,2], z = pred_Knots$PhiAll.test %*% pred_Knots$xiAll.map, pmat = p),
       col = 'brown', pch = 4, lwd = 2)

