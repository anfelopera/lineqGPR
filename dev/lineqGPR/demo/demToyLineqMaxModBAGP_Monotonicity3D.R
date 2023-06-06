library("lineqGPR")
library("DiceDesign")
library("plot3D")
library("viridis")
library("Matrix")

rm(list=ls())
set.seed(7)

#### Synthetic data ####
targetFun <- function(x) {
  return(x[, 1]*x[, 3] + x[, 2])
}
D <- 3 # number of active + inactive input variables
d <- 2 # number of active input variables

# building DoE in dimension D via Latin hypercube sampling (LHS)
nbtrain <- 10*D
xdesign <- lhsDesign(nbtrain, D, seed = 8)$design
xdesign <- maximinSA_LHS(xdesign)$design
ydesign <- targetFun(xdesign)

# building a DoE for assessing the model
ntest <- 40*D
xtest <- lhsDesign(ntest, D, seed = 8)$design
xtest <- maximinSA_LHS(xtest)$design

#### Constrained model ####
# creating the model
nblock <- 2
model <- create(class = "lineqBAGP", x = xdesign, y = ydesign,
                constrType = rep("monotonicity", nblock), 
                partition = list(c(1,3), 2),
                m = c(5, 3, 4))

for (k in 1:nblock)
  model$kernParam[[k]]$par <- c(1, rep(2, model$localParam$dblock[k]))

model$localParam$sampler <- "HMC"
model$nugget <- 1e-5
model$varnoise <- 0.05*sd(ydesign)^2 

model_test <- augment(model)

# model_temp <- lineqGPOptim(model,
#                            additive = TRUE,
#                            partition = list(c(1,3), 2),
#                            estim.varnoise = TRUE, # to add this info at the MaxMod level
#                            bounds.varnoise = c(1e-7, Inf), # to add this info at the MaxMod level
#                            lb = rep(1e-2, model_temp$d+1), ub = c(Inf, rep(1, model_temp$d)) # to add this info at the MaxMod level
# )

# # evaluating the model
# model <- AdditiveMaxMod(model,
#                         xtest,
#                         tol = 1e-4, # numerical stability issue
#                         max_iter = 9,
#                         reward_new_knot = 1e-12,
#                         reward_new_dim = 1e-12,
#                         print_iter = TRUE,
#                         nClusters = 5,
#                         save_history = TRUE)
# message("\nNumber of active dimensions: ", d)
# message("Number of actived dimensions via MaxMod: ", model$d, "\n")
# idxAdd <- unique(model$MaxMod$optDecision)

ntest <- 10
xtestGrid  <- as.matrix(expand.grid(seq(0, 1, length = ntest), seq(0, 1, length = ntest)))

pred <- predict(model, xtestGrid)
model$localParam$sampler <- "HMC"
sim.model <- simulate(model, nsim = 1e2, xtest = xtestGrid)

u <- expand.grid(model$ulist[[1]], model$ulist[[2]])
pred_Knots <- predict(model, as.matrix(u))
PhiAllknots.test <- cbind(pred_Knots$Phi.test[[1]][rep(1:model$localParam$m[1], times = model$localParam$m[2]), ],
                          pred_Knots$Phi.test[[2]][rep(1:model$localParam$m[2], each = model$localParam$m[1]), ])

colormap <- rev(viridis(1e2))
par(mfrow = c(2,2), mar=c(1.5,1.5,1.5,1))
p <- persp3D(x = unique(xtestGrid[, idxAdd[1]]), y = unique(xtestGrid[, idxAdd[2]]),
             z = matrix(targetFun(cbind(xtestGrid[, idxAdd], 1), d), ntest, ntest),
             xlab = paste("x", idxAdd[1], sep = ""), ylab = paste("x", idxAdd[2], sep = ""), zlab = "y(x1,x2)",
             main = "Target function", phi = 20, theta = -30, col = colormap, colkey = FALSE)
points3D(xdesign[,idxAdd[1]], xdesign[,idxAdd[2]], ydesign, pch = 20, cex = 2,  col = "black", add = TRUE)

p <- persp3D(x = unique(xtestGrid[, 1]), y = unique(xtestGrid[, 2]),
             z = matrix(pred$PhiAll.test %*% pred$xiAll.map, ntest, ntest),
             xlab = paste("x", idxAdd[1], sep = ""), ylab = paste("x", idxAdd[2], sep = ""), zlab = "y(x1,x2)",
             main = "MAP solution", phi = 20, theta = -30, col = colormap, colkey = FALSE)
points3D(xdesign[,idxAdd[1]], xdesign[,idxAdd[2]], ydesign, pch = 20, cex = 2,  col = "black", add = TRUE)
points(trans3D(x = u[,1], y = u[,2], z = pred_Knots$PhiAll.test %*% pred_Knots$xiAll.map, pmat = p),
       col = 'brown', pch = 4, lwd = 2)

p <- persp3D(x = unique(xtestGrid[, 1]), y = unique(xtestGrid[, 2]),
             z = matrix(rowMeans(sim.model$PhiAll.test %*% sim.model$xiAll.sim), ntest, ntest),
             xlab = paste("x", idxAdd[1], sep = ""), ylab = paste("x", idxAdd[2], sep = ""), zlab = "y(x1,x2)",
             main = "MCMC solution", phi = 20, theta = -30, col = colormap, colkey = FALSE)
points3D(xdesign[,idxAdd[1]], xdesign[,idxAdd[2]], ydesign, pch = 20, cex = 2,  col = "black", add = TRUE)
points(trans3D(x = u[,1], y = u[,2], z = pred_Knots$PhiAll.test %*% pred_Knots$xiAll.map, pmat = p),
       col = 'brown', pch = 4, lwd = 2)

p <- persp3D(x = unique(xtestGrid[, 1]), y = unique(xtestGrid[, 2]),
             z = matrix(sim.model$PhiAll.test %*% sim.model$xiAll.sim[,1], ntest, ntest),
             xlab = paste("x", idxAdd[1], sep = ""), ylab = paste("x", idxAdd[2], sep = ""), zlab = "y(x1,x2)",
             main = "MCMC sample", phi = 20, theta = -30, col = colormap, colkey = FALSE)
points3D(xdesign[,idxAdd[1]], xdesign[,idxAdd[2]], ydesign, pch = 20, cex = 2,  col = "black", add = TRUE)
points(trans3D(x = u[,1], y = u[,2], z = pred_Knots$PhiAll.test %*% pred_Knots$xiAll.map, pmat = p),
       col = 'brown', pch = 4, lwd = 2)

