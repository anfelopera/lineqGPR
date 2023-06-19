library("lineqGPR")
library("DiceDesign")
library("plot3D")
library("viridis")
# library("Matrix")

rm(list=ls())
set.seed(7)

#### Synthetic data ####
d <- 3 # number of active input variables
partition <- list(c(1,3), 2) # partition of the block structure
nblocks <- length(partition) # nb of blocks

targetFun <- function(x, partition) {
  return(x[, partition[[1]][1]]*x[, partition[[1]][2]] + x[, partition[[2]][1]])
}

# building Latin hypercube sampling (LHS) design 
nbtrain <- 4*d
xdesign <- lhsDesign(nbtrain, d, seed = 0)$design
# xdesign <- maximinSA_LHS(xdesign)$design

# adding extra training points with second block's variable = 0
nbextra <- 4
xdesign2 <- matrix(0, nbextra, d)
xdesign2[, partition[[1]]] <- maximinSA_LHS(lhsDesign(nbextra, 2, seed = 0)$design)$design
xdesign <- rbind(xdesign, xdesign2)
ydesign <- targetFun(xdesign, partition)

# defining the 3D grid for predictions 
n1D <- 20
xbase <- seq(0, 1, length = n1D)
xtest <- as.matrix(expand.grid(xbase, xbase, xbase))
ytest <- targetFun(xtest, partition)

#### Constrained model ####
# creating the model
model <- create(class = "lineqBAGP", x = xdesign, y = ydesign,
                constrType = rep("monotonicity", nblocks), 
                partition = partition,
                m = c(4, 3, 5))

# modifying the covariance parameters of each block
for (k in 1:nblocks)
  model$kernParam$par[[k]] <- c(1, rep(0.5, model$localParam$dim_block[k]))
model$localParam$sampler <- "HMC"
model$nugget <- 1e-5

# model_temp <- lineqGPOptim(model,
#                            additive = TRUE,
#                            block = TRUE,
#                            partition = list(c(1,3), 2),
#                            estim.varnoise = TRUE, # to add this info at the MaxMod level
#                            bounds.varnoise = c(1e-7, Inf), # to add this info at the MaxMod level
#                            lb = rep(1e-2, model_temp$d+1), ub = c(Inf, rep(1, model_temp$d)) # to add this info at the MaxMod level
# )

pred <- predict(model, xtest)
model.sim <- simulate(model, 1e3, seed = 1, xtest)

# Q2 criterion
var_design <- mean((ytest- mean(ytest))^2)
message("Unconstrained GP mean: ", 1 - mean((ytest- pred$y.mean)^2)/var_design)
message("Constrained GP mode: ", 1 - mean((ytest- pred$y.mode)^2)/var_design)
message("Constrained GP mean via MCMC: ", 1 - mean((ytest- rowMeans(model.sim$y.sim))^2)/var_design)


par(mfrow = c(2,2), mar=c(1.5,1.5,1.5,1))
plot(ytest, pred$y.mean, main = paste("Q2 = ", 1 - mean((ytest- pred$y.mean)^2)/var_design))
plot(ytest, pred$y.mode, main = paste("Q2 = ", 1 - mean((ytest- pred$y.mode)^2)/var_design))
plot(ytest, rowMeans(model.sim$y.sim), main = paste("Q2 = ", 1 - mean((ytest- rowMeans(model.sim$y.sim))^2)/var_design))

mean((as.matrix(ydesign) - predict(model, xdesign)$y.mean)^2)
mean((as.matrix(ydesign) - predict(model, xdesign)$y.mode)^2)

colormap <- rev(viridis(1e2))
par(mfrow = c(2,2), mar=c(1.5,1.5,1.5,1))
p <- persp3D(x = xbase, y = xbase,
             z = matrix(ytest[which(xtest[,2] == 0)], n1D, n1D),
             xlab = paste("x", partition[[1]][1], sep = ""), ylab = paste("x", partition[[1]][2], sep = ""), #zlab = "y(x1,x3)",
             main = "Target function", phi = 20, theta = -30, col = colormap, colkey = FALSE)
idxProj <- which(xdesign[,2] < 0.1)
points3D(xdesign[idxProj, partition[[1]][1]],
         xdesign[idxProj, partition[[1]][2]], ydesign[idxProj], pch = 20, cex = 2,  col = "black", add = TRUE)

p <- persp3D(x = xbase, y = xbase,
             z = matrix(pred$y.mean[which(xtest[,2] == 0)], n1D, n1D),
             xlab = paste("x", partition[[1]][1], sep = ""), ylab = paste("x", partition[[1]][2], sep = ""), # zlab = "y(x1,x3)",
             main = "Unconstrained predictive mean", phi = 20, theta = -30, col = colormap, colkey = FALSE)
points3D(xdesign[idxProj, partition[[1]][1]],
         xdesign[idxProj, partition[[1]][2]], ydesign[idxProj], pch = 20, cex = 2,  col = "black", add = TRUE)

p <- persp3D(x = xbase, y = xbase,
             z = matrix(pred$y.mode[which(xtest[,2] == 0)], n1D, n1D),
             xlab = paste("x", partition[[1]][1], sep = ""), ylab = paste("x", partition[[1]][2], sep = ""), # zlab = "y(x1,x3)",
             main = "Constrained Mode", phi = 20, theta = -30, col = colormap, colkey = FALSE)
points3D(xdesign[idxProj, partition[[1]][1]],
         xdesign[idxProj, partition[[1]][2]], ydesign[idxProj], pch = 20, cex = 2,  col = "black", add = TRUE)

p <- persp3D(x = xbase, y = xbase,
             z = matrix(rowMeans(model.sim$y.sim)[which(xtest[,2] == 0)], n1D, n1D),
             xlab = paste("x", partition[[1]][1], sep = ""), ylab = paste("x", partition[[1]][2], sep = ""), # zlab = "y(x1,x3)",
             main = "Constrained HMC Mean ", phi = 20, theta = -30, col = colormap, colkey = FALSE)
points3D(xdesign[idxProj, partition[[1]][1]],
         xdesign[idxProj, partition[[1]][2]], ydesign[idxProj], pch = 20, cex = 2,  col = "black", add = TRUE)



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