library("lineqGPR")
library("DiceDesign")
library("plot3D")
library("viridis")
library("sensitivity") #Caculer les indices de Sobol
rm(list=ls())
set.seed(7)

#### Synthetic data ####
d <- 3 # number of active input variables
partition <- list(c(1,3), 2) # partition of the block structure
nblocks <- length(partition) # nb of blocks

# targetFun <- function(x, partition)
#   return(3*x[, 1]*x[, 3] + sin(3*x[, 2]))
# 
targetFun <- function(x, partition)
   # return (log(1+3*x[, 1]+x[, 3]) + x[, 2]^2)
  return (3*x[, 1]*x[, 3]^2 + x[, 2]^2)

# sobol()
# building Latin hypercube sampling (LHS) design 
nbtrain <- 3*d
xdesign <- lhsDesign(nbtrain, d, seed = 0)$design
xdesign <- maximinSA_LHS(xdesign)$design

# adding extra training points with second block's variable = 0
nbextra <- 4
xdesign2 <- matrix(0, nbextra, d)
xdesign2[, partition[[1]]] <- maximinSA_LHS(lhsDesign(nbextra, 2, seed = 0)$design)$design
xdesign <- rbind(xdesign, xdesign2)
#xdesign <- matrix(runif(3*500, min=0, max=1), ncol=3)
ydesign <- targetFun(xdesign, partition)

# defining the 3D grid for predictions 
n1D <- 20
xbase <- seq(0, 1, length = n1D)
xtest <- as.matrix(expand.grid(xbase, xbase, xbase))
ytest <- targetFun(xtest, partition)

nblock <- length(partition) # nb of blocks
dim_block <- sapply(partition, function(x) length(x))

m <- c(2, 4, 5)
subdivision <- vector("list", nblock)
for (j in 1:nblocks)
  subdivision[[j]] <- lapply(1:dim_block[j],function(k) matrix(seq(0, 1, by = 1/(m[partition[[j]][k]]-1)), ncol = 1))
subdivision2 <- subdivision
subdivision2[[1]][[2]]<-sort(c(subdivision2[[1]][[2]], 0.333))


constrType <- c("monotonicity", "convexity", "convexity")


#### Constrained model ####
# creating the model
model <- create(class = "lineqBAGP", x = xdesign, y = ydesign,
                constrType = constrType, 
                partition = partition,
                subdivision = subdivision
                )

# model <- augment(model)
# pred <- predict(model, xtest)

# # new.model <- BAGPMaxMod(model, max_iter = 7*ncol(model$x),
# #                         reward_new_knot = 1e-4, reward_new_dim = 1e-9,
# #                         print_iter = FALSE, nClusters = 5, tol= 1e-7,
# #                         save_history = FALSE, constrType="none"
# # )[[1]]

#model <- new.model
# modifying the covariance parameters of each block
for (k in 1:nblocks)
  model$kernParam$par[[k]] <- c(1, rep(0.5, model$localParam$dim_block[k]))


# model$localParam$sampler <- "HMC"
model$nugget <- 1e-5
model <- lineqGPOptim(model,
                      additive = TRUE, # if the model is additive
                      block = TRUE, # if the model is additive per blocks
                      estim.varnoise = FALSE, # to add this info at the MaxMod level
                      bounds.varnoise = c(1e-7, Inf), # to add this info at the MaxMod level
                      lb = rep(1e-2, model$d+model$localParam$nblocks),
                      ub = c(Inf, 0.7, 0.7, Inf, 0.7), # to add this info at the MaxMod level
                      # ub = rep(Inf, model$d+model$localParam$nblocks),
                      opts = list(algorithm = "NLOPT_LD_MMA",
                      #algorithm = "NLOPT_LN_COBYLA",
                      print_level = 3,
                      ftol_abs = 1e-3,
                      maxeval = 50,
                      check_derivatives = TRUE)
)

model.sim <- simulate(model, nsim = 1e1, seed = 0, xtest = xtest)

#### Assessment of predictions ####
# Q2 criterion
var_design <- mean((ytest- mean(ytest))^2)
Q2 <- 1 - c(mean((ytest - model.sim$y.mean)^2),
            mean((ytest - model.sim$y.mode)^2),
            mean((ytest - rowMeans(model.sim$y.sim))^2)
            )/var_design
Q2 <- round(Q2*100, digits = 3)

# QQ plots
par(mfrow = c(2,2), mar=c(1.5,1.5,1.5,1))
plot(ytest, model.sim$y.mean,
     main = paste("ytest vs unconstrained predictive mean. Q2", Q2[1]),
     xlab = "ytest", ylab = "Unconstrained predictive mean")
plot(ytest, model.sim$y.mode,
     main = paste("ytest vs Constrained mode. Q2", Q2[2]),
     xlab = "ytest", ylab = "Constrained mode")
plot(ytest, rowMeans(model.sim$y.sim),
     main = paste("ytest vs Constrained HMC Mean. Q2", Q2[3]),
     xlab = "ytest", ylab = "Constrained HMC Mean")

# Plotting predictions when x2 = 0.
colormap <- rev(viridis(1e2))
par(mfrow = c(2,2), mar=c(1.5,1.5,1.5,1))
p <- persp3D(x = xbase, y = xbase,
             z = matrix(ytest[which(xtest[,2] == 0)], n1D, n1D),
             xlab = paste("x", partition[[1]][1], sep = ""),
             ylab = paste("x", partition[[1]][2], sep = ""), zlab = "Y(x1, 0, x3)",
             main = "Target function", phi = 20, theta = -30, col = colormap, colkey = FALSE)
idxProj <- which(xdesign[,2] < 0.1)
points3D(xdesign[idxProj, partition[[1]][1]],
         xdesign[idxProj, partition[[1]][2]], ydesign[idxProj], pch = 20, cex = 2,  col = "black", add = TRUE)

p <- persp3D(x = xbase, y = xbase,
             z = matrix(model.sim$y.mean[which(xtest[,2] == 0)], n1D, n1D),
             xlab = paste("x", partition[[1]][1], sep = ""),
             ylab = paste("x", partition[[1]][2], sep = ""), zlab = "Y(x1, 0, x3)",
             main = paste("Unconstrained predictive mean. Q:", Q2[1]),
             phi = 20, theta = -30, col = colormap, colkey = FALSE)
points3D(xdesign[idxProj, partition[[1]][1]],
         xdesign[idxProj, partition[[1]][2]], ydesign[idxProj], pch = 20, cex = 2,  col = "black", add = TRUE)

p <- persp3D(x = xbase, y = xbase,
             z = matrix(model.sim$y.mode[which(xtest[,2] == 0)], n1D, n1D),
             xlab = paste("x", partition[[1]][1], sep = ""),
             ylab = paste("x", partition[[1]][2], sep = ""), zlab = "Y(x1, 0, x3)",
             main = paste("Constrained mode. Q:", Q2[2]),
             phi = 20, theta = -30, col = colormap, colkey = FALSE)
points3D(xdesign[idxProj, partition[[1]][1]],
         xdesign[idxProj, partition[[1]][2]], ydesign[idxProj], pch = 20, cex = 2,  col = "black", add = TRUE)

p <- persp3D(x = xbase, y = xbase,
             z = matrix(rowMeans(model.sim$y.sim)[which(xtest[,2] == 0)], n1D, n1D),
             xlab = paste("x", partition[[1]][1], sep = ""),
             ylab = paste("x", partition[[1]][2], sep = ""), zlab = "Y(x1, 0, x3)",
             main = paste("Constrained HMC Mean. Q:", Q2[3]),
             phi = 20, theta = -30, col = colormap, colkey = FALSE)
points3D(xdesign[idxProj, partition[[1]][1]],
         xdesign[idxProj, partition[[1]][2]], ydesign[idxProj], pch = 20, cex = 2,  col = "black", add = TRUE)
