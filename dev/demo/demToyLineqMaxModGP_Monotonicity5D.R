library("lineqGPR")
require("DiceDesign")
library("plot3D")
library("viridis")

rm(list=ls())

# Synthetic data: monotonic function
modatan <- function(x, a) return(atan(a*x))
targetFun <- function(x, d) {
  y <- 0
  a <- (1-(1:d)/(d+1))*5
  for (k in 1:d)
    y <- y + modatan(x[, k], a[k])
  return(y)
}

D <- 5 # number of active + inactive input variables
d <- 2 # number of active input variables

# building DoE in dimension D via Latin hypercube sampling (LHS)
nbtrain <- 20*D
xdesign <- lhsDesign(nbtrain, D, seed = 8)$design
xdesign <- maximinSA_LHS(xdesign)$design
ydesign <- targetFun(xdesign, d)

# building a DoE for assessing the model
ntest <- 80*D
xtest <- lhsDesign(ntest, D, seed = 8)$design
xtest <- maximinSA_LHS(xtest)$design

# initializing a 1D GP model with only two knots 
model <-  create(class = 'lineqGP',
                 x = xdesign, y = ydesign,
                 constrType = c("monotonicity"), m = 2)
model$kernParam$type <- "matern52"
model$varnoise <- var(ydesign)
model$kernParam$nugget <- 1e-5

model <- MaxMod(model,
                xtest,
                tol = 1e-5,
                max_iter = 10*model$d,
                reward_new_knot = 1e-6,
                reward_new_dim = 1e-9,
                print_iter = TRUE,
                nClusters = 10,
                save_history = TRUE)

message("\nNumber of active dimensions: ", d)
message("Number of actived dimensions via MaxMod: ", model$d, "\n")

# evaluating the model using an equispaced grid of points
ntest <- 10
xtestGrid  <- as.matrix(expand.grid(seq(0, 1, length = ntest), seq(0, 1, length = ntest)))
ytestGrid <- targetFun(xtestGrid, d)
pred <- predict(model, xtestGrid)

# plotting the MAP estimate
colormap <- rev(viridis(1e2))
par(mfrow = c(1,2))
p <- persp3D(x = seq(0, 1, length = ntest), y = seq(0, 1, length = ntest),
             z = matrix(pred$Phi.test %*% pred$xi.map, nrow = ntest),
             xlab = "x2", ylab = "x1", zlab = "y(x1,x2)",
             main = "target function",
             phi = 20, theta = -30, col = colormap,
             contour = TRUE, colkey=FALSE)
points(trans3D(x = model$x[, 1], y = model$x[, 2], z = ydesign, pmat = p),
       col = 'black', pch = 19)
u <- expand.grid(model$ulist[[1]], model$ulist[[2]])
pred_Knots <- predict(model, as.matrix(u))
points(trans3D(x = u[, 1], y = u[, 2], z = pred_Knots$Phi.test %*% pred_Knots$xi.map, pmat = p),
       col = 'brown', pch = 4, lwd = 2)

diff = ytestGrid - pred$Phi.test %*% pred$xi.map
image2D(abs(matrix(diff, nrow = ntest)), col = rev(colormap),
        main = "absolute error", xlab = "x2", ylab = "x1")
points2D(model$x[, 1], model$x[, 2], add = TRUE, pch = 19, col ='black')
points2D(u[, 1], u[, 2], add = TRUE, pch = 4, lwd = 4, col ='brown')

