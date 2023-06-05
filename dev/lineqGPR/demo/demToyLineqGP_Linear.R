library("lineqGPR")

rm(list=ls())
set.seed(7)

#### Synthetic data ####
targetFun <- function(x) return(1/(1+exp(-10*(x-0.5)))) # sigmoid function
x <- seq(0, 1, by = 0.01)
y <- targetFun(x)
DoE <- splitDoE(x, y, DoE.idx = c(11, 41, 91))

#### Constrained model ####
# creating the model
model <- create(class = "lineqGP", x = DoE$xdesign, y = DoE$ydesign,
                constrType = "linear", m = 20)
model$localParam$sampler <- "HMC"
model$varnoise <- 1e-6
model$kernParam$par <- c(sigma2 = 1^2, theta = 0.2)

# defining the linear system of inequalities
opt <- 1 # or 2, 3, 4
m <- model$localParam$m
if (opt == 1){
  # "Forward" Monotonicity
  plotname <- "Forward Monotonicity"
  model$Lambda <- diag(m)
  diag(model$Lambda[-1,-nrow(model$Lambda)]) <- -1
  model$lb <- c(-Inf, rep(0, m-1))
  model$ub <- c(Inf, rep(0.5, m-1))
} else if (opt == 2){
  # "Backward" Monotonicity
  plotname <- "Backward Monotonicity"
  model$Lambda <- diag(m)
  diag(model$Lambda[-ncol(model$Lambda),-1]) <- -1
  model$lb <- c(rep(-Inf, m-1), -Inf)
  model$ub <- c(rep(0, m-1), Inf)
} else if (opt == 3){
  # "Short Bounded Monotonicity"
  plotname <- "Short Boundedness & Monotonicity"
  bounds <- c(lower = -1.1, upper = 1.1)
  model$Lambda <- diag(m)
  diag(model$Lambda[-1,-nrow(model$Lambda)]) <- -1
  model$Lambda <- rbind(model$Lambda, rev(model$Lambda[1,]))
  model$lb <- c(bounds[1], rep(0, m-1), bounds[1])
  model$ub <- c(bounds[2], rep(Inf, m-1), bounds[2])
} else if (opt == 4){
  # "Large Bounded Monotonicity"
  plotname <- "Large Boundedness & Monotonicity"
  bounds <- c(lower = -1.1, upper = 1.1)
  model$Lambda <- diag(m)
  diag(model$Lambda[-1,-nrow(model$Lambda)]) <- -1
  model$Lambda <- rbind(diag(m), model$Lambda)
  model$lb <- c(rep(bounds[1], m), -Inf, rep(0, m-1))
  model$ub <- c(rep(bounds[2], m), rep(Inf, m))
}

# simulating samples from the model
ptm <- proc.time()
sim.model <- simulate(model, DoE$xtest, nsim = 1e4, seed = 7)
proc.time() - ptm

# plotting the results
layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
plot(sim.model, xlab = "x", ylab = "y(x)",
     main = paste("Constrained GP under ", paste(model$constrType, collapse = " & "),
                  " conditions: ", plotname, sep = ""))

# evaluating the correlation between the MCMC samples at index "idx_acf".
idx_acf <- 21
abline(v = DoE$xtest[idx_acf], lty = 2, col = "darkgreen")
plot(sim.model$ysim[idx_acf,], type = "p", main = "MCMC samples", xlab = "samples",
     ylab = paste("y(",round(DoE$xtest[idx_acf],1),")", sep = ""), col = "darkgreen")
acf(sim.model$ysim[idx_acf,], lag.max = ncol(sim.model$ysim), main = "Correlation of the samples")
par(mfrow = c(1,1))

