library("lineqGPR")

rm(list=ls())
set.seed(7)

#### Synthetic data ####
targetFun <- function(x) return(1/(1+exp(-10*(x-0.5)))) # sigmoid function
x <- seq(0, 1, by = 0.01)
y <- targetFun(x)
DoE <- splitDoE(x, y, DoE.idx = c(11, 41, 81))

#### Constrained model ####
# creating the model
model <- create(class = "lineqGP", x = DoE$xdesign, y = DoE$ydesign,
                constrType = c("boundedness", "monotonicity"), m = 20)
model$localParam$sampler <- "HMC"
model$bounds[1, ] <- c(lower = 0, upper = 1) # re-defining boundedness constraints

# simulating samples from the model
sim.model <- simulate(model, nsim = 1e3, seed = 7, xtest = DoE$xtest)

# plotting the results
layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
plot(sim.model, bounds = model$bounds[1, ], xlab = "x", ylab = "y(x)",
     main = paste("Constrained GP under ", paste(model$constrType, collapse = " & "),
                  " conditions (before covariance parameters estimaton)", sep = ""))
points(DoE$xtest, DoE$ytest, col = "red", pch = 4)

# evaluating the correlation between the MCMC samples at index "idx_acf".
idx_acf <- 21
abline(v = DoE$xtest[idx_acf], lty = 2, col = "darkgreen")
plot(sim.model$ysim[idx_acf,], type = "p", main = "MCMC samples", xlab = "samples",
     ylab = paste("y(",round(DoE$xtest[idx_acf],1),")", sep = ""), col = "darkgreen")
acf(sim.model$ysim[idx_acf,], lag.max = ncol(sim.model$ysim),
    main = "Correlation of the samples")

#### Optimizing the model ####
model2 <- lineqGPOptim(model,
                       x0 = model$kernParam$par,
                       eval_f = "logLik",
                       add.constr = FALSE,
                       mcmc.opts = list(probe = "ExpT", nb.mcmc = 1e4),
                       opts = list(algorithm = "NLOPT_LD_MMA",
                                   print_level = 3,
                                   ftol_abs = 1e-3,
                                   maxeval = 20,
                                   check_derivatives = TRUE,
                                   parfixed = c(FALSE, FALSE)),
                       lb = c(0.1, 0.1),
                       ub = c(2, 0.3))
# evaluating the model
sim.model2 <- simulate(model2, nsim = 1e3, seed = 7, xtest = DoE$xtest)

# plotting the results
layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
plot(sim.model2, bounds = model$bounds, xlab = "x", ylab = "y(x)",
     main = paste("Constrained GP under ", paste(model$constrType, collapse = " & "),
                  " conditions (after covariance parameters estimaton)", sep = ""))
points(DoE$xtest, DoE$ytest, col = "red", pch = 4)

# evaluating the correlation between the MCMC samples at index "idx_acf".
idx_acf <- 21
abline(v = DoE$xtest[idx_acf], lty = 2, col = "darkgreen")
plot(sim.model2$ysim[idx_acf,], type = "p", main = "MCMC samples", xlab = "samples",
     ylab = paste("y(",round(DoE$xtest[idx_acf],1),")", sep = ""), col = "darkgreen")
acf(sim.model2$ysim[idx_acf,], lag.max = ncol(sim.model$ysim),
    main = "Correlation of the samples")
par(mfrow = c(1,1))

