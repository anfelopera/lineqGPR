library("lineqGPR")
library("DiceDesign")
library("plot3D")
library("viridis")
# library("Matrix")

rm(list=ls())
set.seed(7)
#### Synthetic data ####
d <- 8 # number of active input variables
partition1 <- list(7,8,c(1,4), c(2,3))
partition2 <- list(c(2,3,7,8), c(1,4,5))
partition3 <- list(c(1,2), c(3,4,6), c(5,7))
subdivision1 <- list(list(c(0,0.3,1)),list(c(0,0.5,1)), list(c(0,0.4,1),c(0,1)), list(c(0,0.1,1),c(0,1)))
subdivision2 <- list(list(c(0,0.1,0.3,1),c(0,0.3,0.7,1),c(0,0.3,0.4,0.6,0.7,1),c(0,0.5,1)),
                     list(c(0,0.1,0.4,1), c(0,0.3,0.7,1),c(0,1)))
subdivision3 <- list(list(c(0,1),c(0,1)),list(c(0,1),c(0,1),c(0,1)),
                     list(c(0,1), c(0,1)))
partition <- list(2,3,c(7,8), c(1,4))
subdivision <- list(list(c(0,0.3,1)),list(c(0,1)), list(c(0,0.4,1),c(0,1)), list(c(0,0.1,1),c(0,1)))

#sol <- merge_block(partition1,subdivision1,2,3)
condition_changment_basis(partition1,partition2,subdivision1,subdivision2)
#condition_changment_basis(partition1,sol[[1]],subdivision1,sol[[2]])

nblock1 <- length(partition1) #number of blocks in partition1 
nblock2 <- length(partition2) #number of blocks in partition2

dim_block1 <- sapply(partition1, function(x) length(x))
dim_block2 <- sapply(partition2, function(x) length(x))

targetFun <- function(x, partition) {
  return(10*x[, 1]*x[, 2] + 3*x[, 3]*x[, 4]*x[, 6] + x[, 5]*x[, 7])
}

d <- 8 
# building Latin hypercube sampling (LHS) design
#set.seed(1)
#xdesign <- matrix(runif(8*100, min=0, max=1), ncol=8)

nbtrain <- 20*d
xdesign <- lhsDesign(nbtrain, d, seed = 0)$design
xdesign <- maximinSA_LHS(xdesign)$design
ydesign <- targetFun(xdesign,partition2)

ntest <- 300
xtest <- lhsDesign(ntest, d, seed = 0)$design
ytest <- targetFun(xtest,partition2)

# xtest <- maximinSA_LHS(xtest)$design


#### Constrained model ####
# creating the model

model1 <- create(class = "lineqBAGP", x = xdesign, y = ydesign,
                constrType = rep("monotonicity", ncol(xdesign)), 
                partition = partition1,
                subdivision =subdivision1
)
model2 <- create(class = "lineqBAGP", x = xdesign, y = ydesign,
                 constrType = rep("monotonicity", nblock2), 
                 partition = partition2,
                 subdivision =subdivision2
)
model3 <- create(class = "lineqBAGP", x = xdesign, y = ydesign,
                 constrType = rep("monotonicity", length(partition3)), 
                 partition = partition3,
                 subdivision =subdivision3
)
new.model1 <- optimise.parameters(model1)
merged <- merge_block(new.model1$partition, new.model1$subdivision, new.model1$kernParam, 2, 3)

#new.model3 <- optimise.parameters(model3)
#new.model3$kernParam <- model3$kernParam
#new.model3$varnoise <- 1e-5 #model3$varnoise
#pred1 <- predict(model3, xtest)$y.mean
#pred2 <- predict(new.model3, xtest)$y.mean

#mean((pred1-ytest)^2)
#mean((pred2-ytest)^2)
###############################################################################



#criteria <- square_norm_int(model1,new.model1, same_basis = TRUE)

new.model <- BAGPMaxMod(model1, max_iter = 10*ncol(model1$x),
                       reward_new_knot = 1e-4, reward_new_dim = 1e-9,
                       print_iter = FALSE, nClusters = 10,
                       GlobalconstrType= rep("monotonicity", ncol(model1$x))
                       )





# model1 <- lineqGPOptim(model1,additive = TRUE,
#                    block = TRUE,
#                    estim.varnoise = TRUE, # to add this info at the MaxMod level
#                    bounds.varnoise = c(1e-7, Inf), # to add this info at the MaxMod level
#                    lb = rep(1e-2, model1$d+model1$localParam$nblocks),
#                    ub = c(Inf, 0.7, Inf, 0.7, Inf, 0.7, 0.7, Inf, 0.7, 0.7), # to add this info at the MaxMod level
#                                    # ub = rep(Inf, model$d+model$localParam$nblocks),
#                                    opts = list(algorithm = "NLOPT_LD_MMA",
#                                    #algorithm = "NLOPT_LN_COBYLA",
#                                    print_level = 3,
#                                    ftol_abs = 1e-3,
#                                    maxeval = 50,
#                                    check_derivatives = TRUE)
#              )
# 
# model2 <- lineqGPOptim(model2,additive = TRUE,
#              block = TRUE,
#              estim.varnoise = TRUE, # to add this info at the MaxMod level
#              bounds.varnoise = c(1e-7, Inf), # to add this info at the MaxMod level
#              lb = rep(1e-2, model2$d+model2$localParam$nblocks),
#              ub = c(Inf, 0.7, 0.7, 0.7, 0.7, Inf, 0.7, 0.7, 0.7), # to add this info at the MaxMod level
#              # ub = rep(Inf, model$d+model$localParam$nblocks),
#              opts = list(algorithm = "NLOPT_LD_MMA",
#                 #algorithm = "NLOPT_LN_COBYLA",
#                          print_level = 3,
#                          ftol_abs = 1e-3,
#                          maxeval = 50,
#                          check_derivatives = TRUE)
# )

#pred1 <- predict(model1, xtest)
#pred2 <- predict(model2, xtest)

# plot(ytest,pred2$y.mode)


############### Testes pour vérifier que le calculs sont bons : square norm ###########


partition0 <- list(c(1,2))
partition1 <- list(c(1,2),5)
partition2 <- list(c(1,2,5))

subdivision0 <- list(list(c(0,1),c(0,1)))
subdivision1 <- list(list(c(0,1),c(0,1)),list(c(0,1)))
subdivision2 <- list(list(c(0,1),c(0,1),c(0,1)))

model0 <- create(class = "lineqBAGP", x = xdesign, y = ydesign,
                 constrType = rep("monotonicity", 1),
                 partition = partition0,
                 subdivision =subdivision0
)
x <- matrix(runif(8*100000, min=0, max=1), ncol=8)
pred0 <- predict(model0,x)
pred1 <- predict(model1,x)
pred2 <- predict(model2,x)

sum((pred0$y.mode-pred1$y.mode)**2)/100000
square_norm_int(model0,model1)
sum((pred0$y.mode-pred2$y.mode)**2)/100000
square_norm_int(model0,model1)
square_norm_int(model0,model2)

str(new.model3)


