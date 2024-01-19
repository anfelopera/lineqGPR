#  Start by importing the library we need
library("lineqGPR")
library("DiceDesign")
library("plot3D")
library("viridis")
library("ggplot2")
library("scales")

rm(list=ls())
## Functions we will study

fun1 <- function(x){ # 4 dimensional function
  return(3*(x[,1]-x[,4])^2 +sin(x[,2]+x[,3]))
}

fun2 <- function(x){
  return(10*x[,1]*x[,2] + 0.01*x[,3])
}

fun3 <- function(x){
  return(4*x[,1] + x[,2] + 10*x[,3] + 0.1*x[,4])
}

fun4 <- function(x){
  return(atan(4*x[,1] + x[,2] + 10*x[,3] + 0.1*x[,4]))
}

fun5 <- function(x) {
  return(10*x[, 1]*x[, 2] + x[, 3]*x[, 4]*x[, 6] + x[, 5]*x[, 7] + 0.1*x[,9] +sin(x[,10]+x[,11]))
}

fun6 <- function(x) {
  return(10*x[, 1]*x[, 3] + x[, 2]*x[, 4] + atan(0.3*x[5]+0.5*x[, 6]))
}

#Construction of the models

set.seed(1)
n <- c(4,3,4,4,11, 6) #corresponding dim of our function
size <- 12 #number of observations

fun <- c(6)
for (i in fun){# Making observation points
  eval(parse(text = paste("x",i," <- matrix(runif(", n[i]^2*size, ", min=0, max=1), ncol=",n[i],")",sep= "")))
  eval(parse(text = paste("x", i, " <- maximinSA_LHS(x", i , ")$design", sep="")))
  eval(parse(text = paste("y",i," <- fun",i, "(x",i,")",sep= "")))
}

constraints <- c("none", "monotonicity", "monotonicity", "monotonicity", "none", "monotonicity")
for(i in fun){
  eval(parse(text = paste("model",i," <- create(class = 'lineqBAGP', x = x",i,  ", y = y",i ,", constrType ='",
                          constraints[i],"')",sep= "")))
}

#res2 <- BAGPMaxMod(model2, max_iter = 15*3, print_iter = FALSE, nClusters =10,
#GlobalconstrType= c('boundedness',3), Estim_varnoise = TRUE, new_criteria = FALSE)

for (i in fun){# Making observation points
  message("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! FUNCTION " , i, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
  
  eval(parse(text = paste("res",i, " <- BAGPMaxMod(model",i,", max_iter = 10*",n[i],
                          ", print_iter = FALSE, nClusters = 1, GlobalconstrType= rep('",
                          constraints[i] , "',",n[i], "), Estim_varnoise = TRUE, new_criteria = FALSE)", sep=""
  )))
}


plot_history(res1,1)
plot_history(res2,2)
plot_history(res6,6)

par(mfrow = c(2,3), mar=c(1.5,1.5,1.5,1))
for (i in fun){# Making observation points
  eval(parse(text = paste("plot_history(res",i,",",i,")", sep=""
  )))
}



# res1 <- BAGPMaxMod(model1, max_iter = 10*ncol(model1$x),
#                   print_iter = FALSE, nClusters = 10, #min(12,n[1]+(n[1]*(n[1]-1))/2),
#                   GlobalconstrType= rep('boundedness',4),
#                   Block_max_size = 2, New_criteria = FALSE)
# 
# 
# res2 <- BAGPMaxMod(model2, max_iter = 10*ncol(model2$x),
#                    print_iter = FALSE,
#                    nClusters = 10,#min(12,n[2]+(n[2]*(n[2]-1))/2),
#                    tol= 1e-3,
#                    GlobalconstrType= c('monotonicity','monotonicity', 'decreasing'),
#                    Block_max_size = 3, New_criteria = FALSE)
# 

