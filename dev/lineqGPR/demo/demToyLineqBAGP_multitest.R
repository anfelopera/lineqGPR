#  Start by importing the library we need

library("lineqGPR")
library("DiceDesign")
library("plot3D")
library("viridis")


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

fun5 <- function(x, partition) {
  return(10*x[, 1]*x[, 2] + x[, 3]*x[, 4]*x[, 6] + x[, 5]*x[, 7] + 0.1*x[,9] +sin(x[,10]+x[,11]))
}


#Construction of the models

set.seed(10)
n <- c(4,3,4,4,11) #corresponding dim of our function
size <- 40 #number of observations
for (i in 1:5){# Making observation points
  #  eval(parse(text = paste("x", i, " <- maximinSA_LHS(x", i , ")$design", sep="")))
  eval(parse(text = paste("x",i," <- matrix(runif(", n[i]^2*size, ", min=0, max=1), ncol=",n[i],")",sep= "")))
  eval(parse(text = paste("y",i," <- fun",i, "(x",i,")",sep= "")))
}

constraints <- c("none", "monotonicity", "monotonicity", "monotonicity", "boundedness")
for(i in 1:5){
  eval(parse(text = paste("model",i," <- create(class = 'lineqBAGP', x = x",i,  ", y = y",i ,", constrType ='",
                          constraints[i],"')",sep= "")))
}

#res2 <- BAGPMaxMod(model2, max_iter = 15*3, print_iter = FALSE, nClusters =10,
#GlobalconstrType= c('boundedness',3), Estim_varnoise = TRUE, new_criteria = FALSE)
fun <- c(2,3,4)
for (i in fun){# Making observation points
  message("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! FUNCTION " , i, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
  
  eval(parse(text = paste("res",i, " <- BAGPMaxMod(model",i,", max_iter = 15*",n[i],
                          ", print_iter = FALSE, nClusters = 1, GlobalconstrType= rep('",
                          constraints[i] , "',",n[i], "), Estim_varnoise = TRUE, new_criteria = FALSE)", sep=""
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



# res1[[2]]
# 
# 
# 
# res2[[2]]
# 
# 
# res3[[2]]
# 
# 
# res4[[2]]
# 
# 
# res5[[2]]
# 
# 
