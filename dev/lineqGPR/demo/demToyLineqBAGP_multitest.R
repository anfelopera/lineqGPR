#  Start by importing the library we need
library("lineqGPR")
library("DiceDesign")
library("plot3D")
library("viridis")
library("ggplot2")
library("scales")
library("reshape2") 

rm(list=ls())
## Functions we will study

fun1 <- function(x){ # 4 dimensional function
  return(log(0.01+3*x[,1]+2*x[,4]) +((x[,2]-0.5)+(x[,3]-0.5))^2)
}

fun2 <- function(x){
  return(10*x[,1]*x[,2] + 0.5*sin(x[,3]))
}

fun3 <- function(x){
  return(4*x[,1] + x[,2] + 10*x[,3] + 0.1*x[,4])
}

fun4 <- function(x){
  return(atan(2*x[,1] + 0.5*x[,2]) + atan(x[,3] + 0.1*x[,4]))
}

fun5 <- function(x) {
  return(10*x[, 1]*x[, 2] + x[, 3]*x[, 4]*x[, 6] + x[, 5]*x[, 7] + 0.1*x[,9] +sin(0.2*x[,10]+0.2*x[,11]))
}

fun6 <- function(x) {
  return(10*x[, 1]*x[, 3] + sin(x[, 2]*x[, 4]) + atan(3*x[,5]+5*x[, 6]))
}

#Construction of the models

n <- c(4,3,4,4,11, 6) #corresponding dim of our function
multiplier <- c(13, 5, 4, 13, 12, 8) #number of observations per dimension
sizes <- n*multiplier

seeds <- c(0,0,2,0,11,0) #The random seed used
fun <- c(1,2,3,4,5,6)
for (i in fun){# Making observation points
  eval(parse(text = paste("set.seed(",seeds[i],")",sep= "")))
  eval(parse(text = paste("x",i," <- matrix(runif(", n[i]*sizes[i], ", min=0, max=1), ncol=",n[i],")",sep= "")))
  eval(parse(text = paste("x", i, " <- maximinSA_LHS(x", i , ")$design", sep="")))
  eval(parse(text = paste("y",i," <- fun",i, "(x",i,")",sep= "")))
  eval(parse(text = paste("xtest",i," <- matrix(runif(", 10000*n[i], ", min=0, max=1), ncol=",n[i],")",sep= "")))
}

constraints <- c("none", "monotonicity", "monotonicity", "monotonicity", "monotonicity", "monotonicity")
const <- c( list(c("monotonicity", "convexity", "convexity", "monotonicity")), 
               lapply(2:6, function(i) rep(constraints[i], n[i])))
  
for(i in fun){
  eval(parse(text = paste("model",i," <- create(class = 'lineqBAGP', x = x",i,  ", y = y",i ,", constrType ='",
                          constraints[i],"')",sep= "")))
}

tolPrec <- c(1e-4, 1e-7, 1e-7, 5e-5, 1e-7, 1e-7)
tolCriteria <- c(1e-6, 1e-7, 1e-7, 1e-6, 1e-7, 5e-7)
for (i in fun){# Making observation points
  message("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! FUNCTION " , i, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
  eval(parse(text = paste("res",i, " <- BAGPMaxMod(model",i,", max_iter = (7*",n[i],
                          "), print_iter = FALSE, nClusters = 12, GlobalconstrType= const[[",
                          i, "]], Estim_varnoise = TRUE, tolPrecision = ",tolPrec[i], 
                          ", tolCriteria = ", tolCriteria[i]  ,", Block_max_size = 3)", sep=""
  )))
  eval(parse(text = paste("model",i," <- res", i,  "[[1]]", sep= "")))
}

plot_history(res1)
plot_history(res2)
plot_history(res3)
plot_history(res4)
plot_history(res5)
plot_history(res6)

partitions <- list(list(c(1,4),c(2,3)), list(c(1,2), c(3)), list(c(1), c(2), c(3), c(4)),
                   list(c(1,2), c(3,4)), list(c(1,2),c(3,4,6), c(5,7), c(9), c(10,11)),
                   list(c(1,3),c(2,4), c(5,6)))

nblocks <- c(2,2,4,2,5,3)


f11 <- function(x){ # 4 dimensional function
  return(log(0.01+3*x[,1]+2*x[,4])-0.8035706 )
}
f12 <- function(x){ # 4 dimensional function
  return(((x[,2]-0.5)+(x[,3]-0.5))^2-0.1666666)
}

f21 <- function(x){
  return(10*(x[,1]*x[,2]- 1/4))
}
f22 <- function(x){
  return(0.5*(sin(x[,3])-0.4596976))
}

f31 <- function(x){
  return(4*x[,1]-2) 
}
f32 <- function(x){
  return(x[,2]-0.5)
} 
f33 <- function(x){
  return(10*(x[,3]-0.5))
}
f34 <- function(x){
  return(0.1*(x[,4]-0.5))
}

f41 <- function(x){
  return(atan(2*x[,1] + 0.5*x[,2]) - 0.8266132)
}
f42 <- function(x){
  return(atan(x[,3] + 0.1*x[,4])-0.4772659)
}

f51 <- function(x) {
  return(10*x[, 1]*x[, 2]-2.5)
}
f52 <- function(x){
  return(x[, 3]*x[, 4]*x[, 6]-1/8)
}
f53 <- function(x){
  return(x[, 5]*x[, 7]- 1/4)
}
f54 <- function(x){
  return(0.1*(x[, 9]-1/2))
}
f55 <- function(x){
  return(sin(0.2*x[,10]+0.2*x[,11])-0.1980067)
}

f61 <- function(x) {
  return(10*x[, 1]*x[, 3]-2.5 )
}
f62 <- function(x) {
  return(sin(x[, 2]*x[, 4])-0.24)
}
f63 <- function(x) {
  return(atan(3*x[,5]+5*x[, 6])-1.265833)
}

# N <- 1000
# fun <- c(2)#1,2,3,4,5,6)
# for (k in fun){
#   for (j in 1:length(partitions[[k]])){
#     eval(parse(text = paste("s",k,j, " <- 0", sep= "")))
#     for (i in 1:N){
#       eval(parse(text = paste("s",k,j," <- s",k,j, " + mean(f",k,j,
#     "(matrix(runif(100000*", n[k], ", min=0, max=1), ncol=",
#     n[k],")))",sep= "")))
#     }
#     eval(parse(text = paste("s",k,j," <- s",k,j, "/N",sep= "")))
#   }
# }

xplot <- 0:100/100
yplot <- 0:100/100
xyplot <- as.matrix(expand.grid(xplot,yplot))
fun <- c(1,2,3,4,6)
for (k in fun){
  eval(parse(text = paste("b_fun", k," <- block_fun(model", k,  ")", sep= "")))
  xy <- matrix(0, nrow=10201, ncol=, n[k])
  for (j in 1:length(partitions[[k]])){
    if (length(partitions[[k]][[j]])==2){
      xy[,partitions[[k]][[j]]] <- xyplot
    } else if (length(partitions[[k]][[j]])==1){
      xy[,partitions[[k]][[j]]] <- xyplot[,1]
    }
    eval(parse(text = paste("z",k,j," <-  f", k, j, "(xy)", sep= "")))
    eval(parse(text = paste("y",k,j," <-  f", k, j, "(x", k, ")", sep= "")))
  }
  eval(parse(text = paste("z",k," <- Compute_fun(model",k,", b_fun",k,", xy)", sep= "")))
}

################################# Function 1 ########################################
par(mfrow = c(1,2), mar=c(1,1,1,1))

persp3D(xplot, yplot, matrix(z11, ncol=101),
        xlab="x", ylab="y", zlab="z", theta = 30, phi = 10, expand = 0.5,
        main = "Target function f11")
points3D(x1[, partitions[[1]][[1]][1]], x1[, partitions[[1]][[1]][2]],
         y11, pch = 20, cex = 2,  col = "black", add = TRUE)
persp3D(xplot, yplot, matrix(z1[[2]], ncol=101),
        xlab="x", ylab="y", zlab="z", theta = 30, phi = 10, expand = 0.5,
        main = "Predictor function of f11")
points3D(x1[, partitions[[1]][[1]][1]], x1[, partitions[[1]][[1]][2]],
         y11, pch = 20, cex = 2,  col = "black", add = TRUE)


persp3D(xplot, yplot, matrix(z12, ncol=101),
        xlab="x", ylab="y", zlab="y", theta = 30, phi = 10, expand = 0.5,
        main = "Target function f11")
points3D(x1[, partitions[[1]][[2]][1]], x1[, partitions[[1]][[2]][2]],
         y12, pch = 20, cex = 2,  col = "black", add = TRUE)
persp3D(xplot, yplot, matrix(z1[[1]], ncol=101),
        xlab="x", ylab="y", zlab="y", theta = 30, phi = 10, expand = 0.5,
        main = "Predictor function of f12")
points3D(x1[, partitions[[1]][[2]][1]], x1[, partitions[[1]][[2]][2]],
         y12, pch = 20, cex = 2,  col = "black", add = TRUE)

#plot(xplot,matrix(z1[[1]], ncol=101)[100,] )
############################ Function 2 ###########################################

par(mfrow = c(1,2))#, mar=c(1,1,1,1))

persp3D(xplot, yplot, matrix(z21, ncol=101),
        xlab="x", ylab="y", zlab="z", theta = 30, phi = 10, expand = 0.5,
        main = "Target function f21")
points3D(x2[, partitions[[2]][[1]][1]], x2[, partitions[[2]][[1]][2]],
         y21, pch = 20, cex = 2,  col = "black", add = TRUE)
persp3D(xplot, yplot, matrix(z2[[1]], ncol=101),
        xlab="x", ylab="y", zlab="z", theta = 30, phi = 10, expand = 0.5,
        main = "Predictor function of f21")
points3D(x2[, partitions[[2]][[1]][1]], x2[, partitions[[2]][[1]][2]],
         y21, pch = 20, cex = 2,  col = "black", add = TRUE)


plot(x2[, partitions[[2]][[2]][1]], y22, pch=19, main = "Target function f22")
lines(xplot, matrix(z22, ncol=101)[,1], col="red")

plot(x2[, partitions[[2]][[2]][1]], y22, pch=19, main = "Predictor function of f22")
lines(xplot, matrix(z2[[2]], ncol=101)[,1], col="blue")



######################## Function 4 #####################################

par(mfrow = c(1,2), mar=c(1,1,1,1))

persp3D(xplot, yplot, matrix(z41, ncol=101),
        xlab="x", ylab="y", zlab="z", theta = 30, phi = 10, expand = 0.5,
        main = "Target function f41")
points3D(x4[, partitions[[4]][[1]][1]], x4[, partitions[[4]][[1]][2]],
         y41, pch = 20, cex = 2,  col = "black", add = TRUE)
persp3D(xplot, yplot, matrix(z4[[2]], ncol=101),
        xlab="x", ylab="y", zlab="z", theta = 30, phi = 10, expand = 0.5,
        main = "Predictor function of f41")
points3D(x4[, partitions[[4]][[1]][1]], x4[, partitions[[4]][[1]][2]],
         y41, pch = 20, cex = 2,  col = "black", add = TRUE)


persp3D(xplot, yplot, matrix(z42, ncol=101),
        xlab="x", ylab="y", zlab="y", theta = 0, phi = 10, expand = 0.5,
        main = "Target function f42")
points3D(x4[, partitions[[4]][[2]][1]], x4[, partitions[[4]][[2]][2]],
         y42, pch = 20, cex = 2,  col = "black", add = TRUE)
persp3D(xplot, yplot, matrix(z4[[1]], ncol=101),
        xlab="x", ylab="y", zlab="y", theta = 0, phi = 10, expand = 0.5,
        main = "Predictor function of f42")
points3D(x4[, partitions[[4]][[2]][1]], x1[, partitions[[4]][[2]][2]],
         y42, pch = 20, cex = 2,  col = "black", add = TRUE)


################################## Function 6 ########################################

par(mfrow = c(1,2), mar=c(1,1,1,1))

persp3D(xplot, yplot, matrix(z61, ncol=101),
        xlab="x", ylab="y", zlab="z", theta = 60, phi = 15, expand = 0.5,
        main = "Target function f61")
points3D(x6[, partitions[[6]][[1]][1]], x6[, partitions[[6]][[1]][2]],
         y61, pch = 20, cex = 2,  col = "black", add = TRUE)
persp3D(xplot, yplot, matrix((z6[[3]]), ncol=101),
        xlab="x", ylab="y", zlab="z", theta = 60, phi = 15, expand = 0.5,
        main = "Predictor function of f61")
points3D(x6[, partitions[[6]][[1]][1]], x6[, partitions[[6]][[1]][2]],
         y61, pch = 20, cex = 2,  col = "black", add = TRUE)

persp3D(xplot, yplot, matrix(z62, ncol=101),
        xlab="x", ylab="y", zlab="z", theta = 30, phi = 10, expand = 0.5,
        main = "Target function f62")
points3D(x6[, partitions[[6]][[2]][1]], x6[, partitions[[6]][[2]][2]],
         y62, pch = 20, cex = 2,  col = "black", add = TRUE)
persp3D(xplot, yplot, matrix(z6[[2]], ncol=101),
        xlab="x", ylab="y", zlab="z", theta = 30, phi = 10, expand = 0.5,
        main = "Predictor function of f62")
points3D(x6[, partitions[[6]][[2]][1]], x6[, partitions[[6]][[2]][2]],
         y62, pch = 20, cex = 2,  col = "black", add = TRUE)

persp3D(xplot, yplot, matrix(z63, ncol=101),
        xlab="x", ylab="y", zlab="z", theta = 30, phi = 10, expand = 0.5,
        main = "Target function f63")
points3D(x6[, partitions[[6]][[3]][1]], x6[, partitions[[6]][[3]][2]],
         y63, pch = 20, cex = 2,  col = "black", add = TRUE)
persp3D(xplot, yplot, matrix(z6[[1]], ncol=101),
        xlab="x", ylab="y", zlab="z", theta = 30, phi = 10, expand = 0.5,
        main = "Predictor function of f63")
points3D(x6[, partitions[[6]][[3]][1]], x6[, partitions[[6]][[3]][2]],
         y63, pch = 20, cex = 2,  col = "black", add = TRUE)

