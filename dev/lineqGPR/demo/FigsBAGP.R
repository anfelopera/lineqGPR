library("lineqGPR")
library("DiceDesign")
library("plot3D")
library("viridis")
library("ggplot2")
library("scales")
library("reshape2") 

rm(list=ls())
#####  Size block 2 ######

dim <- c(5, 10, 20, 40, 60)#, 100)
nknot <- 5
nobs <- 3
for (d in dim){
  #Creation of the functions
  D<-2*d
  eval(parse(text = paste(
    "fun", D, " <- function(x){
    sum <- ", paste("atan(", (5*(1-(1:d)/(d+1))), "*(x[,", 2*(1:d)-1, "]+ 2*x[,", 2*(1:d), "]))",
                    collapse="+"),"\n",
    "return(sum)
    }"
    , sep="")))
  #Creation of the partitions
  eval(parse(text = paste(
    "partition", D, " <- list(",
    paste("c(", 2*(1:d)-1, ",", 2*(1:d), ")",
          collapse=","),
    ")",
    sep = "")))
  #Creation of the subdivisions
  eval(parse(text = paste(
    "subdivision", D, " <- lapply(1:d, function(i) list((0:nknot)/nknot,(0:nknot)/nknot))",
    sep = "")))
}

# print("Creation Models")
# N<-c(1000, 100, 50, 20, 10)
# KernParam <- lapply(1:5, function(i) lapply(1:dim[i], function(x) c(0,0,0)))
# varnoise <- rep(0, 5)
# for(i in 1:5){
#   print(paste("######### Model", i," ##########", sep="")) 
#   D<-2*dim[i]
#   for (n in 1:N[i]){
#     x <- lhsDesign(10*dim[i], D, seed=n)$design
#     eval(parse(text = paste("y <- fun", D, "(x)",sep= "")))
#     eval(parse(text = paste("model <- create(class = 'lineqBAGP', x = x, y = y, partition = partition", D, 
#                             ", subdivision = subdivision", D, ", constrType = rep('monotonicity',", D, "))",
#                             sep= "")))
#     print(paste("Optimisation", n, sep= ""))
#     model <- optimise.parameters(model, model$varnoise, TRUE)
#     for(j in 1:dim[i]){
#       KernParam[[i]][[j]]<-KernParam[[i]][[j]]+model$kernParam$par[[j]]/N[i]
#       varnoise[i]<-varnoise[i]+model$varnoise/N[i]
#     }
#   }
# }

#Res <- list(KernParam, varnoise)
#save(Res, file = "ResBAGP")
load("ResBAGP")
KernParam <- Res[[1]]
varnoise <- c(1e-3,1e-3,1e-2,1e-2, 1e-1)
#unlink("/home/mderonzi/Documents/Mathis/thèse/lineqGPR/dev/lineqGPR/Res/ResBAGP")
print("Fin de l'évaluation des paramètres! ")

for (i in 1:5){
  D<-2*dim[i]
  N <- 10^3
  ########################################### Values of interest ######################################
  eval(parse(text = paste("Q2_sim", i, " <- c()", sep="")))
  eval(parse(text = paste("Q2_c", i, " <- c()", sep="")))
  eval(parse(text = paste("Q2_m", i, " <- c()", sep="")))
  eval(parse(text = paste("time_m", i, " <- c()", sep="")))
  eval(parse(text = paste("time_s", i, " <- c()", sep="")))
  ######################################## Design of experiment Test #################################
  eval(parse(text = paste("xtest",i," <- lhsDesign(", N,",", D, ", seed = 0)$design", sep="")))
  eval(parse(text = paste("ytest",i," <- fun", D, "(xtest",i,")", sep="")))
}

for (j in 1:10){
        print("################################### Start of Computation of  Q2 ########################################################")
  print(paste("#####################################      J = ", j,"     ###############################################################"))
  for(i in 1:5){
    N <- 10^3
    print(paste("######## Prediction", i, "############", sep=""))
    D <- 2*dim[i]
    #################################################### DOE Creation ################3#######################################
    eval(parse(text = paste("x", D, " <- lhsDesign(", nobs*D ,",", D, ", seed = ", j, ")$design", sep="")))
    #eval(parse(text = paste("x", D, " <- maximinSA_LHS(x", D , ")$design", sep="")))
    eval(parse(text = paste("y", D," <- fun", D, "(x",D,")", sep= "")))
    ################################################### Model creation #######################################################
    eval(parse(text = paste("model",i," <- create(class = 'lineqBAGP', x = x", D,  ", y = y",D ,", partition = partition", D, 
                            ", subdivision = subdivision", D, ", constrType = rep('monotonicity',", D, "))",sep= "")))
    eval(parse(text = paste("model",i,"$kernParam$par <- KernParam[[i]]", sep= "")))
    eval(parse(text = paste("model",i,"$varnoise <- varnoise[i]", sep= "")))
    eval(parse(text = paste("model",i,"$localParam$sampler <- 'HMC'", sep="")))
    #################################################  Computation of results #####################################################
    
    print(paste("#### begin sim ", i, "#####", sep=""))
    eval(parse(text = paste("sim",i," <-simulate(model", i,", nsim = ", max(10^(5-i),10) ,", seed = 0, xtest = xtest",i,")", sep="")))
    eval(parse(text = paste("Q2_sim",i," <-c(Q2_sim",i, ",1-norm(as.matrix(rowMeans(sim",i,"$y.sim)-ytest",i,
                            "), type = '2')^2/(N*var(ytest",i,")))", sep="")))
    print("end sim")
    print("Computation of mean and mode")
    eval(parse(text = paste("Q2_c",i," <-c(Q2_c",i, ",1-norm(sim",i,"$y.mode-ytest",i,
                            ", type = '2')^2/(N*var(ytest",i,")))", sep="")))
    eval(parse(text = paste("Q2_m",i," <-c(Q2_m",i, ",1-norm(sim",i,"$y.mean-ytest",i,
                            ", type = '2')^2/(N*var(ytest",i,")))", sep="")))
    eval(parse(text = paste("time_m",i," <-c(time_m",i, ",sim",i,"$predtime)", sep="")))
    eval(parse(text = paste("time_s",i," <-c(time_s",i, ",sim",i,"$simtime)", sep="")))
  }
}

for (i in 1:5){
  eval(parse(text = paste("RQ2",i," <-list(c(mean(Q2_c",i, "), sqrt(var(Q2_c",i, "))), ",
                          "c(mean(Q2_m",i,"), sqrt(var(Q2_m",i,"))), ",
                          "c(mean(Q2_sim",i,"), sqrt(var(Q2_sim",i,"))))",
                          sep="")))
  eval(parse(text = paste("TQ2",i," <-list(c(mean(time_m", i, "[(5*(0:9)+3)]), sqrt(var(time_m", i,"[(5*(0:9)+3)]))),",
                          "c(mean(time_s", i, "[(5*(0:9)+3)]), sqrt(var(time_s", i,"[(5*(0:9)+3)]))))",
                          sep="")))
}
