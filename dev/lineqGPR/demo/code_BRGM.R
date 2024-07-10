library("lineqGPR")
library("DiceDesign")
library("plot3D")
library("viridis")
library("ggplot2")
library("scales")
library("reshape2") 
library("profExtrema")  # package for flooding dataset
library("RColorBrewer")
library("gridExtra")

library(tikzDevice)
library(lhs)
library(RColorBrewer)
#### loading data ####

location <- paste("/home/mderonzi/Documents/Mathis/thèse/Figures/BAGP/figs/")

rm(list=ls())
D <- 5

maxiter <- 10
#### loading data ####
data <- as.matrix(coastal_flooding)
x <- data[ , -6] # considering all the input variables
x[,3] <- cos(2*pi*x[,3])
colnames(x) <- c("Tide", "Surge", "phi", "tm", "tp")
# transforming the output in order to have small ranges.
y <- log10(data[, 6]) # area flooded
# this assumption is also made to account positivity conditions a priori.

# transforming input space into the unit square
for (k in 1:ncol(x))
  x[,k] <- (x[,k] - min(x[,k]))/max(x[,k] - min(x[,k]))

sample_size <- 80

constraints1 <- c("none","none","none","none","none","none")
constraints2 <- c("monotonicity","monotonicity","none","none","none","none")

nbTrial <- 20
Hist1 <- list()
Hist2 <- list()

# unlink("/home/mderonzi/Documents/Mathis/thèse/lineqGPR/dev/lineqGPR/Res/ResMaxModBRGM")
# Hist <- Res[[1]]
# models <- Res[[2]]

Q21 <- matrix(0, nrow = nbTrial, ncol = maxiter)
Q22 <- matrix(0, nrow = nbTrial, ncol = maxiter)

for (i in 1:nbTrial){
  message(paste("
              #######################################################################################

              ############################### New test , i = ",i,"#####################################
              
              #######################################################################################
              "))
  set.seed(10*i)
  train_i <- sample(nrow(data), sample_size)
  xtrain <- x[train_i,]
  xtest <- x[-train_i,]
  colnames(xtrain) <- c("Tide", "Surge", "phi", "tm", "tp")
  ytrain <- y[train_i]
  ytest <- y[-train_i]
  N <- length(ytest)
  
  model1 <- create(class = 'lineqBAGP', x = xtrain, y = ytrain, constrType = constraints1)
  res1 <- BAGPMaxMod(model1, max_iter = maxiter, print_iter = FALSE, nClusters = 12,
                  GlobalconstrType= constraints1, Estim_varnoise = TRUE, 
                  tolPrecision = 1e-2, tolCriteria = 1e-6, Block_max_size = 3, alpha = 1/2,
                  save_models=TRUE)
  model2 <- create(class = 'lineqBAGP', x = xtrain, y = ytrain, constrType = constraints2)
  res2 <- BAGPMaxMod(model2, max_iter = maxiter, print_iter = FALSE, nClusters = 12,
                     GlobalconstrType= constraints2, Estim_varnoise = TRUE, 
                     tolPrecision = 1e-2, tolCriteria = 1e-6, Block_max_size = 3, alpha = 1/2, beta = 1.4,
                     save_models=TRUE)
  
  Hist1[[i]] <- res1$history
  models1 <- res1$hist_models
  
  Hist2[[i]] <- res2$history
  models2 <- res2$hist_models
  for (j in 1:maxiter){
      Q21[i,j] <- 1-norm(predict(models1[[j]], xtest)$y.mod-ytest, type = "2")^2/(var(ytest)*N)
      Q22[i,j] <- 1-norm(predict(models2[[j]], xtest)$y.mod-ytest, type = "2")^2/(var(ytest)*N)
  }
}


Res <- list(Hist1, Hist2, Q21, Q22)
save(Res, file = "ResMaxModBRGM2")
#save(Res, file = "ResMaxModBRGM")
#unlink("/home/mderonzi/Documents/Mathis/thèse/lineqGPR/dev/lineqGPR/Res/ResMaxModBRGM")

#load("ResMaxModBRGM")


viridisPalette <- viridis_pal(option= "A")(maxiter)
viridisPalette[idxQ2medSort] <- viridisPalette

############################ Unconstrained case ################################
ind <- which(Q22[,10]>0.6)

Q2 <- Q21

iterSeq <- seq(maxiter)
iterSeq <- factor(iterSeq, levels = iterSeq)
Q2median <- apply(Q2, 2, median, na.rm = TRUE)
idxQ2medSort <- sort(Q2median, index.return = TRUE)$ix
Q2median <- factor(Q2median, levels = Q2median)
ggData <- data.frame(iter = rep(iterSeq, each = nbTrial),
                     Q2 = c(Q2),
                     n = factor(D+1, levels = D+1),
                     Q2med = rep(Q2median, each = nbTrial),
                     method = "No Constraint")


Q2 <- Q22

Q2median <- apply(Q2, 2, median, na.rm = TRUE)
idxQ2medSort <- sort(Q2median, index.return = TRUE)$ix
Q2median <- factor(Q2median, levels = Q2median)
ggData2 <- data.frame(iter = rep(iterSeq, each = nbTrial),
                     Q2 = c(Q2),
                     n = factor(D+1, levels = D+1),
                     Q2med = rep(Q2median, each = nbTrial),
                     method = "Constraint")

ggData3 <- rbind(ggData, ggData2)
filename <- paste(location,"ExBRGMQ2.tex", sep = "")
tikz(filename, standAlone = TRUE, width = 5, height = 4)
ggplot(ggData3, aes(x=iter, y=Q2, fill = method)) +
  geom_boxplot(width =0.5, outlier.alpha = 0) +
  geom_jitter(width=0.1, alpha=0.5, size = 0.2) +
  #scale_fill_manual(values=viridisPalette) + 
  scale_x_discrete(name ="MaxMod iteration", position = "top") + 
  theme(legend.position = "bottom")+
  coord_cartesian(ylim = c(0, 0.65))+
  theme_bw()
dev.off()
# tools::texi2dvi(filename,pdf=T,clean=TRUE)


####### CHOICES OF MAXMOD #######

nbTrial <- 20
optDecision1 <- matrix(0, nbTrial, maxiter)
optDecision2 <- matrix(0, nbTrial, maxiter)
optDecision3 <- matrix(0, nbTrial, maxiter)

for(i in 1:nbTrial){
  for(j in 1:min(length(Hist2[[i]]), maxiter)){
    if (length(Hist2[[i]][[j]])==1){
      optDecision1[i,j] <- Hist2[[i]][[j]] 
    } else{
      optDecision2[i,j] <- Hist2[[i]][[j]][1]
      optDecision3[i,j] <- Hist2[[i]][[j]][2]
    }
  }
}


#Plotting the Tab
maxiter<-10
iterdf <- rev(seq(maxiter))
iterdf <- factor(iterdf, levels = iterdf)
df <- data.frame(iter = rep(iterdf, each = nbTrial),
                 activedim = c(optDecision1[, rev(seq(maxiter))]),
                 activeblock1 = c(optDecision2[, rev(seq(maxiter))]),
                 activeblock2 = c(optDecision3[, rev(seq(maxiter))]),
                 reps = c(1:maxiter))
#bluePalette <- brewer.pal(9, "PuBu")
bluePalette <- viridis_pal(option= "D")(40)[rev(5:(5+maxiter))]
redPalette <- viridis_pal(option= "B")(50)[rev(25:(25+maxiter))]
Palette <- c(bluePalette, redPalette)
#show_col(bluePalette)
#show_col(redPalette)

fig1 <- ggplot(data=df, aes(activedim, iter)) +
  geom_count(aes(size=after_stat(..n..), color=after_stat(..n..)), shape = 15) +
  geom_count(aes(activeblock1, iter, size=after_stat(..n..)), color="red", shape = 1) +
  geom_count(aes(activeblock2, iter, size=after_stat(..n..)), color="red", shape = 1) +
  scale_x_continuous(breaks = seq(D), position = "top", limits = c(1,D)) +
  #theme(axis.text.x=element_text(angle = 45, hjust = 1)) +
  labs(title="", y = "MaxMod Iteration", x = "Input Dimension", colour = "reps", size = "reps", 
       side= "left") +
  scale_size(breaks = seq(maxiter)) +
  scale_color_continuous(breaks = seq(maxiter), low = bluePalette[1], high = rev(bluePalette)[1]) +
  guides(color = guide_legend(), size=guide_legend(), shape = guide_legend()) +
  theme_bw()
filename <- paste(location,"MaxModChoicesBRGM.tex", sep="")
tikz(filename, standAlone = TRUE, width = 5, height = 4)
plot(fig1)
dev.off()
#tools::texi2dvi(filename,pdf=T,clean=TRUE)



########################################## Comparaison Energie #########################################

maxiter <- 30
constraints <- c("monotonicity","monotonicity","none","none","none","none")

model <- create(class = 'lineqBAGP', x = x, y = y, constrType = constraints1)
res <- BAGPMaxMod(model1, max_iter = maxiter, print_iter = FALSE, nClusters = 12,
                   GlobalconstrType= constraints2, Estim_varnoise = TRUE, 
                   tolPrecision = 1e-2, tolCriteria = 1e-6, Block_max_size = 3, alpha = 1,
                   save_models=TRUE)

iteration <- 1:maxiter
energy <- res$hist_SMSE*sum((y-mean(y))^2)/sum(y^2)
nknots <- sapply(res$hist_models, function(mod) mod$nknots)


max <- max(energy)

history <- lapply(res$history, function(x) f(x))
pos <- rep(max, maxiter)
pos[4] <- 0.95*max

df <- data.frame(iter=iteration, Energy = energy)


fig <- ggplot() +
  geom_line(data = df, aes(x = nknots, y = energy), color = "blue") +
  geom_point(data = df, aes(x = nknots, y = energy), color = "blue", size = 2) + 
  geom_label(data = df, aes(x = nknots, y = pos, label = history), size=3) +
  scale_x_continuous(breaks = nknots, position = "bottom", limits = c(min(nknots), max(nknots))) +
  labs(title = "", y = "$E_n(Y,\\widehat{Y})$", x = "Number of knots") +
  theme_bw()
plot(fig)
filename <- paste(location,"MaxModIterationsBRGM.tex", sep="")
tikz(filename, standAlone = TRUE, width = 11, height = 6)
plot(fig)
dev.off()

