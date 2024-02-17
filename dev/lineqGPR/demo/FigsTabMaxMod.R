library("lineqGPR")
library("DiceDesign")
library("plot3D")
library("viridis")
library("ggplot2")
library("scales")
library("reshape2")
library("Matrix")
library("RColorBrewer")
library("gridExtra")

library(tikzDevice)
library(lhs)
library(RColorBrewer)

rm(list = ls())

#setwd(dirname(rstudioapi::getSourceEditorContext()$path))

fun <- function(x) {
  return(2*x[, 1]*x[, 3] + sin(x[, 2]*x[, 4]) + atan(3*x[,5]+5*x[, 6]))
}

D <- 6
multiplier <- 7
size <- D*multiplier
tolPrec <- 1e-4
tolCriteria <- 5e-7
alpha <- 1/2

constraint <- rep("monotonicity", D)

nbTrial <- 20
nbExp <- 1
maxiter <- 22

Hist <- list()
models <- list()
for (i in 1:nbTrial){
  xdesign <- lhsDesign(size,D, seed = 2*i)$design
  xdesign <- maximinSA_LHS(xdesign)$design
  ydesign <- fun(xdesign)
  model <- create(class = 'lineqBAGP', x = xdesign, y = ydesign, constrType = constraint)
  res <- BAGPMaxMod(model, max_iter = maxiter, print_iter = FALSE, nClusters = 12,
          GlobalconstrType= constraint, Estim_varnoise = TRUE, tolPrecision = tolPrec,
          tolCriteria = tolCriteria, Block_max_size = 3, save_models = TRUE, alpha = alpha)
  Hist[[i]] <- res$history
  models[[i]] <- res$hist_models 
}

#Computation of Q2 values

N <- 10^5
maxiter <- 17
Q2 <- matrix(0, nrow = nbTrial, ncol= maxiter) 
xtest <- lhsDesign(N, D, seed = i)$design
ytest <- fun(xtest)

for (i in 1:nbTrial){
  L <- length(models[[i]])
  for (j in 1:maxiter){
    if (j<=L){
    Q2[i,j] <- 1-norm(predict(models[[i]][[j]], xtest)$y.mod-ytest, type = "2")^2/(var(ytest)*N)
    } else {
      Q2[i,j] <- Q2[i,j-1]
    }
  }
}


Res <- list(Hist, models, Q2)
save(Res, file = "ResMaxMod")
unlink("/home/mderonzi/Documents/Mathis/thèse/lineqGPR/dev/lineqGPR/Res/ResMaxMod")

maxiter<- 17
Q2 <- Q2[,(1:maxiter)]

iterSeq <- seq(maxiter)
iterSeq <- factor(iterSeq, levels = iterSeq)
Q2median <- apply(Q2, 2, median, na.rm = TRUE)
idxQ2medSort <- sort(Q2median, index.return = TRUE)$ix
Q2median <- factor(Q2median, levels = Q2median)
ggData <- data.frame(iter = rep(iterSeq, each = nbTrial),
                     Q2 = c(Q2),
                     n = factor(D+1, levels = D+1),
                     Q2med = rep(Q2median, each = nbTrial))


viridisPalette <- viridis_pal(option= "A")(maxiter)
viridisPalette[idxQ2medSort] <- viridisPalette

#magmaPalette <- magma_pal()(maxiter)
viridisPalette[idxQ2medSort] <- viridisPalette

filename <- paste("Example1Q2boxplot.tex", sep = "")
tikz(filename, standAlone = TRUE, width = 6, height = 4)
ggplot(ggData, aes(x=iter, y=-log10(1-Q2), fill = Q2med)) +
  geom_boxplot(width =0.5, outlier.alpha = 0) +
  theme_bw() +
  theme(legend.position="bottom") +
  geom_jitter(width=0.1, alpha=0.5, size = 0.2) +
  scale_fill_manual(values=viridisPalette) + 
  scale_x_discrete(name ="MaxMod iteration") + 
  theme(legend.position = "none") +
  scale_y_continuous(name ="$-\\log(1-Q^2)$", breaks = 
  -round(log10(1-c(seq(0.2, 0.8, 0.2), 0.9, 0.95, 0.975, 0.99, 0.995, 0.999)), digits = 2),
  sec.axis = sec_axis(~ -10^(-.)+1 , name = "$Q^2$", 
  breaks = c(seq(0.2, 0.8, 0.2), 0.9, 0.95, 0.975, 0.99, 0.995, 0.999)))
dev.off()
tools::texi2dvi(filename,pdf=T,clean=TRUE)



# Installing the package


optDecision1 <- matrix(0, nbTrial, maxiter)
optDecision2 <- matrix(0, nbTrial, maxiter)
optDecision3 <- matrix(0, nbTrial, maxiter)

for(i in 1:nbTrial){
  for(j in 1:min(length(Hist[[i]]), maxiter)){
    if (length(Hist[[i]][[j]])==1){
      optDecision1[i,j] <- Hist[[i]][[j]] 
    } else{
      optDecision2[i,j] <- Hist[[i]][[j]][1]
      optDecision3[i,j] <- Hist[[i]][[j]][2]
    }
  }
}


#Plotting the Tab

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
  guides(color= guide_legend(), size=guide_legend()) +
  theme_bw()
filename <- "MaxModChoices.tex"
tikz(filename, standAlone = TRUE, width = 5, height = 5)
plot(fig1)
dev.off()
tools::texi2dvi(filename,pdf=T,clean=TRUE)
plot(fig1)

# fig2 <- ggplot(data=df, aes(activedim, iter)) +
#   geom_count(aes(activeblock1, iter, size=after_stat(..n..), color=after_stat(..n..)), shape = 8) +
#   geom_count(aes(activeblock2, iter, size=after_stat(..n..), color=after_stat(..n..)), shape = 8) +
#   scale_x_continuous(breaks = seq(D), position = "top", limits = c(1,D)) +
#   # theme(axis.text.x=element_text(angle = 45, hjust = 1)) +
#   #labs(title="", y = "MaxMod Iteration", x = "Input Dimension", colour = "reps", size = "reps", 
#   #     side= "left") +
#   labs(title="", y = "MaxMod Iteration", x = "Input Dimension", colour = "reps", size = "reps", 
#        side= "left") +
#   scale_size(breaks = seq(maxiter)) +
#   scale_color_continuous(breaks = seq(maxiter), low = redPalette[1], high = rev(redPalette)[1]) +
#   guides(color= guide_legend(), size=guide_legend()) +
#   theme_bw()
# plot(fig2)
