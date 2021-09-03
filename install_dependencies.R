#### run to install all the dependencies of the package lineqGPR ####
# list of dependencies
listPackages <- c("nloptr", "broom", "tmg", "mvtnorm","MASS", "quadprog",
                  "Matrix", "restrictedMVN", "TruncatedNormal", "ggplot2",
                  "viridis", "purrr", "plot3D")
# installing dependencies
install.packages(listPackages)
