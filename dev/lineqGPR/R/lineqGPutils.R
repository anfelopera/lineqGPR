#' @title Model Creations
#' @description Wrapper function for creations of model functions.
#' The function invokes particular methods which depend on the \code{class}
#' of the first argument.
#' 
#' @param class a character string corresponding to the desired class.
#' @param ... further arguments passed to or from other methods.
#' (see, e.g., \code{\link{create.lineqGP}})
#' 
#' @return A model object created according to its \code{class}.
#'
#' @seealso \code{\link{augment}}, \code{\link{predict}},
#'          \code{\link{simulate}}
#'          
#' @author A. F. Lopez-Lopera
#'
#' @examples
#' ## Not run:
#' model <- list()
#' model2 <- create(class = "ClassName", model)
#' model2
#' ## End(Not run)
#'
#' @export
create <- function(class, ...) {
  fun <- paste("create.", class, sep = "")
  fun <- try(get(fun))
  if (class(fun) == "try-error") {
    class(model) <- class
    warning('class "', class, '" is not supported')
  } else {
    model <- fun(...)
    class(model) <- class
  }
  return(model)
}

#' @title Sampling Methods of Truncated Multivariate Normal Distributions
#' @description Wrapper function with a collection of Monte Carlo and Markov Chain
#' Monte Carlo samplers for truncated multivariate normal distributions.
#' The function invokes particular samplers which depend on the class of the first argument.
#' 
#' @param object an object with:
#'        \code{mu} (mean vector), \code{Sigma} (covariance matrix),
#'        \code{lb} (lower bound vector), \code{ub} (upper bound vector).
#' @param nsim an integer corresponding to the number of simulations.
#' @param ... further arguments passed to or from other methods.
#' 
#' @return A matrix with the sample path. Samples are indexed by columns.
#'
#' @seealso \code{\link{tmvrnorm.RSM}}, \code{\link{tmvrnorm.HMC}}, 
#'          \code{\link{tmvrnorm.ExpT}}
#'          
#' @author A. F. Lopez-Lopera
#'
#' @export
tmvrnorm <- function(object, nsim, ...)
  UseMethod("tmvrnorm")

#' @title Training/test data generator according to a given Design of Experiment (DoE)
#' @description Split the data in training/test sets according to a given DoE.
#' 
#' @param x a vector (or matrix) with the input locations.
#' @param y a vector with the output observations.
#' @param DoE.idx the numeric indices of the training data used in the design.
#' @param DoE.type if \code{is.null(DoE.idx)}, a character string
#' corresponding to the type of DoE.
#' Options: \code{rand} (random desings), \code{regs} (regular-spaced desings).
#' @param ratio  if \code{is.null(DoE.idx)}, a number with the ratio
#' nb_train/nb_total (by default, ratio = 0.3).
#' @param seed an optional value corresponding to the seed for random methods.
#' 
#' @return A list with the DoE: \code{list(xdesign, ydesign, xtest, ytest)}.
#' 
#' @section Comments:
#' This function is in progress. Other types of DoEs will be considered using
#' the \code{DiceDesign} package.
#'
#' @author A. F. Lopez-Lopera
#'
#' @examples
#' # generating the toy example
#' x <- seq(0, 1, length = 100)
#' y <- sin(4*pi*x)
#'
#' # regular DoE
#' DoE <- splitDoE(x, y, DoE.type = "regs", ratio = 0.3)
#' plot(x,y)
#' points(DoE$xdesign, DoE$ydesign, col = "red", pch = 20)
#' points(DoE$xtest, DoE$ytest, col = "blue", pch = 20)
#' legend("topright", c("training data", "test data"),
#'        pch = rep(20, 2), col = c("red", "blue"))
#'
#' # random DoE
#' DoE <- splitDoE(x, y, DoE.type = "rand", ratio = 0.3, seed = 1)
#' plot(x,y)
#' points(DoE$xdesign, DoE$ydesign, col = "red", pch = 20)
#' points(DoE$xtest, DoE$ytest, col = "blue", pch = 20)
#' legend("topright", c("training data", "test data"),
#'        pch = rep(20, 2), col = c("red", "blue"))
#'
#' @export
splitDoE <- function(x, y, DoE.idx = NULL,
                     DoE.type = c("rand", "regs"),
                     ratio = 0.3, seed = NULL) {
  if (!is.null(DoE.idx) && max(DoE.idx) > length(y))
    stop('max(idx.DoE) cannot be larger than length(y)')
  x <- as.matrix(x)
  if (is.null(DoE.idx)) {
    nb_obs <- length(y)
    nb_DoE <- round(ratio*nb_obs)
    set.seed(seed)
    DoE.type <- match.arg(DoE.type)
    switch(DoE.type,
           rand = {
             temp <- sample(nb_obs)
             DoE.idx <- temp[1:nb_DoE]},
           regs = {
             DoE.idx <- round(seq(1, nb_obs, length.out = nb_DoE))}
    )
  }
  return(list(xdesign = x[DoE.idx, ], ydesign = y[DoE.idx],
              xtest = x[-DoE.idx, ], ytest = y[-DoE.idx]))
}

#' @title Error Measures for GP Models.
#' @description Compute error measures for GP models:
#' mean absulte error (\code{"mae"}), mean squared error (\code{"mse"}),
#' standardised mse (\code{"smse"}), mean standardised log loss (\code{"msll"}),
#' Q2 (\code{"q2"}), predictive variance adequation (\code{"pva"}),
#' confidence interval accuracy (\code{"cia"}).
#' 
#' @param y a vector with the output observations used for training.
#' @param ytest a vector with the output observations used for testing.
#' @param mu a vector with the posterior mean.
#' @param varsigma a vector with the posterior variances.
#' @param type a character string corresponding to the type of the measure.
#' @param control an optional list with parameters to be passed (e.g. cia: "nsigma").
#' 
#' @return The values of the error measures.
#'
#' @seealso \code{\link{errorMeasureRegressMC}}
#' 
#' @author A. F. Lopez-Lopera
#'
#' @references C. E. Rasmussen, and C. K. I. Williams (2005),
#' "Gaussian Processes for Machine Learning (Adaptive Computation and Machine Learning)".
#' \emph{The MIT Press}.
#'
#' F. Bachoc (2013),
#' "Cross validation and maximum likelihood estimations of hyper-parameters of
#' Gaussian processes with model misspecification".
#' \emph{Computational Statistics & Data Analysis},
#' 66:55-69.
#'
#' @examples
#' # generating the toy example
#' n <- 100
#' w <- 4*pi
#' x <- seq(0, 1, length = n)
#' y <- sin(w*x)
#'
#' # results with high-level noises generating the toy example
#' nbsamples <- 100
#' set.seed(1)
#' ynoise <- y + matrix(rnorm(n*nbsamples, 0, 10), ncol = nbsamples)
#' mu <- apply(ynoise, 1, mean)
#' sigma <- apply(ynoise, 1, sd)
#' matplot(x, ynoise, type = "l", col = "gray70")
#' lines(x, y, lty = 2, col = "red")
#' lines(x, mu, col = "blue")
#' lines(x, mu+1.98*sigma, lty = 2)
#' lines(x, mu-1.98*sigma, lty = 2)
#' legend("topright", c("target", "mean", "confidence", "samples"),
#'        lty = c(2,1,2,1), col = c("red", "blue", "black", "gray70"))
#' t(errorMeasureRegress(y, y, mu, sigma^2))
#'
#' # results with low-level noises generating the toy example
#' set.seed(1)
#' ynoise <- y + matrix(rnorm(n*nbsamples, 0, 0.05), ncol = nbsamples)
#' mu <- apply(ynoise, 1, mean)
#' sigma <- apply(ynoise, 1, sd)
#' matplot(x, ynoise, type = "l", col = "gray70")
#' lines(x, y, lty = 2, col = "red")
#' lines(x, mu, col = "blue")
#' lines(x, mu+1.98*sigma, lty = 2)
#' lines(x, mu-1.98*sigma, lty = 2)
#' legend("topright", c("target", "mean", "confidence", "samples"),
#'        lty = c(2,1,2,1), col = c("red", "blue", "black", "gray70"))
#' t(errorMeasureRegress(y, y, mu, sigma^2))
#'
#' @export
errorMeasureRegress <- function(y, ytest, mu, varsigma,  type = "all",
                                control = list(nsigma = 1.96)) {
  # computing the errors
  maerror <- abs(ytest - mu)
  mserror <- (ytest - mu)^2
  std_mserror <- mserror/var(ytest)
  q2 <- 1 - sum(mserror)/sum((ytest - mean(ytest))^2)
  pva <- mserror/varsigma
  sserror <- ((ytest - mu)^2)/(2*varsigma)
  logprob <- 0.5*log(2*pi*varsigma) + sserror -
             0.5*log(2*pi*var(y)) - ((ytest - mean(y))^2)/(2*var(y))
  logprob <- logprob[complete.cases(logprob)]
  cia <- ytest >= (mu-control$nsigma*sqrt(varsigma)) &
         ytest <= (mu+control$nsigma*sqrt(varsigma))

  error <- c(mae = mean(maerror), mse = mean(mserror),
             smse = mean(std_mserror), msll = mean(logprob),
             q2 = q2, pva = abs(log(mean(pva))), cia = mean(cia))
  switch(type,
    mae = {return(error["mae"])},
    mse = {return(error["mse"])},
    smse = {return(error["smse"])},
    msll = {return(error["msll"])},
    q2 = {return(error["q2"])},
    pva = {return(error["pva"])},
    cia = {return(error["cia"])},
    {return(as.matrix(error, nrow = 1))}
  )
}

#' @title Error Measures for GP Models using Monte Carlo Samples.
#' @description Compute error measures for GP models using Monte Carlo samples:
#' mean absulte error (\code{"mae"}), mean squared error (\code{"mse"}),
#' standardised mse (\code{"smse"}), Q2 (\code{"q2"}), predictive variance
#' adequation (\code{"pva"}), confidence interval accuracy (\code{"cia"}).
#' 
#' @param y a vector with the output observations used for training.
#' @param ytest a vector with the output observations used for testing.
#' @param ysamples a matrix with posterior sample paths. Samples are indexed by columns.
#' @param type a character string corresponding to the type of the measure.
#' @param control an optional list with parameters to be passed (cia: "probs").
#' 
#' @return The values of the error measures.
#'
#' @seealso \code{\link{errorMeasureRegress}}
#' 
#' @author A. F. Lopez-Lopera
#'
#' @examples
#' # generating the toy example
#' n <- 100
#' w <- 4*pi
#' x <- seq(0, 1, length = n)
#' y <- sin(w*x)
#'
#' # results with high-level noises generating the toy example
#' nbsamples <- 100
#' set.seed(1)
#' ynoise <- y + matrix(rnorm(n*nbsamples, 0, 10), ncol = nbsamples)
#' matplot(x, ynoise, type = "l", col = "gray70")
#' lines(x, y, lty = 2, col = "red")
#' legend("topright", c("target", "samples"), lty = c(2,1), col = c("red", "gray70"))
#' t(errorMeasureRegressMC(y, y, ynoise))
#'
#' # results with low-level noises generating the toy example
#' set.seed(1)
#' ynoise <- y + matrix(rnorm(n*nbsamples, 0, 0.05), ncol = nbsamples)
#' matplot(x, ynoise, type = "l", col = "gray70")
#' lines(x, y, lty = 2, col = "red")
#' legend("topright", c("target", "samples"), lty = c(2,1), col = c("red", "gray70"))
#' t(errorMeasureRegressMC(y, y, ynoise))
#'
#' @importFrom stats quantile
#' @export
errorMeasureRegressMC <- function(y, ytest, ysamples, type = "all",
                                  control = list(probs = c(0.05,0.95))) {
  # precomputing some values
  mu <- apply(ysamples, 1, mean)
  varsigma <- apply(ysamples, 1, sd)^2
  qtls <- apply(ysamples, 1, quantile, probs = control$probs)

  # computing the errors
  maerror <- abs(ytest - mu)
  mserror <- (ytest - mu)^2
  std_mserror <- mserror/var(ytest)
  q2 <- sum(mserror)/sum((ytest - mean(ytest))^2)
  pva <- mserror/varsigma
  cia <- ytest >=  qtls[1,] & ytest <= qtls[2,]

  error <- c(mae = mean(maerror), mse = mean(mserror),
             smse = mean(std_mserror),
             q2 = 1 - q2, pva = abs(log(mean(pva))), cia = mean(cia))
  switch(type,
         mae = {return(error["mae"])},
         mse = {return(error["mse"])},
         smse = {return(error["smse"])},
         q2 = {return(error["q2"])},
         pva = {return(error["pva"])},
         cia = {return(error["cia"])},
         {return(as.matrix(error, nrow = 1))}
  )
}

#' @title Linear Systems of Inequalities
#' @description Build the linear system of inequalities given specific bounds.
#' 
#' @param d the number of linear inequality constraints.
#' @param l the value (or vector) with the lower bound.
#' @param u the value (or vector) with the upper bound.
#' @param A a matrix containing the structure of the linear equations.
#' @param lineqSysType a character string corresponding to the type of the
#' linear system. Options: \code{twosides}, \code{oneside}. \cr
#'     - \code{twosides} : Linear system given by
#'       \deqn{\boldsymbol{l} \leq \boldsymbol{A x} \leq \boldsymbol{u}.}{l \le A x \le u.}
#'
#'     - \code{oneside} : Extended linear system given by
#'       \deqn{\boldsymbol{M x} + \boldsymbol{g} \geq \boldsymbol{0} \quad \mbox{with} \quad \boldsymbol{M} = [\boldsymbol{A}, -\boldsymbol{A}]^\top
#'        \quad \mbox{and} \quad \boldsymbol{g} = [-\boldsymbol{l}, \boldsymbol{u}]^\top.}{M x + g \ge 0 with M = [A, -A]^T and g = [-l, u]^T.}
#' @param rmInf If \code{TRUE}, inactive constraints are removed
#' (e.g. \eqn{-\infty \leq x \leq \infty}{-Inf \le x \le Inf}).
#' 
#' @return  A list with the linear system of inequalities:
#' \code{list(A,l,u)} (\code{twosides}) or \code{list(M,g)} (\code{oneside}).
#'
#' @author A. F. Lopez-Lopera
#'
#' @examples
#' n <- 5
#' A <- diag(n)
#' l <- rep(0, n)
#' u <- c(Inf, rep(1, n-1))
#' bounds2lineqSys(n, l, u, A, lineqSysType = "twosides")
#' bounds2lineqSys(n, l, u, A, lineqSysType = "oneside", rmInf = FALSE)
#' bounds2lineqSys(n, l, u, A, lineqSysType = "oneside", rmInf = TRUE)
#'
#' @export
bounds2lineqSys <- function(d = nrow(A), l = 0, u = 1, A = diag(d),
                            lineqSysType = "twosides", rmInf = TRUE) {
  if (!is.element(length(l), c(1, d)) || !is.element(length(u), c(1, d)))
    stop('The bounds have to be d-dimensionals')
  if (nrow(A) != d)
    stop('nrow(A) has to be equal to d')
  if (any(l >= u))
    stop('The elements from "u" has to be greater than the elements from "l"')

  if (length(l) == 1) l <- rep(l, d)
  if (length(u) == 1) u <- rep(u, d)

  switch(lineqSysType,
         twosides = {
           if (rmInf && any(l == -Inf & u == Inf)) {
             idx <- which(l == -Inf & u == Inf)
             l <- l[-idx]
             u <- u[-idx]
             A <- A[-idx, ]
           }
           return(list(A = A, l = l, u = u))
         }, oneside = {
           g <- c(-l, u)
           M <- rbind(A, -A)
           if (rmInf && any(g == Inf)) {
             idx <- which(g == Inf)
             g <- g[-idx]
             M <- M[-idx, ]
             # if (length(g) == 0) stop('bounds not defined')
           }
           if (!is.matrix(M)) ##
             M <- matrix(M, nrow = 1) ##
           return(list(M = M, g = g))
         }, {
           stop('type of linear system "', lineqSysType, '" is not supported')
         }
  )
}

#' @title Plot for the \code{"lineqGP"} S3 Class
#' @description Plot for the \code{"lineqGP"} S3 class.
#' 
#' @param x an object with \code{"lineqGP"} S3 class.
#' @param y not used.
#' @param ytest the values of the test observations. If \code{!is.null(ytest)}, ytest is drawn.
#' @param probs the values of the confidence intervals evaluated at probs.
#' @param bounds the values of the bounds of a constrained model. If \code{!is.null(bounds)}, bounds are drawn.
#' @param addlines pptional Logical. If \code{TRUE}, some samples are drawn.
#' @param nblines if \code{addlines}. The number of samples to be drawn.
#' @param ... further arguments passed to or from other methods.
#' 
#' @return Plot with the \code{"lineqGP"} model.
#'
#' @seealso \code{\link{plot}}, \code{\link{ggplot.lineqGP}}
#' 
#' @author A. F. Lopez-Lopera
#'
#' @import graphics stats grDevices
#' @export
plot.lineqGP <- function(x, y, ytest = NULL,
                         probs = c(0.025, 0.975), bounds = NULL,
                         addlines = TRUE, nblines = 5, ...) {
  model <- x
  xtrain <- model$x
  ytrain <- model$y
  xtest <- model$xtest
  ysim <- model$ysim

  # precomputing some terms
  qtls <- apply(ysim, 1, quantile, probs = probs)
  ymean <- apply(ysim, 1, mean)

  # drawing the model according to the method
  plot(c(xtest, rev(xtest)), c(qtls[2, ], rev(qtls[1, ])),
       col = "gray90", type = "l", lty = 1, ...)
  polygon(c(xtest, rev(xtest)), c(qtls[2, ], rev(qtls[1, ])),
          col = "gray90", border = NA)
  if (addlines){
    set.seed(nblines)
    matplot(xtest, ysim[, sample(ncol(ysim),nblines)],
            type ='l', lty = 3, add = TRUE,
            col = rainbow(nblines, 0.7), ...)
  }
  lines(xtest, ymean, lty = 1, col = "darkgreen")
  points(xtrain, ytrain, pch = 20, ...)
  if(!is.null(bounds) && length(bounds) == 2)
    abline(h = bounds, lty = 2)
  if(!is.null(ytest))
    points(xtest, ytest, pch = 4, col = "red", ...)
  fig <- recordPlot()
  return(fig)
}

#' @title GGPlot for the \code{"lineqGP"} S3 Class
#' @description GGPlot for the \code{"lineqGP"} S3 class.
#' 
#' @param data an object with \code{"lineqGP"} S3 class.
#' @param mapping not used.
#' @param ... further arguments passed to or from other methods.
#' @param ytest the values of the test observations. \code{If !is.null(ytest)}, ytest is drawn.
#' @param probs the values of the confidence intervals evaluated at probs.
#' @param bounds the values of the bounds of a constrained model. If \code{!is.null(bounds)}, bounds are drawn.
#' @param addlines an optional Logical. If \code{TRUE}, some samples are drawn.
#' @param nblines if \code{addlines}. The number of samples to be drawn.
#' @param fillbackground an optional logical. If \code{TRUE}, fill gray background.
#' @param alpha.qtls a number indicating the transparency of the quantiles.
#' @param xlab a character string corresponding to the title for the x axis.
#' @param ylab a character string corresponding to the title for the y axis.
#' @param main a character string corresponding to the overall title for the plot.
#' @param xlim the limit values for the x axis.
#' @param ylim the limit values for the y axis.
#' @param lwd a number indicating the line width.
#' @param cex a number indicating the amount by which plotting text and symbols should be scaled.
#' 
#' @return GGPlot with the \code{"lineqGP"} model.
#'
#' @seealso \code{\link{ggplot}}, \code{\link{plot.lineqGP}}
#' 
#' @author A. F. Lopez-Lopera
#'
#' @import ggplot2
#' @export
ggplot.lineqGP <- function(data, mapping, ytest = NULL,
                           probs = c(0.05, 0.95), bounds = NULL,
                           addlines = TRUE, nblines = 5,
                           fillbackground = TRUE, alpha.qtls = 0.4,
                           xlab = "", ylab = "", main = "",
                           xlim = NULL, ylim = NULL,
                           lwd = 1, cex = 1.5, ...) {
  # alpha.qtls <- toString(alpha.qtls) # a bug for updated ggplot2 package
  model <- data
  xtrain <- model$x
  ytrain <- model$y
  xtest <- model$xtest
  ysim <- model$ysim

  # precomputing some terms
  qtls <- apply(ysim, 1, quantile, probs = probs)
  ymean <- apply(ysim, 1, mean)

  # drawing the model according to the method
  if (is.null(xlim))
    xlim <- range(c(xtrain, xtest))
  if (is.null(ylim))
    ylim <- range(qtls)
  fig <- ggplot() + labs(title = main) + 
    scale_x_continuous(name = xlab) + scale_y_continuous(name = ylab) +
    coord_cartesian(xlim = xlim, ylim = ylim)
  if (!fillbackground)
    fig <- fig + theme_bw()
  fig <- fig + theme(axis.text = element_text(size = 20),
                     axis.title = element_text(size = 25),
                     legend.position = "none") +
    geom_ribbon(aes(x = xtest, ymin = qtls[1, ], ymax = qtls[2, ]),
                colour ="gray90", fill = "gray60", alpha = alpha.qtls)
  if (addlines) {
    set.seed(nblines)
    ggData <- data.frame(x = xtest,
                         y = matrix(ysim[, seq(nblines)], ncol = 1),
                         label = rep(seq(nblines), each = nrow(ysim)))
      fig <- fig + geom_line(mapping = aes(x = ggData$x, y = ggData$y,
                                           group = ggData$label,
                                           colour = as.factor(ggData$label)),
                           linetype = "longdash", size = 0.8*lwd)
  }
  fig <- fig + geom_line(aes(x = xtest, y = ymean),
                         colour = "#3366FF", size = lwd) +
    geom_point(aes(xtrain, ytrain),
               colour = "black", shape = 20, size = cex, stroke = lwd)
  if(!is.null(bounds) && length(bounds) == 2)
    fig <- fig + geom_hline(yintercept = as.numeric(bounds),
                            linetype = "dashed", size = lwd)
  if(!is.null(ytest)) {
    fig <- fig + geom_point(aes(x = xtest, y = ytest),
                            colour = "red", shape = 4, size = cex, stroke = lwd)
  }
  par(cex.axis = 1.5, cex.lab = 2.0, lwd = lwd)
  plot(fig)
  return(fig)
}

# #' @title GGPlot for the \code{"lineqDGP"} S3 Class
# #' @description GGPlot for the \code{"lineqDGP"} S3 class.
# #' See \code{\link{ggplot.lineqGP}} for more details.
# #' 
# #' @param data an object with \code{lineqDGP} S3 class.
# #' @param mapping not used.
# #' @param ... further arguments passed to or from other methods.
# #' 
# #' @return GGPlot with the \code{"lineqDGP"} model.
# #'
# #' @seealso \code{\link{ggplot.lineqGP}}, \code{\link{ggplot}}
# #' 
# #' @author A. F. Lopez-Lopera
# #'
# #' @method ggplot lineqDGP
# #' @export
# ggplot.lineqDGP <- function(data, mapping, ...)
#   ggplot.lineqGP(data, mapping, ...)

#' @title Plot for the \code{"lineqAGP"} S3 Class
#' @description Plot for the \code{"lineqAGP"} S3 class.
#' See \code{\link{plot.lineqGP}} for more details.
#' 
#' @param x an object with \code{"lineqAGP"} S3 class.
#' @param y not used.
#' @param ... further arguments passed to or from other methods.
#' 
#' @return Plot with the \code{"lineqAGP"} model.
#'
#' @seealso \code{\link{ggplot.lineqGP}}, \code{\link{plot}}
#' 
#' @author A. F. Lopez-Lopera
#'
#' @method plot lineqAGP
#' @export
plot.lineqAGP <- function(x, y, ...)
  plot.lineqGP(x, y, ...)

#' @title List of blocks
#' @description Create a list according to a given a partition of set $\{1,cdots,D\}$.
#' 
#' @param partition a list containing a partition of set $\{1,cdots,D\}$
#'
#' @return a list structure according to the partition.
#'          
#' @author A. F. Lopez-Lopera
#'
#' @examples
#' partition <- list(c(1,3), 2)
#' listBlocks(partition)
#' @export
listBlocks <- function(partition) {
  nblocks <- length(partition)
  dblock <- sapply(partition, length)
  object <- vector("list", nblocks)
  names(object) <- paste("block", 1:nblocks, sep = "")
  
  for (j in 1:nblocks) {
    object[[j]] <- vector("list", dblock[j])
    names(object[[j]]) <- paste("x", partition[[j]], sep = "")
  }
  return(object)
}




