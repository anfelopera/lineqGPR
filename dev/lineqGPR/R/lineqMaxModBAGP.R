#' @title MaxMod algorithm for \code{"lineqGP"} and \code{"lineqBAGP"} Models
#' @description A wrapper function for the maximum Modification algorithm.
#' 
#' @param model an object with class \code{lineqGP} or \code{lineqBAGP}
# #' @param xtest test data for assessing the modification of the MAP estimate
#' @param tol a number corresponding to the tolerance of algorithm.
#' The algorithm stops if the (MaxMod criterion)/norm(pred_old) < tol
#' @param max_iter an integer corresponding to number of iterations
#' @param reward_new_knot a number corresponding to the reward of adding a new knot in an existing dimension
#' @param reward_new_dim a number corresponding to the reward of adding a new dimension
#' @param print_iter a logical variable to print results at each iteration
#' @param nClusters an integer corresponding to the number of clusters
#' @param save_history a logical variable to save the model at each iteration
#' @param xtest a vector of element we want to try our predictions
#' @param constrType a string that indicate the constraint we want that our predictor satisfies 
#' @param Block_max_size a number indicating the maximum block size
#'
# #' @param MCtest a logical variable to approximate MaxMod criterion via MC
#' 
#' @return an object with class \code{lineqGP} or \code{lineqBAGP} containing the resulting model
#'
#' @author M. Deronzier and A. F. Lopez-Lopera 
#' 
#' @export
BAGPMaxMod <- function(model, xtest=0,  max_iter = 5*ncol(model$x),
                       reward_new_knot, reward_new_dim = 1e-9,
                       print_iter = FALSE, nClusters = 20, tol = 5e-8,
                       save_history = FALSE, constrType, #contrType not handled for the moment,
                       Block_max_size = 4){
  D <- ncol(model$x)
  #Initialisation of the variables 
  partition <-list()
  nblock <- length(partition)
  subdivision <- list()
  activeVar <- rep(FALSE, D)
  pred <- NULL
  iter<- 1
  history <- vector("list", max_iter-1)
  #History of the choices and result of the square norm criterion
  #MaxMod_iter_values <- vector("")
  #rownames(MaxMod_iter_values) <- c("MaxMod criterion", "knot's position", "decision")
  #colnames(MaxMod_iter_values) <- paste("dim", 1:D)
  
  if (nClusters > 1) {
    if (!requireNamespace("foreach", quietly = TRUE))
      stop("Package \"foreach\" not found")
    if (!requireNamespace("doParallel", quietly = TRUE))
      stop("Package \"doParallel\" not found")
    
    nClusters_max <- max(parallel::detectCores(), 1)
    if (nClusters > nClusters_max) {
      warning("The maximum number of clusters has been exceeded.
      Further computations will consider only ", nClusters_max, " clusters")
      nClusters <- nClusters_max
    }
    cl <- parallel::makeCluster(nClusters)
    doParallel::registerDoParallel(cl)
  }
  
  #Start of the while loop 
  while(iter < max_iter) { #Number of iterations we'll do
    #Construction of the different choices
    #print("iteration:",iter, "partition:",partition)
    message("####################### ITERATION: ", iter, " ###########################")
    inactiveVar <- which(activeVar==FALSE)
    all_active <- eval(parse(text = paste("activeVar[", 1:D,"]", sep = "", collapse = "&&")))
    if(nblock>1){
      grid <- expand.grid(c(1:nblock),c(1:nblock))
      options <- grid[which((grid[,2]-grid[,1])>0),]
      options_expand <- append(as.list(c(1:D)), lapply(1:dim(options)[1], function(x) as.matrix(options[x,])))
    } else{
      options_expand <- as.list(c(1:D))
    }
    Change <- FALSE #Change is indicating if we modify the partition
    if (nClusters==1) {
      maximum <- -1
      for (i in 1:length(options_expand)) {
        MaxModcriterion_temp <- MaxModCriterionBAGP(model,iter,
                                                    options_expand,i,
                                                    constrType,
                                                    reward_new_dim,
                                                    reward_new_knot,
                                                    activeVar,
                                                    Block_max_size)
                    
        if ((MaxModcriterion_temp[[2]][1] > maximum) && (MaxModcriterion_temp[[2]][2]>=0)) {
            imax <- i
            maximum <- MaxModcriterion_temp[[2]][1]
            normmax <- MaxModcriterion_temp[[2]][2]
            #print("Max for the subdivision Subdivision")
            #print(MaxModcriterion_temp[[1]]$subdivision)
            Change <- TRUE
          }
        }
    } else if (nClusters >1) {
      Criteria <-try(foreach::"%dopar%"(foreach::foreach(i = (1:length(options_expand)),
                     .combine = cbind,.errorhandling ='remove'),{
          MaxModCriterionBAGP(model,iter,
                              options_expand, i,
                              constrType, reward_new_dim,
                              reward_new_knot, activeVar,
                              Block_max_size)[[2]]
        }))
      
      if (sum(Criteria[2,]>=0)>0){
        admitted <- which(Criteria[2,]>=0)
        maximum <- max(Criteria[1,admitted])
        imax <- which(Criteria[1,] == maximum)
        normmax <- Criteria[2,imax]
        Change <- TRUE
      }
      print("Criteria : ")
      print(Criteria)
    }
    # Update of variables
    if (!Change){
      message("Couldn't improve the prediction over observations")
      if (nClusters > 1)
        parallel::stopCluster(cl)
      return (list(model,history))
    }
    model <- MaxModCriterionBAGP(model,iter,options_expand,imax,constrType,
                 reward_new_dim,reward_new_knot,activeVar,Block_max_size)[[1]] #Pred can be removed
    #pred <- predict(model,0)
    print(model$kernParam$l)
    partition <- model$partition
    subdivision <- model$subdivision
    nblock <- model$localParam$nblocks
    activeVar[options_expand[[imax]][1]] <- TRUE
    if (maximum<tol){
      message("precision reached")
      if (nClusters > 1)
        parallel::stopCluster(cl)
      return(list(model, history))
    }
    history[[iter]] <- partition
    if(imax<=D){
      print(paste("activation variable ", imax, sep=))
    }
    else{
      print(paste("merging block", options_expand[[imax]][1], " and block",  options_expand[[imax]][2], sep = ""))
    }
    print(paste("imax :", imax, "  Criteria = ", maximum , "  diff_norm = ",
                normmax, sep = ""))
    #print(subdivision)
    print(partition)
    iter <- iter+1
  }
  if (nClusters > 1)
    parallel::stopCluster(cl)
  return(list(model, history)) 
}



#' @title MaxMod criterion for \code{"lineqBAGP"} Models
#' @description MaxMod criterion used to add a new knot or a new active dimension (see Bachoc et al., 2020)
#' 
#' @param model an object with class \code{lineqBAGP}
#' @param reward_new_dim a number corresponding to the reward of adding a new dimension
#' @param reward_new_knot a number corresponding to the reward of adding a new knot in an existing dimension
#' @param iter an integer corresponding to the iteration of the MaxMod algorithm 
#' @param constrType string that indicates the type of constraint that our predictor should satisfy
#' @param options_expand is a list of elements of type numeric that indicate a variable or a couple (variable,block)
#' @param i an integer, options_expand[i] indicating the choice we are making to increment the basis
#' @param activeVar a sequence of boolean, activeVar[k] indicates if varialbe k is activated or not 
#' @param Block_max_size a number indicating the maximum block size 
# #' @param MCtest a logical variable to approximate MaxMod criterion via MC
#' 
#' @return the value of the knot that minimize the MaxMod criterion and the 
#' value of the objective function 
#'
#' @author M. Deronzier
#'
#' @references F. Bachoc, A. F. Lopez-Lopera, and O. Roustant (2022),
#' "Sequential construction and dimension reduction of block additive Gaussian processes under inequality constraints".
#'
#' @importFrom utils tail
#' @export

MaxModCriterionBAGP <- function(model, iter, options_expand, i, constrType="none",
                                reward_new_dim, reward_new_knot,activeVar, Block_max_size) {
  D <- ncol(model$x)
  if (i<=D){ #Incrementing a variable 
    option_name <- "case1"
  } else{ #Merging two blocks
    option_name <- "case2"
  }
  option <- options_expand[[i]]
  if (iter == 1) { #iter==1 mean that the partition is empty
    model_update <- create(class = "lineqBAGP", x = model$x, y = model$y,
                           constrType ="none", 
                           partition = list(option),
                           subdivision = list(list(c(0,1))))
    #model_update <- optimise.parameters(model_update)
    Xi <- predict(model_update,0,0)$xi.mod
    return(list(model_update,c(Xi[1]*Xi[2] + ((Xi[1]-Xi[2])**2)/3,0)))
  } else {
    if (option_name=="case1"){# Creating a new block or a adding a knot in an already active variable
      if (activeVar[option[1]]){# Checking if the variable is already active
        new_t <- optimize(f = construct_t, interval = c(0, 1), model,
                          option[1], reward_new_knot = 1e-6#, Nscale = 1
        )[[1]]
        new.subdivision <- model$subdivision
        pos <- bijection(model$partition, option[1])
        new.subdivision[[pos[1]]][[pos[2]]] <- sort(c(model$subdivision[[pos[1]]][[pos[2]]], new_t))
        model_update <- create(class = "lineqBAGP", x = model$x, y = model$y,
                               constrType = model$constrType,
                               partition = model$partition, #To verify depends on model
                               subdivision = new.subdivision)
        model_update <- optimise.parameters(model_update)
      } else {
        new.partition <- model$partition
        new.subdivision <- model$subdivision
        new.partition[[model$nblock + 1]] <- option
        new.subdivision[[model$nblock + 1]] <- vector("list",1)
        new.subdivision[[model$nblock + 1]][[1]] <- c(0,1)
        model_update <- create(class = "lineqBAGP", x = model$x, y = model$y,
                               constrType = rep(model$constrType[1], length(new.partition)), # A voir plus tard
                               partition = new.partition,
                               subdivision = new.subdivision)
        model_update <- optimise.parameters(model_update)
      }
    } else if (option_name=="case2"){#We're merging two blocks
      size_basis <- lapply(lapply(model$subdivision, function(j) sapply(j, function(k) length(k))),
                           function(j) prod(j))
      block_size <- lapply(model$partition, function(x) length(x))
      if ( option[1]==option[2] || (block_size[[option[1]]]+block_size[[option[2]]])>Block_max_size){ 
      # if ((size_basis[option[[1]]]==2 || size_basis[[option[2]]]==2) || option[1]==option[2] ||
      #     (block_size[[option[1]]]+block_size[[option[2]]])>Block_max_size){ 
      #   #Checking if the blocks are allowed to ber merged
        return(list(model, as.matrix(c(0, -1))))
      } else{
        new.param <- merge_block(model$partition,model$subdivision, option[1],option[2])
        new.partition <- new.param[[1]]
        new.subdivision <- new.param[[2]]
        model_update <- create(class = "lineqBAGP", x = model$x, y = model$y,
                               constrType = rep(model$constrType[1], length(new.partition)),
                               partition = new.partition,
                               subdivision = new.subdivision
        )
        model_update <- optimise.parameters(model_update)
      }
    }
    criteria <- square_norm_int(model, model_update)
    diff_norm <-  norm(predict(model_update, model$x)$y.mod-model$y)/norm(mean(model$y)-model$y)
  }
  return(list(model_update, as.matrix(c(criteria, diff_norm))))
}



#' @title wrapper function 
#' @description Transfome model so it fits with the elements exepected to optimize the parameters
#' 
#' @param model an object with clall \code{lineqBAGP}

# #' @param NscaLe the number of knots we want to add
#' @return The wew model fitted to be optimised
#'
#' @author M. Deronzier 
#'
#' @references F. Bachoc, A. F. Lopez-Lopera, and O. Roustant (2020),
#' "Sequential construction and dimension reduction of Gaussian processes under inequality constraints".
#' \emph{ArXiv e-prints}
#' <arXiv:2009.04188>
#'
#' @export

wrapper <- function(model) {
  # pred <- predict(model, xtest = xtest)\
  new.partition <- lapply(model$partition, function(b) (1:length(b)))
  if(length(new.partition)>1){
    for (j in (2:length(new.partition))){
      new.partition[[j]] <- new.partition[[j]]+new.partition[[j-1]][length(new.partition[[j-1]])]
    }
  }
  new.x <- model$x[,unlist(model$partition)] 
  new.model <- create(class = "lineqBAGP", x = new.x, y = model$y,
                      constrType = rep(model$constrType[1], length(new.partition)), 
                      partition = new.partition,
                      subdivision = model$subdivision
  )
  return(new.model)
}

#' @title Optimising parameters 
#' @description Function that optimise parameters
#' 
#' @param model an object with clall \code{lineqBAGP}

# #' @param NscaLe the number of knots we want to add
#' @return The wew model fitted to be optimised
#'
#' @author M. Deronzier 
#'
#' @references F. Bachoc, A. F. Lopez-Lopera, and O. Roustant (2020),
#' "Sequential construction and dimension reduction of Gaussian processes under inequality constraints".
#' \emph{ArXiv e-prints}
#' <arXiv:2009.04188>
#'
#' @export

optimise.parameters <- function(model) {
  # pred <- predict(model, xtest = xtest)\
  new.model <- wrapper(model)
  new.model <- lineqGPOptim(new.model,
                            additive = TRUE, # if the model is additive
                            block = TRUE, # if the model is additive per blocks
                            estim.varnoise = TRUE, # to add this info at the MaxMod level
                            bounds.varnoise = c(1e-3, Inf), # to add this info at the MaxMod level
                            lb = rep(1e-2, new.model$d+new.model$localParam$nblocks),
                            #ub = c(Inf, 0.7, 0.7, Inf, 0.7), # to add this info at the MaxMod level
                            ub = rep(Inf, new.model$d+new.model$localParam$nblocks),
                            opts = list(#algorithm = "NLOPT_LD_MMA",
                                        algorithm = "NLOPT_LN_COBYLA",
                                        print_level = 0,
                                        ftol_abs = 1e-3,
                                        maxeval = 30,
                                        check_derivatives = TRUE)
                          )
  new.model$partition <- model$partition 
  lapply(1:length(new.model$partition), function(j) names(new.model$partition[[j]]) <- names(model$partition[[j]])) 
  lapply(1:length(new.model$subdivision), function(j) names(new.model$partition[[j]]) <- names(model$subdivision[[j]])) 
  new.model$x <- model$x
  new.model$d <- model$d
  return(new.model)
}


#' @title Choice on the posittion of the node 
#' @description Compute the modification of the MAP estimate according to the MaxMod criterion
#' 
#' @param model an object with clall \code{lineqBAGP}
#' @param choice the choice made for the variable where there will be a refinement of the subdivision   
#' @param t the position of the refinement
#' @param reward_new_knot a number corresponding to the reward of adding a new knot in an existing dimension
# #' @param NscaLe the number of knots we want to add
#' @return the modification of the MAP estimate after adding a new knot
#'
#' @author M. Deronzier
#'
#' @references F. Bachoc, A. F. Lopez-Lopera, and O. Roustant (2020),
#' "Sequential construction and dimension reduction of Gaussian processes under inequality constraints".
#' \emph{ArXiv e-prints}
#' <arXiv:2009.04188>
#'
#' @export

construct_t <- function(t, model, choice,
                        reward_new_knot = 1e-6 #Nscale = 1
) {
  subdivision2 <- model$subdivision
  pos <- bijection(model$partition, choice[1])
  subdivision2[[pos[1]]][[pos[2]]] <- sort(c(model$subdivision[[pos[1]]][[pos[2]]],t))
  model2 <- create(class = "lineqBAGP", x = model$x, y = model$y,
                   constrType = model$constrType, 
                   partition = model$partition,
                   subdivision = subdivision2)
  model2$varnoise <- model$varnoise
  model2$kernParam <- model$kernParam
  #optimise.parameters(model2)
  return(-square_norm_int(model,model2))
}



