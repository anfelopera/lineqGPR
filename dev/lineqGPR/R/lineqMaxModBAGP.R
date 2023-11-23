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
#' @author A. F. Lopez-Lopera
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
  criteria<-100
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
  while((iter < max_iter)&& criteria>1e-7) { #Number of iterations we'll do
    #Construction of the different choices
    #print("iteration:",iter, "partition:",partition)
    print(paste("####################### ITERATION : ", iter, " ###########################", sep=""))
    inactiveVar <- which(activeVar==FALSE)
    all_active <- eval(parse(text = paste("activeVar[", 1:D,"]", sep = "", collapse = "&&")))
    if(nblock>1){
      grid <- expand.grid(c(1:nblock),c(1:nblock))
      options <- grid[which((grid[,2]-grid[,1])>0),]
      options_expand <- append(as.list(c(1:D)), lapply(1:dim(options)[1], function(x) as.matrix(options[x,])))
    } else{
      options_expand <- as.list(c(1:D))
    }
    if (nClusters==1) {
      for (i in 1:length(options_expand)) {
        print(paste("iter:", iter, " i:", i, " choice:", options_expand[[i]]))
        MaxModcriterion_temp <- MaxModCriterionBAGP(
          model,iter,pred,
          options_expand,i,
          constrType,
          reward_new_dim,
          reward_new_knot,
          activeVar,
          Block_max_size)
        print(paste("Criteria = ", MaxModcriterion_temp[[2]], sep = ""))
        if (i == 1) {
          #Initialisation of the minimum value
          imax <- 1
          maximum <- MaxModcriterion_temp[[2]]
          print("Max for the subdivision Subdivision")
          print(MaxModcriterion_temp[[1]]$subdivision)
        } else {
          #The value has been initialized check if the new model is better or not
          if (MaxModcriterion_temp[[2]] > maximum) {
            imax <- i
            maximum <- MaxModcriterion_temp[[2]]
            print("Max for the subdivision Subdivision")
            print(MaxModcriterion_temp[[1]]$subdivision)
          }
        }
      }
    } else if (nClusters >1) {
      Criteria <-
        try(foreach::"%dopar%"(foreach::foreach(i = (1:length(options_expand)),
          .combine = cbind,.errorhandling ='remove'),
        {
          MaxModCriterionBAGP(
            model,iter,pred,
            options_expand, i,
            constrType, reward_new_dim,
            reward_new_knot, activeVar,
            Block_max_size)[[2]]
        }))
      maximum <- max(Criteria)
      imax <- which(Criteria == maximum)
      
    }
    # Update of variables
    model <- MaxModCriterionBAGP(model,iter,pred,options_expand,imax,constrType,
                 reward_new_dim,reward_new_knot,activeVar,Block_max_size)[[1]]
    #pred <- predict(model,0)
    partition <- model$partition
    subdivision <- model$subdivision
    nblock <- model$localParam$nblocks
    activeVar[options_expand[[imax]][1]] <- TRUE
    if (maximum<tol){
      return(model)
    }
    iter <- iter+1
    print(model$subdivision)
  }
  
  #Here we should update the history of the choices
  
  return(model) 
}


#' #' @title MaxMod criterion for \code{"lineqGP"} Models
#' #' @description MaxMod criterion used to add a new knot or a new active dimension (see Bachoc et al., 2020)
#' #'
#' #' @param model an object with class \code{lineqGP}
#' #' @param ulist a list with the location of the knots
#' #' @param xtest test data for assessing the modification of the MAP estimate
#' #' @param activeDim a sequence containing the active dimensions
#' #' @param activeDim_IdxSeq a sequence containing the indices of the actived input dimensions
#' #' @param idx_add index of the dimension to be activated
#' #' @param reward_new_dim a number corresponding to the reward of adding a new dimension
#' #' @param reward_new_knot a number corresponding to the reward of adding a new knot in an existing dimension
#' #' @param iter an integer corresponding to the iteration of the MaxMod algorithm
#' #' @param pred an object with class \code{lineqGP} containing the predictive model
#' #'
#' #' @return the value of the knot that minimize the MaxMod criterion and the
#' #' value of the objective function
#' #'
#' #' @author A. F. Lopez-Lopera
#' #'
#' #' @references F. Bachoc, A. F. Lopez-Lopera, and O. Roustant (2020),
#' #' "Sequential construction and dimension reduction of Gaussian processes under inequality constraints".
#' #' \emph{ArXiv e-prints}
#' #' <arXiv:2009.04188>
#' #'
#' #' @export
#' #'
#' MaxModCriterion <- function(model, iter, pred , option_name , option,
#'                             reward_new_dim, reward_new_knot
#' ) {
#'   class <- class(model)
#'   fun <- paste("MaxModCriterion.", class, sep = "")
#'   fun <- try(get(fun))
#'   if (class(fun) == "try-error") {
#'     stop('class "', class, '" is not supported')
#'   } else {
#'     model <- fun(model, iter, pred , option_name , option,
#'                  reward_new_dim, reward_new_knot
#'     )
#'   }
#'   return(model)
#' }

#' @title MaxMod criterion for \code{"lineqBAGP"} Models
#' @description MaxMod criterion used to add a new knot or a new active dimension (see Bachoc et al., 2020)
#' 
#' @param model an object with class \code{lineqBAGP}
#' @param reward_new_dim a number corresponding to the reward of adding a new dimension
#' @param reward_new_knot a number corresponding to the reward of adding a new knot in an existing dimension
#' @param iter an integer corresponding to the iteration of the MaxMod algorithm 
#' @param pred an object with class \code{lineqBAGP} containing the predictive model
#' @param constrType string that indicates the type of constraint that our predictor should satisfy
#' @param option_name string that indicates in which case we are for the construction of the partition
#' @param option is an element of type numeric that indicate a variable or a couple (variable,block)
#' @param activeVar a sequence of boolean, activeVar[i] indicates if i is activated or not 
#' @param Block_max_size a number indicating the maximum block size 
# #' @param MCtest a logical variable to approximate MaxMod criterion via MC
#' 
#' @return the value of the knot that minimize the MaxMod criterion and the 
#' value of the objective function 
#'
#' @author A. F. Lopez-Lopera
#'
#' @references F. Bachoc, A. F. Lopez-Lopera, and O. Roustant (2022),
#' "Sequential construction and dimension reduction of block additive Gaussian processes under inequality constraints".
#'
#' @importFrom utils tail
#' @export

MaxModCriterionBAGP <- function(model, iter, pred , options_expand, i, constrType="none",
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
    Xi <- predict(model_update,0)$xi.mean
    return(list(model_update,Xi[1]*Xi[2] + ((Xi[1]-Xi[2])**2)/3))
  } else {
    Xi <- pred$xi.mod #This will be our comparaison point
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
                               partition = model$partition,
                               subdivision = new.subdivision)
      } else {
        new.partition <- model$partition
        new.subdivision <- model$subdivision
        new.partition[[model$nblock + 1]] <- option
        new.subdivision[[model$nblock + 1]] <- vector("list",1)
        new.subdivision[[model$nblock + 1]][[1]] <- c(0,1)
        model_update <- create(class = "lineqBAGP", x = xdesign, y = ydesign,
                               constrType = rep(constrType, (model$nblock+1)),
                               partition = new.partition,
                               subdivision = new.subdivision)
      }
    } else if (option_name=="case2"){#We're merging two blocks
      size_basis <- lapply(lapply(model$subdivision, function(j) sapply(j, function(k) length(k))),
                           function(j) prod(j))
      block_size <- lapply(model$partition, function(x) length(x))
      if ((size_basis[option[[1]]]==2 && size_basis[[option[2]]]==2) || option[1]==option[2] ||
          (block_size[[option[1]]]+block_size[[option[2]]])>Block_max_size){ 
        #Checking if the blocks are allowed to ber merged
        return(list(model, 0))
      } else{
        new.param <- merge_block(model$partition,model$subdivision, option[1],option[2])
        new.partition <- new.param[[1]]
        new.subdivision <- new.param[[2]]
        model_update <- create(class = "lineqBAGP", x = xdesign, y = ydesign,
                               constrType = model$constrType,
                               partition = new.partition,
                               subdivision = new.subdivision
        )
      }
    }
    criteria <- square_norm_int(model, model_update)
  }
  return(list(model_update, criteria))
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
#' @author A. F. Lopez-Lopera
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
  # pred <- predict(model, xtest = xtest)
  subdivision2 <- model$subdivision
  pos <- bijection(model$partition, choice[1])
  subdivision2[[pos[1]]][[pos[2]]] <- sort(c(model$subdivision[[pos[1]]][[pos[2]]],t))
  model2 <- create(class = "lineqBAGP", x = model$x, y = model$y,
                   constrType = model$constrType, 
                   partition = model$partition,
                   subdivision = subdivision2)
  return(-square_norm_int(model,model2))
}
