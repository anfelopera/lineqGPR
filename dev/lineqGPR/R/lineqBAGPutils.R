# Useful functions for \code{"lineqBAGP"} Models.

#' @title Block operation wrapper function
#' @description The additivity gives us block matrices, hence we can accelerate
#' the computations this is the role of algorithm in this section to facilitate
#' the manipulation of matrices of the type \code{list["vector"]} 
#' \eqn{\{1, \cdots, D \}} induces a natural bijection from
#' \eqn{\{ (j,k), j \in \{1, \cdots, J \}, k \in \{1,\cdots, \mathcal{J}_b\} \} }
#' into \eqn{\{1, \cdots, D\}} there is an straight one, \code{partition[[j]][[k]]}
#' gives the element \eqn{i} and ...  
#' 
#' @param A A list containing matrices.
#' @param type A string describing the type of operation to be applied for all the blocks.
#' 
#' @return A matrix.
#'
#' @author M. Deronzier and A. F. Lopez-Lopera
#' 
#' @importFrom Matrix bdiag
#' @export
block_to_matrix <- function(A, type = c("bdiag", "cbind", "rbind", "sum")) {
  type <- match.arg(type)
  if (type == "sum") {
    return (Reduce("+", A))
  } else {
    return (eval(parse(text = paste(type, "(", 
                                    paste("A[[", 1:length(A), "]]", sep = "", collapse = ","),
                                    ")", sep = ""))))
  }
}

# #' We can extract evaluation of our GP with a product of matrices $\Phi and \xi$, however, how to 
# #' arrange our function the tensorisation will give us a natural way to arrange things
# 
# #' This function gives the matrice $\Phi$ as a list of matrices per blocks 
# 
# #' Take the vector subdivision and $x$ and get out the object 
# #' $(phi_{{j,k}})(x_{j})_{j,k in}$
# #' @return Phi_per_var: such that Phi_per_var[[j]][[k]][i] give a seq
# #' of the hat function $\phi_{j,k,i}(X_i)$ where $(X_i)=(x^{(1)}_i,\cdots x^{(n)_i})$ and 
# #' \phi_{j,k,i} is the hat function given by the subdivision[[j]][[k]]  
# #' @return Phi_var_to_tensor: such that Phi_var_to_tensor[[j]] is the 
# #' tensorisation of vectors Phi_per_var[[j]][[k]] for k in $\{1,...,dim_block[j]\}$ 
# #' @return Big_Phi is the vector per block defined in the article

#' @title Hat basis functions per variable
#' @description Compute the hat basis functions for each input variable. 
#' 
#' @param subdivision A list containing the subdivisions for each variable:
#' \code{subdivision[[j]][[k]]} denotes the subdivision of the block \eqn{j} and
#' local input variable \eqn{k}.
#' @param partition A list containing the partition of set \eqn{\{1, \cdots, d\}}.
#' @param x An \eqn{n \times d} matrix with the input data.
#' 
#' @return A list containing the basis functions for each variable. The matrices
#' are indexed per blocks.
#'
#' @author M. Deronzier and A. F. Lopez-Lopera
#' 
#' @export
Phi_per_var <- function(subdivision, partition, x){
  dim_block <- sapply(subdivision, function(x) length(x))
  Phi.perVar <- vector("list", length(subdivision))
  for (j in 1:length(dim_block)) {
    Phi.perVar[[j]] <- lapply(1:dim_block[j],
                              function(k) basisCompute.lineqGP(x[ , partition[[j]][k]],
                                                               subdivision[[j]][[k]]))
  } 
  return(Phi.perVar)
}

#' @title Hat basis functions per block
#' @description Compute the hat basis functions for each block.
#' 
#' @param Phi.perVar A list containing the hat basis functions for each variable.
#' The matrices are indexed per blocks.
#' @return A list containing the hat basis functions for each block.
#'
#' @author M. Deronzier and A. F. Lopez-Lopera
#' 
#' @export
Phi_var_to_tensor <- function(Phi.perVar){
  # This function take the phi per variables and get back the family $\phi_{\el_j}(x_{j})$
  nblocks <- length(Phi.perVar)
  dim_block <- sapply(Phi.perVar, function(x) length(x))
  
  Phi.perBlock <- list("vector", nblocks)
  nbdesign <- nrow(Phi.perVar[[1]][[1]])
  for(j in 1:nblocks){
    Phi.perBlock[[j]] <- matrix(0, nrow = nbdesign, prod(sapply(Phi.perVar[[j]], ncol)))
    for (i in 1:nbdesign){
      Phi.perBlock[[j]][i, ] <- eval(parse(text = paste("Phi.perVar[[j]][[", 1:dim_block[j],
                                                        "]][", i, ", ]", sep = "", collapse = "%x%")))
    }
  }
  return(Phi.perBlock)
}

#' @title Hat basis functions per variable
#' @description Compute the hat basis functions for each input variable. 
#' 
#' @param subdivision A list containing the subdivisions for each variable:
#' \code{subdivision[[j]][[k]]} denotes the subdivision of the block \eqn{j} and
#' local input variable \eqn{k}.
#' @param partition A list containing the partition of set \eqn{\{1, \cdots, d\}}.
#' @param x An \eqn{n \times d} matrix with the input data.
#' 
#' \code{subdivision[[j]][[k]]} denotes the subdivision of the block \eqn{j} and
#' local input variable \eqn{k}.
#' 
#' @return A list containing the hat basis functions for each block.
#'
#' @author M. Deronzier and A. F. Lopez-Lopera
#' 
#' @export
Block_Phi<-function(subdivision, partition, x){
  return (Phi_var_to_tensor(Phi_per_var(subdivision, partition, x)))
}

########################### Creation of the kernel matrix ################################
#' @title Covariance matrix per variable
#' @description Compute the covariance matrix for each input variable.
#' 
#' @param subdivision A list containing the subdivisions for each variable:
#' \code{subdivision[[j]][[k]]} denotes the subdivision of the block \eqn{j} and
#' local input variable \eqn{k}.
#' @param par A list containing the values of the kernel parameters (variance, lengthscales).
#' \code{par[[j]][l]} denotes the \eqn{l}-th parameter of the block \eqn{j}.
#' @param type A character string corresponding to the type of the kernel.
#' Options: "gaussian", "matern32", "matern52", "exponential".
#' 
#' @return A list containing the covariance matrix for each variable. The matrices
#' are indexed per blocks.
#'
#' @author M. Deronzier and A. F. Lopez-Lopera
#' 
#' @export
Gamma_var <- function(subdivision, par, type) {
  nblocks <- length(subdivision)
  dim_block <- sapply(subdivision, function(x) length(x))
  
  Gamma.perVar <- vector("list", nblocks)
  names(Gamma.perVar) <- paste("block", 1:nblocks, sep = "")
  if (is.list(type)) {
    for (j in 1:nblocks) {
      Gamma.perVar[[j]] <- lapply(1:dim_block[j],
                                  function(k)
                                    kernCompute(subdivision[[j]][[k]],
                                                type = type[[j]][k],
                                                par = par[[j]][c(1, k+1)]))
      names(Gamma.perVar[[j]]) <- names(subdivision[[j]])
      attr(Gamma.perVar[[j]], "par") = par[[j]]
    }
  } else if (length(type) == length(subdivision)) {
    for (j in 1:nblocks) {
      Gamma.perVar[[j]] <- lapply(1:dim_block[j],
                                  function(k)
                                    kernCompute(subdivision[[j]][[k]],
                                                type = type[j],
                                                par = par[[j]][c(1, k+1)]))
      names(Gamma.perVar[[j]]) <- names(subdivision[[j]])
      attr(Gamma.perVar[[j]], "par") = par[[j]]
    }
  } else if (length(type) == 1) {
    for (j in 1:nblocks) {
      Gamma.perVar[[j]] <- lapply(1:dim_block[j],
                                  function(k)
                                    kernCompute(subdivision[[j]][[k]],
                                                type = type,
                                                par = par[[j]][c(1, k+1)]))
      names(Gamma.perVar[[j]]) <- names(subdivision[[j]])
      attr(Gamma.perVar[[j]], "par") = par[[j]]
    }
  } else {
    stop("The argument 'type' must be a list with the same structure of the
         'subdivision', a sequence with legth(type) = length(subdivision), or
         a character string")
  }
  names(Gamma.perVar) <- paste("block", 1:length(subdivision), sep = "")
  return(Gamma.perVar)
}

#' @title Covariance matrix per block
#' @description Compute the covariance matrix for each block.
#' 
#' @param Gamma.perVar A list containing the covariance matrix for each variable.
#' The matrices are indexed per blocks.
#' @return A list containing the covariance matrix for each block.
#'
#' @author M. Deronzier and A. F. Lopez-Lopera
#' 
#' @export
Gamma_var_to_tensor <- function(Gamma.perVar) {#, dim_block, J){
  nblocks <- length(Gamma.perVar)
  dim_block <- sapply(Gamma.perVar, function(x) length(x))
  Gamma.perBlock <- list("vector", nblocks)
  names(Gamma.perBlock) <- paste("block", 1:nblocks, sep = "")
  for(j in 1:nblocks) {
    params <- attr(Gamma.perVar[[j]], "par")
    Gamma.perBlock[[j]] <- eval(parse(text = paste("(Gamma.perVar[[j]][[",
                                                   1:dim_block[j],
                                                   "]]/params[1])", sep = "", collapse = "%x%")))
    Gamma.perBlock[[j]] <- params[1]*Gamma.perBlock[[j]]
    
    if (dim_block[[j]] > 1) {
      expr <- matrix(paste("Gamma.perVar[[j]][[", seq(dim_block[[j]]),"]]", sep = ""),
                     dim_block[[j]], dim_block[[j]], byrow = TRUE)
      diag(expr) <- paste("attr(Gamma.perVar[[j]][[", seq(dim_block[[j]]),
                          "]], 'gradient')[[2]]", sep = "")
      for (i in seq(dim_block[[j]]+1)) {
        if (i == 1) {
          gradGammaTemp <- Gamma.perBlock[[j]]/params[1]
          attr(Gamma.perBlock[[j]], "gradient")$sigma2 <- gradGammaTemp
        } else {
          expr2 <- paste(expr[i-1, ], collapse = " %x% ")
          gradGammaTemp <- eval(parse(text = expr2))
          attr(Gamma.perBlock[[j]], "gradient")[[i]] <- gradGammaTemp/params[1]^(dim_block[[j]]-1)
        }
      }
      names(attr(Gamma.perBlock[[j]], "gradient")) <- c("sigma2", paste("theta", seq(dim_block[[j]]), sep = ""))
    }
    attr(Gamma.perBlock[[j]], "derivative") <- NULL
    attr(Gamma.perBlock[[j]], "par") <- params
  }
    
    
  
  
  
  
  return(Gamma.perBlock)
}

#### some basic operations for block lists ####

#' @title Computation of elementwise operations of lists of matrices
#' @description Compute an elementwise operation to a list (or two lists) of matrices.
#' 
#' @param A A list containing matrices.
#' @param operation A character string with the operation to be applied.
#' Options: "sum", "prod", "transpose", "chol", "inv", "chol2inv", "scalarMatMul".
#' @param B A list containing matrices.
#' @param alpha A real number.
#' @param beta A real number.
#' @return A list after applying the elementwise operation.
#'
#' @author M. Deronzier and A. F. Lopez-Lopera
#' 
#' @export
block_compute <- function(A,
                          operation = c("sum", "prod", "transpose",
                                        "chol", "inv", "chol2inv",
                                        "scalarMatMul"),
                          B = NULL, alpha = 1, beta = 1) {
  operation <- match.arg(operation)
  switch (operation,
    sum = {
      return(lapply(1:length(A), function(x) alpha*A[[x]] + beta*B[[x]]))
    }, prod = {
      return (lapply(1:length(A), function(x) A[[x]] %*% B[[x]]))
    }, transpose = {
      return(lapply(A, function(x) t(x)))
    }, chol = {
      return (lapply(A, function(x) chol(x)))
    }, chol2inv = {
      return (lapply(A, function(x) chol2inv(x)))
    }, inv = {
      return (lapply(A, function(x) chol2inv(chol(x))))
    }, scalarMatMul = {
      return(lapply(A, function(x) alpha*x)) 
    }
  )
}

#### backup of unused fonctions ####

#' @title Bijection operator
#' @description A wrapper function used in \code{"lineqBAGP"} models to obtain the indices
#' of the corresponding block and the local input variable of a global variable. 
#' The partition \eqn{\{\mathcal{J}_1, \cdots, \mathcal{J}_J \}} of
#' \eqn{\{1, \cdots, D \}} induces a natural bijection from
#' \eqn{\{ (j,k), j \in \{1, \cdots, J \}, k \in \{1,\cdots, \mathcal{J}_b\} \} }
#' into \eqn{\{1, \cdots, D\}} there is an straight one, \code{partition[[j]][[k]]}
#' gives the element \eqn{i} and ... (anfe: to be complete!!)  
#' 
#' @param partition A list containing the partion of set \eqn{\{1, \cdots, D\}}
#' @param i The targeted global input variable 
#' 
#' @return a sequence containing the indices of the block and the input variable. 
#'
#' @author M. Deronzier
#'
# #' @examples
# #' partition <- list(c(1,3), c(2,4,5))
# #' bijection(partition, 3)
#' 
# #' @export
bijection <- function(partition, i){
  for (j in 1:length(partition)){
    if (i %in% partition[[j]]) 
      return(c(j, which(partition[[j]] == i)))
  }
}

#' @title Name
#' @description description  
#' 
#' @param A A list containing matrices
#' @param B A list containing matrices
#' 
#' @return A matrix object 
#'
#' @author M. Deronzier and A. F. Lopez-Lopera
#' 
# #' @export
block_scalar <- function(A, B){
  scalar <- eval(parse(text = paste("cbind(",
                                    paste("A[[1]]%*%B[[", 1:length(B), "]]",
                                          sep = "", collapse = ","),
                                    ")", sep = "")))
  for (i in 2:length(A)){
    scalar <- rbind(scalar, eval(parse(text = paste("cbind(",
                                                    paste("A[[", i, "]]%*%B[[", 1:length(B), "]]",
                                                          sep = "", collapse = ","),
                                                    ")", sep = ""))))
  }
  return(scalar)
}

#' @title Matrix Inversion Lemma
#' @description Compute \eqn{(Z + U W V^\top)^{-1}} via the matrix Inversion Lemma
#' (Rasmussen and Williams, 2015; Appendix A.3).
#' 
#' @param U An \eqn{n \times m} matrix 
#' @param V An \eqn{n \times m} matrix 
#' @param W An \eqn{m \times m} inversible matrix 
#' @param invZ An \eqn{n \times n} inversible matrix 
#' 
#' @return \eqn{(Z + U W V^\top)^{-1}}
#'
#' @author M. Deronzier and A. F. Lopez-Lopera
#' 
#' @references C. E. Rasmussen and C. K. I. Williams (2005),
#' "Gaussian Processes for Machine Learning".
#' \emph{The MIT Press}
#' 
# #' @export
inv_lemma <- function(U, V, W, invZ){
  # Computation of the conditional covariance matrice using the inversion lemma formulation:
  # $(Z+UWV^{\top})^{-1}=Z^{-1}-Z^{-1}U(W^{-1}+ V^{t}Z^{-1}U)^{-1}V^{\top}Z^{-1}$ it is useful
  # when $U, V\in M_{n,m}, W \in M_{m,m} $ and $Z \in M_{n,n}$. the factorisation form is useful
  # when $m<<n$
  # 
  # As it is useful only if we know invZ we take it in parameters 
  # !!!!!! Carefull works only with W symetric matrix !!!!!!! inversion cholesky  !!!!!!!!! 
  
  inv_W <- chol2inv(chol(W))
  return(invZ - invZ %*% U %*% chol2inv(chol(inv_W%*%t(V)%*%invZ%*%U)) %*% t(V) %*% invZ)
}

# # Block case, be careful about matrices we use
# block_inv_lemma <- function(W, invZ, U, V = NULL){
#   if(is.null(V)){
#     inv_ZU <- inv_ZV <- block_compute(invZ, "prod", U)
#   }
#   else{
#     inv_ZU <- block_compute(invZ, U, "prod")
#     inv_ZV <- block_compute(invZ, V, "prod")
#   }
#     t_inv_ZV <- block_compute(inv_ZV, "transpose")
#   return(block_to_matrix(invZ, "diag") - block_to_matrix(WU, "cbind")%*% ## %-%?
#          chol2inv(chol(block_to_matrix(block_compute(t_inv_ZV, "prod", U)), "sum"))%*%
#          block_to_matrix(t_inv_ZV, "rbind"))
# }
