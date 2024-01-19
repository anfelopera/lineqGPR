# Useful functions for \code{"lineqBAGP"} Models.
# If there an issure remove .Rproj.user

#' @title Block operation wrapper function
#' @description The additivity gives us block matrices, hence we can accelerate
#' the computations this is the role of algorithm in this section to facilitate
#' the manipulation of matrices of the type \code{list["vector"]} 
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

#' @title Block operation wrapper function
#' @description The additivity gives us block matrices, hence we can accelerate
#' the computations this is the role of algorithm in this section to facilitate
#' the manipulation of matrices of the type \code{list["vector"]} 
#' 
#' @param M a matrix
#' @param subdivision a subdivision of [0,1]
#' @param type is a string defining how we want to cut the blocks
#' 
#' @return A matrix per block, meaning a list of matrices
#'
#' @author M. Deronzier and A. F. Lopez-Lopera
#' 
#' @importFrom Matrix bdiag
#' @export

matrix_to_block <- function(M, subdivision, type = c("bdiag", "cbind", "rbind")) {
  J <- length(subdivision)
  block_size <- c(0, sapply(lapply(subdivision, function(x) sapply(x,function(y) length(y))),
                       function(x) prod(x)))
  block_M <- vector("list", length(subdivision))
  if (type=="diag"){
    for (j in 1:J){
      n1 <- sum(block_size[1:j])+1
      n2 <- n1+ block_size[j+1]-1
      block_M[[j]] <- M[n1:n2,n1:n2]
    }
  } else if (type=="cbind"){
    for (j in 1:J){
      n1 <- sum(block_size[1:j])+1
      n2 <- n1+ block_size[j+1]-1
      block_M[[j]] <- as.matrix(matrix(M[,n1:n2], ncol=1))
    }
  } else if (type=="rbind"){
    for (j in 1:J){
      n1 <- sum(block_size[1:j])+1
      n2 <- n1+ block_size[j+1]-1
      block_M[[j]] <- as.matrix(M[n1:n2,])
    }
  }
  return(block_M)
}
# #' We can extract evaluation of our GP with a product of matrices Phi and xi, however, how to 
# #' arrange our function the tensorisation will give us a natural way to arrange things
# 
# #' This function gives the matrice Phi as a list of matrices per blocks 
# 
# #' Take the vector subdivision and x and get out the object 
# #' (phi_{{j,k}})(x_{j})_{j,k in}
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
  nblocks <- length(Phi.perVar)
  dim_block <- sapply(Phi.perVar, function(x) length(x))
  
  Phi.perBlock <- vector("list", nblocks)
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
  Gamma.perBlock <- vector("list", nblocks)
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
                                        "scalarMatMul", "det"),
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
      return (lapply(A, function(x) chol(x+1e-9*diag(ncol(x)))))#########################Potential instability issue Here
    }, chol2inv = {
      return (lapply(A, function(x) chol2inv(x)))#+1e-9*diag(ncol(x)))))
    }, inv = {
      return (lapply(A, function(x) chol2inv(chol(x))))#+1e-8*diag(ncol(x))))))
    }, scalarMatMul = {
      return(lapply(A, function(x) alpha*x)) 
    }, det = {
      return (sapply(A, function(x) det(x)))
    }
  )
}

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
#' @author M. Deronzier and A. F. Lopez-Lopera
#'
#'@examples
#' partition <- list(c(1,3), c(2,4,5))
#' bijection(partition, 3)
#' @export 

bijection <- function(partition, i){
  for (j in 1:length(partition)){
    if (i %in% partition[[j]]) 
      return(c(j, which(partition[[j]] == i)))
  }
}

#' @title Block product matrices
#' @description Algorithm returning for two blocks matrices A, B the matrice per block A[i]*B[j] 
#' 
#' @param A A list containing matrices
#' @param B A list containing matrices
#' 
#' @return A matrix object 
#'
#' @author M. Deronzier and A. F. Lopez-Lopera
#' 
#' @export
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


#' @title Algorithm that check if you can make a changment basis.
#' @description This algorithm check if we are in the hypotheses of the proposition
#' to make a changment basis. 
#' 
#' @param partition1,partition2 two partitions  
#' @param subdivision1,subdivision2 two subdivisions adapted to the partitions we need 
#' 
#' @return boolean indicating if we can make a changment basis
#'
#' @author M. Deronzier 
#' 
#' @export

condition_changment_basis <- function(partition1,partition2,subdivision1,subdivision2){
  set1 <- eval(parse(text = paste("sort(c(",
                                   paste("c(partition1[[", 1:length(partition1), "]])",
                                         sep = "", collapse = ","),
                                   "))", sep = "")))
  set2 <- eval(parse(text = paste("sort(c(",
                                  paste("c(partition2[[", 1:length(partition2), "]])",
                                        sep = "", collapse = ","),
                                  "))", sep = "")))
  if(isFALSE(setequal(set1, intersect(set1,set2)))){
    return(FALSE)
  }else{
    for (i in set1){
      bij1 <- bijection(partition1,i)
      bij2 <- bijection(partition2,i)
      if(isFALSE(setequal(subdivision1[[bij1[1]]][[bij1[2]]], 
                          intersect(subdivision1[[bij1[1]]][[bij1[2]]],
                          subdivision2[[bij2[1]]][[bij2[2]]])))){
        return(FALSE)
      }
    }
  }
  return(TRUE)
}

#' @title Creation of changing basis matrices in case 1 of the article
#' @description This algorithm make the changement basis corresponding of the case 1
#' of the article.
#' 
#' @param subdivision1 a subdivision
#' @param subdivision2 a subdivision
#' 
#' @return Boolean indicating if we can make a changment basis
#'
#' @author M. Deronzier 
#' 
#' @export

changment_basis_knots <- function(subdivision1,subdivision2){
  n <- length(subdivision1)
  submatrix <- vector("list", length(subdivision2))
  for (j in 1:length(subdivision1)){
    submatrix[[j]] <- basisCompute.lineqGP(subdivision2[[j]], subdivision1[[j]])
  } 
  return(eval(parse(text = paste("submatrix[[", 1:n,"]]", sep = "", collapse = "%x%"))))
}


#' @title Creation of changing basis matrices in case 2 of the article
#' @description This algorithm make the changment basis corresponding to the case 2
#' of the article.
#' 
#' @param partition1 a partition 
#' @param partition2 a partition
#' @param subdivision2  a subdivision 
#' 
#' @return Boolean indicating if we can make a changment basis
#'
#' @author M. Deronzier 
#' 
#' @export

changment_basis_inclusion <- function(partition1, partition2, subdivision2){
  n1 <- length(partition1)
  n2 <- length(partition2)
  submatrix <- vector("list", n1)
  for(i in 1:n2){
    if (inclusion(partition2[i],partition1)){
      submatrix[[i]] <- diag(length(subdivision2[[i]]))
    } else {
      submatrix[[i]] <- rep(1,length(subdivision2[[i]]))
    }
  }
  return(eval(parse(text = 
                      paste("submatrix[[", 1:n2,"]]", collapse = "%x%", sep = "")
  )))
}

#' @title Indicate if set1 is in set2
#' @description Give a boolean saying if set1 is in set2
#' 
#' @param set1,set2 two vectors representing sets  
#' 
#' @return TRUE if set1 is included in set2  FALSE otherwise
#'
#' @author M. Deronzier 
#' 
#' @export

inclusion <- function(set1,set2){
  return (setequal(set1,intersect(set1,set2)))
}


#' @title Creation of changing basis matrices in the general case
#' @description The algorithm explained in the proof of the article.
#' 
#' @param partition1 a partition 
#' @param partition2 a partition
#' @param subdivision1 a subdivision
#' @param subdivision2  a subdivision 
#' 
#' @return P a changment basis matrix from the basis generated by (partition1, subdivision1) to 
#' the basis generated by (partition2, subdivision2)
#'
#' @author M. Deronzier 
#' 
#' @export

changement_basis_matrix <- function(partition1,partition2,subdivision1,subdivision2){
  if (isFALSE(condition_changment_basis(partition1,partition2,subdivision1,subdivision2))){
    stop("Changement base can't be computed")
  }
  
  # For each block in partition1, the block in which it is included in partition2
  container <- sapply(partition1, function(X) which(sapply(partition2, function(x) inclusion(X, x))==TRUE))
  # Useful variables
  
  n1 <- length(partition1)
  n2 <- length(partition2)
  
  subdivision_size1 <- lapply(subdivision1, function(x) sapply(x,function(y) length(y)))
  subdivision_size2 <- lapply(subdivision2, function(x) sapply(x,function(y) length(y)))
  block_size1 <- sapply(subdivision_size1, function(x) prod(x))
  block_size2 <- sapply(subdivision_size2, function(x) prod(x))
  
  # The first modification we'll do corresponding to case 1
  P1_block <- vector("list", n1)
  new_subdivision <- vector("list",n1)
  for (j in 1:n1){
    new_subdivision[[j]] <- lapply(partition1[[j]], function(x) 
      subdivision2[[container[j]]][[which(x==partition2[[container[j]]])]])
    P1_block[[j]] <- changment_basis_knots(subdivision1[[j]],new_subdivision[[j]])
  }
  P2_block <- vector("list", n1)
  for (j in 1:n1){
    P2_block[[j]] <- changment_basis_inclusion(partition1[[j]], partition2[[container[[j]]]],
                                                subdivision2[[container[[j]]]])
  }
  P_block <- lapply(1:n1, function(x) P2_block[[x]]%*%P1_block[[x]])
  
  if(n2==1){
    return(block_to_matrix(P_block, "cbind"))
  }
  #### Algorithm to agglomerate blocks 
  P <- vector("list", n1)
  for(j in 1:n1){
    if(container[j]==1){
      P[[j]] <- rbind(P_block[[j]], matrix(0, ncol= block_size1[j], nrow=(sum(block_size2[2:n2]))))
    } else if (container[j]==n2){
      P[[j]] <- rbind(matrix(0, ncol=block_size1[j], nrow=(sum(block_size2[1:(n2-1)]))), P_block[[j]])
    } else {
      P[[j]] <- rbind(matrix(0, ncol=block_size1[j], nrow=(sum(block_size2[1:(container[j]-1)]))), P_block[[j]],
                      matrix(0, ncol=block_size1[j], nrow=(sum(block_size2[(container[j]+1):n2]))))
    }
  }
  return(block_to_matrix(P,"cbind"))
}

############################## Function to calculate the square norm criteria ##########################################

#' @title integral in the 1 dimension case
#' @description gives the value of integral of  the hat basis functions constructed 
#' from a 1 dimensional subdivision.
#'
#' @param subdivision of type sequence, representing a subdivision of [0,1]  
#' 
#' @return the integral fo each hat function described by the subdivision
#'
#' @author M. Deronzier 
#' 
#' @export

seq_int <- function(subdivision){
  J <- length(subdivision)
  subdivision1 <- c(0, subdivision, 1)
  1/2*sapply(2:(J+1), function(i) subdivision1[i+1]-subdivision1[i-1])
}


#' @title integral in the general case
#' @description Give the vector E defined in the article.
#'
#' @param subdivision of type list(list(sequence)) the subdivisions of [0,1] 
#' @return the integral fo each hat function described by the subdivision
#'
#' @author M. Deronzier 
#' 
#' @export

Mat_E <- function(subdivision){
  J <- length(subdivision)
  E <- vector("list", J)
  for (j in 1:J){
    E[[j]] <- eval(parse(text = 
              paste("seq_int(subdivision[[",j,"]][[", 1:length(subdivision[[j]]),"]])", collapse = "%x%", sep = "")
    ))
  }
  return(E)
}


#' @title Gram matrix for the 1 dimensional basis functions (\code{"lineqBAGP"})
#' @description Compute the Gram matrix of the basis functions for \code{"lineqBAGP"} models
#' 
#' @param subdivision a subdivision of [0,1]
#' @param n the size of the subdivision
#' 
#' @return Gram matrix of the basis functions
#'
#' @author M. Deronzier 
#'
#' @references F. Bachoc, A. F. Lopez-Lopera, and O. Roustant (2020),
#' "Sequential construction and dimension reduction of Gaussian processes under inequality constraints".
#' \emph{ArXiv e-prints}
#' <arXiv:2009.04188>
#'
#' @export

sub_Gram <- function(subdivision,n=length(subdivision)){
  n <- length(subdivision)
  sub <- c(0,subdivision,1)
  dsub <- diff(sub)
  mat <- diag(1/3*(dsub[1:n]+dsub[2:(n+1)]))
  for (i in 1:(n-1)){
    mat[i,i+1] <- mat[i+1,i] <- (1/6)*dsub[i+1]
  }
  return(mat)
}

#' @title Gram matrix for the 1 dimensional basis functions (\code{"lineqBAGP"})
#' @description Compute the Gram matrix of the basis functions for \code{"lineqBAGP"} models
#' 
#' @param subdivision the multiples subdivisions on the variable
#' @param J number of blocks
#' @return Gram matrix of the basis functions
#'
#' @author M. Deronzier 
#'
#' @references F. Bachoc, A. F. Lopez-Lopera, and O. Roustant (2020),
#' "Sequential construction and dimension reduction of Gaussian processes under inequality constraints".
#' \emph{ArXiv e-prints}
#' <arXiv:2009.04188>
#'
#' @export

Gram <- function(subdivision, J =length(subdivision)){
  gram <- vector("list", J)
  for(j in 1:J){
    gram[[j]] <- eval(parse(text = 
                              paste("sub_Gram(subdivision[[",j,"]][[",
                                    1:length(subdivision[[j]]),"]])",
                                    collapse = "%x%", sep = "")
    ))
  }  
  return(gram)
}


#' @title Gram matrix of the basis functions (\code{"lineqBAGP"})
#' @description Compute the Gram matrix of the basis functions for \code{"lineqBAGP"} models
#' 
#' @param model1 a model of BAGP
#' @param model2 a more precise model of BAGP
#' 
#' @return Gram matrix of the basis functions
#'
#' @author M. Deronzier 
#'
#' @references F. Bachoc, A. F. Lopez-Lopera, and O. Roustant (2020),
#' "Sequential construction and dimension reduction of Gaussian processes under inequality constraints".
#' \emph{ArXiv e-prints}
#' <arXiv:2009.04188>
#'
#' @export

square_norm_int <- function(model1,model2){
  Gram_block <- Gram(model2$subdivision)
  E_block <- Mat_E(model2$subdivision)
  mat <-changement_basis_matrix(model1$partition, model2$partition, model1$subdivision, model2$subdivision)
  Xi <- matrix_to_block((mat%*%predict(model1,0)$xi.mod),
                        model2$subdivision, type="rbind")
  new.Xi <- matrix_to_block(predict(model2,0)$xi.mod, model2$subdivision, type = "rbind")
  eta <- block_compute(Xi,"sum", new.Xi,1,-1)
  criteria <- block_to_matrix(block_compute(block_compute(block_compute(
    eta,"transpose"),"prod", Gram_block),"prod",eta),"sum")
  products <- as.matrix(unlist(block_compute(block_compute(eta,"transpose"),"prod",E_block)))
  if (length(model1$partition)>length(model2$partition)){
    return((criteria + sum(products%x%products) - sum(sapply(products, function(x) x^2)))
           /(model2$nknots-model1$nknots+2)^(1.4))
  } else{
    return((criteria + sum(products%x%products) - sum(sapply(products, function(x) x^2)))
           /(model2$nknots-model1$nknots)^(1.4))
  }
}

#' @title Gram matrix of the basis functions (\code{"lineqBAGP"})
#' @description Compute the Gram matrix of the basis functions for \code{"lineqBAGP"} models
#' 
#' @param model1 a model of BAGP
#' @param model2 a more precise model of BAGP
#' @param alpha a float indicating the relation between the two criteria
#' 
#' @return Gram matrix of the basis functions
#'
#' @author M. Deronzier 
#'
#' @references F. Bachoc, A. F. Lopez-Lopera, and O. Roustant (2020),
#' "Sequential construction and dimension reduction of Gaussian processes under inequality constraints".
#' \emph{ArXiv e-prints}
#' <arXiv:2009.04188>
#'
#' @export

square_norm_int_2 <- function(model1,model2, alpha=0.2){
  Gram_block <- Gram(model2$subdivision)
  E_block <- Mat_E(model2$subdivision)
  mat <- changement_basis_matrix(model1$partition, model2$partition, model1$subdivision, model2$subdivision)
  Xi <- matrix_to_block((mat%*%predict(model1,0)$xi.mod), model2$subdivision, type="rbind")
  new.Xi <- matrix_to_block(predict(model2,0)$xi.mod, model2$subdivision, type = "rbind")
  eta <- block_compute(Xi,"sum", new.Xi,1,-1)
  criteria <- block_to_matrix(block_compute(block_compute(block_compute(
    eta,"transpose"),"prod", Gram_block),"prod",eta),"sum")
  criteria1 <- block_to_matrix(block_compute(block_compute(block_compute(
    Xi,"transpose"),"prod", Gram_block),"prod",Xi),"sum")
  criteria2 <- block_to_matrix(block_compute(block_compute(block_compute(
    new.Xi,"transpose"),"prod", Gram_block),"prod",new.Xi),"sum")
  products <- as.matrix(unlist(block_compute(block_compute(eta,"transpose"), "prod", E_block)))
  products1 <- as.matrix(unlist(block_compute(block_compute(Xi,"transpose"), "prod", E_block)))
  products2 <- as.matrix(unlist(block_compute(block_compute(new.Xi,"transpose"), "prod", E_block)))
  pred <- predict(model2, model2$x)$y.mod
  diff_norm <- norm(pred-model2$y)/(norm(pred)+norm(model2$y))
  globcrit <- (criteria + sum(products%x%products) - sum(sapply(products, function(x) x^2)))/
              (criteria1 + sum(products1%x%products1) - sum(sapply(products1, function(x) x^2))+
               criteria2 + sum(products2%x%products2) - sum(sapply(products2, function(x) x^2)))
  res <- alpha*(1-diff_norm)+globcrit
  if (length(model1$partition)>length(model2$partition)){
    return(res/(model2$nknots-model1$nknots+1)^1)
  } else{
    return(res/(model2$nknots-model1$nknots)^1)
  }
}

#' @title Merge block (\code{"lineqBAGP"})
#' @description take a partition and subdivision and  the numero of blocks to merge
#' 
#' 
#' @param partition a partition
#' @param subdivision a subdivision
#' @param kernParam list of kernel parameters
#' @param n1 an integer corresponding to the first block to merge
#' @param n2 an integer corresponding to the second block to merge
#' 
#' @return partition and subdivision of the merged model
#'
#' @author M. Deronzier 
#'
#' @references F. Bachoc, A. F. Lopez-Lopera, and O. Roustant (2020),
#' "Sequential construction and dimension reduction of Gaussian processes under inequality constraints".
#' \emph{ArXiv e-prints}
#' <arXiv:2009.04188>
#'
#' @export

merge_block <- function(partition, subdivision, kernParam, n1, n2){
  partition1 <- vector("list", length(partition)-1)
  subdivision1 <- vector("list", length(partition)-1)
  kernParam1 <- vector("list",2)
  names(kernParam1) <- names(kernParam)[1:2]
  kernParam1$type <- kernParam$type
  kernParam1$par <- vector("list", length(partition)-1)
  if(length(partition)>2){
    for (i in 2:(length(partition)-1)){
      partition1[[i]] <- partition[[setdiff(1:length(partition),c(n1,n2))[i-1]]]
      subdivision1[[i]] <- subdivision[[setdiff(1:length(partition),c(n1,n2))[i-1]]]
      kernParam1$par[[i]] <- kernParam$par[[setdiff(1:length(partition),c(n1,n2))[i-1]]]
    }
  }
  partition1[[1]] <- c(partition[[n1]],partition[[n2]])
  subdivision1[[1]] <- c(subdivision[[n1]],subdivision[[n2]])
  kernParam1$par[[1]] <- c(kernParam$par[[n1]][1]+kernParam$par[[n2]][1],
                           kernParam$par[[n1]][-1], kernParam$par[[n2]][-1])
  new.partition <- lapply(partition1, function(x) sort(x))
  new.subdivision <- subdivision1
  new.kernParam <- kernParam
  for (i in 1:length(partition1[[1]])){
    pos<-bijection(partition1, new.partition[[1]][i])
    new.subdivision[[1]][[i]]<- subdivision1[[pos[1]]][[pos[2]]]
    new.kernParam$par[[1]][i+1] <- kernParam1$par[[pos[1]]][pos[2]+1]
  }
  return(list(new.partition, new.subdivision, new.kernParam))
}



#' @title Rename
#' @description rename the parameters of the model  
#' 
#' 
#' @param model a model
#' 
#' @return the model with the variables renamed
#'
#' @author M. Deronzier 
#'
#' @references F. Bachoc, A. F. Lopez-Lopera, and O. Roustant (2020),
#' "Sequential construction and dimension reduction of Gaussian processes under inequality constraints".
#' \emph{ArXiv e-prints}
#' <arXiv:2009.04188>
#'
#' @export

rename <- function(model){
  #names(model$subdivision) <- names(model$partition) <- paste("block", 1:length(model$partition), sep = "")
  for (j in 1:length(model$partition)) {
    names(model$partition[[j]]) <- names(model$subdivision[[j]]) <- paste("var", model$partition[[j]], sep = "")
    names(model$kernParam$par[[j]]) <- c(paste("sigma", j, sep = ""), 
                                         paste("var", model$partition[[j]], sep = ""))
    
  }
  return(model)
}

################################################ PLOT FUNCTIONS ################################################

#' @title function
#' @description rename the parameters of the model  
#' 
#' @param x list of history
#' 
#' @return The choices that we make for iterations
#'
#' @author M. Deronzier 
#'
#' @references F. Bachoc, A. F. Lopez-Lopera, and O. Roustant (2020),
#' "Sequential construction and dimension reduction of Gaussian processes under inequality constraints".
#' \emph{ArXiv e-prints}
#' <arXiv:2009.04188>
#'
#' @export

f <- function(x){
  if (length(x)==1){
    return(paste("(",x,")", sep=""))
  }else{
    return(paste("[", x[1], ", ", x[2],"]", sep=""))
  }
}

#' @title Plot history
#' @description Function that plot the evolution of the criteria over iterations   
#' 
#' @param res a list containing the model and the history of choices returned by MaxMod
#' @param i an integer 
#' 
#' @return the model with the variables renamed
#'
#' @author M. Deronzier 
#'
#' @references F. Bachoc, A. F. Lopez-Lopera, and O. Roustant (2020),
#' "Sequential construction and dimension reduction of Gaussian processes under inequality constraints".
#' \emph{ArXiv e-prints}
#' <arXiv:2009.04188>
#'
#' @export
plot_history <- function(res, i){
  history <- lapply(res$history, function(x) f(x))
  n <- length(res$history)
  iteration <- 1:(n-1)
  model <- res$model
  hist <- history[2:n]
  histC <- res$hist_Criteria[2:n]
  histR <- res$hist_diffnorm[2:n]
  histCR <- histC/sqrt(histC)
  choices <- res$history
  df <- data.frame(iteration, histC,histR)
  scaleFac <- max(histC)/max(histR)
  
  ggplot(df, aes(iteration, histC))+
    geom_line(aes(y = histC), color = "red")+
    geom_line(aes(y = histR*scaleFac), color = "blue")+
    geom_label(aes(x = iteration, y = histC, label = hist))+
    scale_y_continuous(name="L2 Square diff", sec.axis=sec_axis(~./scaleFac, name="R-square")) +
    theme(
      axis.title.y.left=element_text(color="red"),
      axis.text.y.left=element_text(color="red"),
      axis.title.y.right=element_text(color="blue"),
      axis.text.y.right=element_text(color="blue"))+
    ggtitle(paste("Evolution of criteria for function ",i, sep=""))
}
