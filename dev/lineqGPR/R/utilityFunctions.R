#' @title Useful functions for \code{"lineqGP"} Models.
#' @description Manipulation of blocks matrices for fast computation \code{"lineqGP"} models.
#' 
#' @param sub_partition :list[seq] A sub-partition for first input and locations 
#' @param subdivision :list[list[seq]] The subdivisions for each axe subdivision[[j]][[k]]
#' @param subdivision_taille list[seq] subdivision_size[[j]][k] =length(subdivision[[j]][[k]])
#' @param dim_block :seq, the size of the blocks, dim_block[j]=length(partition[[j]])
#'  will give the subdivision of $[0,1]$ for the variable partition[[j]][[k]]    
#' @param X matrix such that X[,i] is an observation $x^{(i)}$
#' @param x float supposed to be a real number in $[0,1]$
#' @param Y the vector of observation Y[i]=f(X[,i]) 
#' @param A,B list["vector"] Matrices per block of type 
#' @param d a number corresponding to the dimension of the input space.
#' 
#' @return Kernel matrix \eqn{K(x_1,x_2)}{K(x1,x2)}
#' (or \eqn{K(x_1,x_1)}{K(x1,x1)} if \eqn{x_2}{x2} is not defined).
#'
#' @author A. F. Lopez-Lopera, M. Deronzier
#'
#' @examples
#' x <- seq(0, 1, 0.01)
#' K <- kernCompute(x, type = "gaussian", par =  c(1, 0.1))
#' image(K, main = "covariance matrix")
#'
#' @seealso \code{\link{k1gaussian}}, \code{\link{k1matern52}}, \code{\link{k1matern32}}, \code{\link{k1exponential}} 
#' 
#' @export

#' The partition $\{\mathcal{J}_1, \cdots, \,\mathcal{J}_J\} of $\{1, \cdots,\, D\}$ induces a
#' natural bijection from $\{(j,k), j \in \{1, \cdots, J \},\, k \in \{1,\cdots, \,|\mathcal{J}_b\}\|\}\}$
#' into \{1, \cdots, D\} there is an straight one, partition[[j]][[k]] gives the element i and 
#' @return Bijection will give the other one 

Bijection <- function(partition, i){
  for (j in 1:length(partition)){
    if (i %in% partition[[j]]){
      return (c(j,which(partition[[j]]==i)))
    }
  }
}

#' The additivity gives us block matrices, hence we can accelerate the computations this is the role of
#' algorithm in this section to facilitate the manipulation of matrices of the type list["vector"]

#' The following functions are meant to make matricial computations with matrices per blocks. 
#' @return block_to_vect input: matrix represented per blocks output: matrix


block_to_vect <- function(A, type="diagonal"){
  if (type=="diagonal"){
    return (eval(parse(text = paste("bdiag(",
                                    paste("A[[", 1:length(A), "]]", sep = "", collapse = ","),
                                    ")", sep = ""))))
  }
  else if (type=="line"){
    return (eval(parse(text = paste("cbind(",
                                    paste("A[[", 1:length(A), "]]", sep = "", collapse = ","),
                                    ")", sep = ""))))
  }
  else if (type=="column"){
    return (eval(parse(text = paste("rbind(",
                                    paste("A[[", 1:length(A), "]]", sep = "", collapse = ","),
                                    ")", sep = ""))))
  }
  else if (type=="sum"){
    return (Reduce("+",A))
  }
  else{
    stop("type invalid : type has to be in ['diagonal','line', 'column','sum']")
  }
}


block_scalar <- function(A,B){
  scalar <- eval(parse(text = paste("cbind(",
                                    paste("A[[1]]%*%B[[",1:length(B), "]]", sep = "", collapse = ","),
                                    ")", sep = "")))
  for (i in 2:length(A)){
    scalar <- rbind(scalar,eval(parse(text = paste("cbind(",
                                                   paste("A[[",i,"]]%*%B[[",1:length(B), "]]", sep = "", collapse = ","),
                                                   ")", sep = ""))))
  }
  return(scalar)
}


#' Here are some basic operations on the block matrices
#' @return block_operation input: matrix represented per blocks output:operation on matrix

  
block_prod <- function(A,B)
  return (lapply(1:length(A), function(x) A[[x]]%*%B[[x]]))

block_chol <- function(A)
  return (lapply(A, function(x) chol(x)))

block_chol2inv <- function(A)
  return (lapply(A, function(x) chol2inv(x)))

block_inv <- function(A)
  return(block_chol2inv(block_chol(A)))

block_transpose <- function(A)
  return(lapply(A, function(x) t(x)))

block_sum <- function(A,B,alpha=1,beta=1)
  return(lapply(1:length(A), function(x) alpha*A[[x]]+(beta*B[[x]])))

block_mul <- function (alpha,A)
  return(lapply(A, function(x) alpha*x))  

  
#' We can extract evaluation of our GP with a product of matrices $\Phi and \xi$, however, how to 
#' arrange our function the tensorisation will give us a natural way to arrange things

#' This function gives the matrice $\Phi$ as a list of matrices per blocks 

#' Take the vector subdivision and $x$ and get out the object 
#' $(phi_{{j,k}})(x_{j})_{j,k in}$
#' @return Phi_per_var: such that Phi_per_var[[j]][[k]][i] give a seq
#' of the hat function $\phi_{j,k,i}(X_i)$ where $(X_i)=(x^{(1)}_i,\cdots x^{(n)_i})$ and 
#' \phi_{j,k,i} is the hat function given by the subdivision[[j]][[k]]  
#' @return Phi_var_to_tensor: such that Phi_var_to_tensor[[j]] is the 
#' tensorisation of vectors Phi_per_var[[j]][[k]] for k in $\{1,...,dim_block[j]\}$ 
#' @return Big_Phi is the vector per block defined in the article


Phi_per_var <- function(subdivision,sub_partition,dim_block,X){
  Phi_per_var <- vector("list", length(subdivision))
  for (j in 1:length(dim_block)) {
    Phi_per_var[[j]] <- lapply(1:dim_block[j], 
                               function(k) basisCompute.lineqGP(X[,sub_partition[[j]][k]],
                                                                subdivision[[j]][[k]]))
  } 
  return(Phi_per_var)
}

#' This function take the phi per variables and get back the family $\phi_{\el_j}(x_{j})$
Phi_var_to_tensor <- function(Phi_per_var,subdivision,subdivision_size,dim_block){
  Phi_tensor <- list("vector", length(subdivision))
  n_row <- nrow(Phi_per_var[[1]][[1]])
  for(j in 1:length(subdivision)){
    Phi_tensor[[j]] <- matrix(0,nrow=n_row,prod(subdivision_size[[j]]))
    for (i in 1:n_row){
      Phi_tensor[[j]][i,] <- eval(parse(text = paste("Phi_per_var[[j]][[", 1:dim_block[j],
                                                     "]][", i, ", ]", sep = "", collapse = "%x%")))
    }
  }
  return(Phi_tensor)
}

Block_Phi<-function(X,sub_partition,subdivision,subdivision_size, dim_block){
  return (Phi_var_to_tensor(Phi_per_var(subdivision,sub_partition,dim_block,X),subdivision,subdivision_size, dim_block))
}



########################### Creation of the kernel matrice ################################
#' @param subdivision
#' @param Parameters list of parameters that indicate the length correlation of the kernels and the 
#' globale variance we consider 
#' @param type single parameter or a list of parameters that indicate the type of kernel we use, 
#' type[[j]][[k]] take value in: "gaussian", "matern32", "matern52", "exponential".


Gamma_var <- function(subdivision, Parameters, type, dim_block, nblock){
  Gam_var <- vector("list", nblock)
  names(Gam_var) <- paste("block", 1:nblock, sep = "")
  if(is.list(type)){
    for (j in 1:nblock) {
      Gam_var[[j]] <- lapply(1:dim_block[j],function(k) kernCompute(
        subdivision[[j]][[k]],subdivision[[j]][[k]],type[[j]][k],Parameters[[j]][c(1, k+1)]))
      names(Gam_var[[j]]) <- paste("x", partition[[j]], sep = "")
    }
  }
  else if(length(type)==length(subdivision)){
    for (j in 1:nblock) {
      Gam_var[[j]] <- lapply(1:dim_block[j],function(k) kernCompute(
        subdivision[[j]][[k]],subdivision[[j]][[k]],type[j],Parameters[[j]][c(1, k+1)]))
      names(Gam_var[[j]]) <- paste("x", partition[[j]], sep = "")
    }
  }
  else if (length(type)==1){
    for (j in 1:nblock) {
      Gam_var[[j]] <- lapply(1:dim_block[j],function(k) kernCompute(
        subdivision[[j]][[k]],subdivision[[j]][[k]],type,Parameters[[j]][c(1, k+1)]))
      names(Gam_var[[j]]) <- paste("x", partition[[j]], sep = "")
    }
  }
  else{
    stop("type is not valid")
  }
  names(Gam_var) <- paste("block", 1:length(subdivision), sep = "")
  return(Gam_var)
}

Gamma_var_to_tensor <- function(GammaVar, dim_block, J){
  Gamma_tensor <- list("vector",J)
  for(j in 1:J){
    Gamma_tensor[[j]] <- eval(parse(text = paste("GammaVar[[j]][[", 1:dim_block[j],
                                                 "]]", sep = "", collapse = "%x%")))
  }
  return(Gamma_tensor)
}

block_Gamma <- function(subdivision, Parameters, type, dim_block, J){
  return(Gamma_var_to_tensor(GammaVar(subdivision, Parameters, type), dim_block, J))
}

###############################################################################################
#' We will introduce here the functions that will be able to constructe and manipulate the kernel matrix
#' The kernel is obtained by tensorisation


#' We will create functions for different matrices we will neeed

#' Computation of the conditional covariance matrice using the inversion lemma formulation:
#' $(Z+UWV^{\top})^{-1}=Z^{-1}-Z^{-1}U(W^{-1}+ V^{t}Z^{-1}U)^{-1}V^{\top}Z^{-1}$ it is useful
#' when $U, V\in M_{n,m}, W \in M_{m,m} $ and $Z \in M_{n,n}$. the factorisation form is useful
#' when $m<<n$

# !!!!!! Carefull works only with W symetric matrix !!!!!!! inversion cholesky  !!!!!!!!! 

#As it is useful only if we know inv_Z we take it in parameters 
inv_lemma <- function(U,V,W,inv_Z){
  inv_W <- chol2inv(chol(W))
  return(inv_Z%-%inv_Z%*%U%*%chol2inv(chol(inv_W%*%t(V)%*%inv_Z%*%U))%*%t(V)%*%inv(Z))
}

#' The type in this function is to deal with a particular case where $Z=\tau^2 I_n$
#' To make computations faster

#' Block case, be careful about matrices we use
block_inv_lemma <- function(W,inv_Z,U,V=NULL){
  if(is.null(V)){
    inv_ZU <- inv_ZV <- block_prod(inv_Z,U)
  }
  else{
    inv_ZU <- block_prod(inv_Z,U)
    inv_ZV <- block_prod(inv_Z,V)
  }
    t_inv_ZV <- block_transpose(inv_ZV)
  return(block_to_vect(inv_Z)%-%block_to_vect(WU,type="line")%*%
         chol2inv(chol(block_to_vect(block_prod(t_inv_ZV,U)),type="sum"))%*%
         block_to_vect(t_inv_ZV,type="column"))
}


#'Fast computation when m>n

create_conditional_mean_1 <- function(K,Phi,Y,Z){
  Kt_Phi <- block_prod(K,transpose(Phi))
  return(block_to_vect(kt_Phi,type="column")%*%
         chol2inv(chol(block_to_vect(Z)%+%block_to_vect(prod(Phi,Kt_Phi)),type="sum"))%*%Y)
}

#'Fast computation when m<n and matrix is identity

create_conditional_mean_2 <- function(K,Phi,Y,tau){
  t_Phi<-block_transpose(Phi)
  Kt_Phi <- block_prod(K,t_Phi)
  inv_K <- block_inv(K)
  t_PhiPhi <- block_to_vect(block_prod(t_Phi,Phi),type="sum")
  return (block_to_vect(Kt_Phi,type="column")%*%
           (diag(nrow(Y))-chol2inv(chol(block_to_vect(inv_K)+block_scalar(phi,Kt_Phi)))%*%Y))
}


