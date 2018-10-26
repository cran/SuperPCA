library(matlab)
#' Compute tensor product over multiple dimensions using Khatri-Rao product
#'
#' @param L List of length K, with matrix entries L{k}: m_k X R;
#' @param range Column indices over which to apply tensor product.Default is all columns (range = [1:R])
#'
#' @return Ans Array of size m_1 X m_2 ... X m_K, where the [i1,...iK]; entry is the sum of product L{k}(i1,r)*...*L{K}(iK,r) over all r in range.
#' @import matlabr
#' @export TensProd
#'
#' @examples
#' L <- list(matrix(rnorm(10*10),10,10),matrix(rnorm(10*10),10,10),matrix(rnorm(1000*10),1000,10))
#' TensProd(L)
#'
TensProd <- function(L,range=NA){
  K <- length(L)
  m <- rep(1,K)
  for(i in 1:K){
    m[i] = nrow(L[[i]])
    R = ncol(L[[i]])
  }
  if(is.na(range)=="FALSE"){
    range <- seq(1:range)
    if(max(range)>R){
      stop("Range exceeds rank!")
    }
    # customize
    newL <- L
    for(i in 1:K){
      newL[[i]] <- L[[i]][,max(range)]
    }
    L <- newL

  }
  tempL <- L[K:2]
  xx <- length(tempL)
  if(length(range)==1 && is.na(range)=="FALSE"){
    for(i in 1:xx){
      tempL[[i]] <- matrix(tempL[[i]])
    }
  }
  matX <- as.matrix(L[[1]])%*%t(kr(tempL)) #X_(1)
  Ans <- array(matX,m)
  return(Ans)
}
