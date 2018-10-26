#library(matlab)
#library(pracma)

#' Performs parafac factorization via ALS
#'
#' @param Y Array of dimension m1*m2*...*mK
#' @param R Desired rank of the factorization. Defalust is all columns(range=1:R)
#'
#' @return list with components
#' \item{U:}{List of basis vectors for parafac factorization, U[k]: mk*R for k=1,...,K. Columns of U[K] have norm 1 for k>1. }
#' \item{SqError:}{Vector of squared error for the approximation at each iteration (should be non-increasing)}
#' @export Parafac
#'
#' @examples
#' A <- array(stats::rnorm(100*10*10), dim=c(100,10,10))
#' Parafac(A,4)
Parafac <- function(Y,R){
  m <- dim(Y)
  L <- length(m)
  U <- list()
  for(l in 2:L){
    U[[l]] <- normc(matrix(stats::rnorm(m[l]*R),m[l],R))
  }
  Index <- 1:L
  SqError <- NULL
  i <- 0
  thresh <- 10^(-1)
  SqErrorDiff <- thresh+1
  Yest <- array(0,m)
  iters <- 1000

  while(i<iters && SqErrorDiff>thresh){

    i<-i+1
    Yest_old <- Yest
    for(l in 1:L){
      newIndex <- c(Index[-match(l,Index)],match(l,Index))
      permuteY <- aperm(Y,newIndex)
      ResponseMat <- array(permuteY,c(prod(dim(Y))/m[l],m[l]))
      PredMat <- matrix(0,prod(m[Index[-match(l,Index)]]),R)

      for(r in 1:R){
        Temp <- TensProd(U[Index[-match(l,Index)]],r)
        PredMat[,r] = array(Temp,c(prod(dim(Temp)),1))
      }
      U[[l]] <- t(ResponseMat)%*%PredMat%*%solve(crossprod(PredMat,PredMat))
      if(l>1){
        U[[l]] <- normc(U[[l]])
      }
    }
    Yest <- TensProd(U)

    Error <- Y-Yest
    SqError[i] <- sum(matrix(array(Error)^2))
    ErrorDiff <- Yest-Yest_old
    SqErrorDiff <- sum(matrix(array(ErrorDiff)^2))
  }
  #print(count)
  return(list(U=U,SqError=SqError))
}
