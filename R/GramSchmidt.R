#pracma::gramSchmidt()
#library(pracma)
GramSchmidt <- function(A1, B){
  #print(match.call())
  nargin <- length(as.list(match.call())) -1
  #print(nargin)

  m <- nrow(A1)
  n <- ncol(A1)
  if(nargin == 2){
    m1 <- nrow(B)
    n2 <- ncol(B)
    Q <- matrix(0,m1,n+n2)
    R <- matrix(0,n+n2,n+n2)
    if(m != m1){
      return(list(Q1=Q1,Q=Q,R=R))
    }
    A <- cbind(B, A1)
    for(j in 1:(n+n2)){
      v <- A[,j,drop=F]
      for(i in seq(1,(j-1))){
        if (j==1)break;
        R[i,j] <- crossprod(Q[,i], A[,j])
        v <- v - R[i,j]*Q[,i,drop = FALSE]
      }
      R[j,j] <- norm(as.matrix(v),type="f")
      Q[,j] <- as.matrix(v)/R[j,j]
    }
    Q1 <- Q[,(n2+1):(n+n2)]
  }
  else{
    G <- pracma::gramSchmidt(A1)
    Q <- G$Q
    R <- G$R
    Q1 <- Q
  }
  return(list(Q1=Q1,Q=Q,R=R))
}
