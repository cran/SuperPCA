

#' Compute a string of Khatri-Rao products
#'
#' The function compute a string of Khatri-Rao products. A*B*C*..., "*" denotes the Khatri-Rao product. If the list only contains two matrices, kr() return the Khatri-Rao product of two matrices.
#'
#' @param x a list of matrices
#'
#' @return result of a string of Khatri-Rao products
#' @import Matrix
#' @export kr
#'
#' @examples
#' #ex1
#' m1 <- matrix(1:9,3,3)
#' m2 <- matrix(1:12,4,3)
#' m3 <- matrix(13:27,5,3)
#' l1 <- list(m1,m2,m3)
#' kr(l1)
#' #ex2
#' l2 <- list(m1,m3)
#' kr(l2)
kr<-function(x){
  if(length(x)==1){
    out <- x[[1]]
  }else{
    out <- Matrix::KhatriRao(x[[1]],x[[2]])
    if(length(x) > 2){
      for(i in 3:length(x)){
        out <- Matrix::KhatriRao(out,x[[i]])
      }
    }
  }
  return(out)
}
