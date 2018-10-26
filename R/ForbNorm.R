#pracma::bsxfun
#pracma::repmat
ForbNorm <- function(U_est, U){
  ind <- t(sign(diag(crossprod(U_est,U))))
  U_est <- pracma::bsxfun("*", U_est, pracma::repmat(ind,nrow(U_est),1))

  loss <- norm(U_est - U, type = "F")
  return (loss)
}
