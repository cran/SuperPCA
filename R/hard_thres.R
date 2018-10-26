hard_thres <- function(a,lam){
  ind <- (abs(a)<lam)
  a[ind] <- 0
  out <- a
  return(out)
}