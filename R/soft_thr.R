soft_thr<-function(In, lambda){
  out <- sign(In)*pmax(abs(In) - lambda, 0)
  return(out)
}
