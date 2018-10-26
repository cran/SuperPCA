PrinAngle <- function(V1, V2, ind=1){
  # if ind=1(default), calc max principal angle
  # if ind=2, calculates all principal angles between column space of V1 and V2

  p1 <- nrow(V1)
  r1 <- ncol(V1)
  p2 <- nrow(V2)
  r2 <- ncol(V2)
  if(p1 != p2){
    stop("Input must be matched")
  }
  S1 <- svd(V1)
  V1 <- S1$u
  S2 <- svd(V2)
  V2 <- S2$u
  if(ind == 1){
    angle <- 180/pi*acos(min(svd(crossprod(V1,V2))$d))
  }
  else if(ind == 2){
    angle <- 180/pi*acos(svd(crossprod(V1,V2))$d)
  }
  return(angle)
}
