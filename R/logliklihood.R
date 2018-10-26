#matlab::ones()
#matlab::zeros
#pracma::blkdiag
#pracma::arrayfun
loglikelihood <- function(X,Y,B0,B,V_joint,V_ind,se2,Sf0,Sf){
  #This funciton calc sum of log likelihood of Y
  #The input is data and the output of SIPCA_LG_v
  #The output is proportional to the sum of log likelihood (without the constant term)
  #
  # The larger the better !!!!

  #get dimension
  n <- nrow(X)
  q <- ncol(X)
  K <- length(Y)
  p <- matlab::zeros(1,K)
  for(k in 1:K){
    n_temp <- nrow(Y[[k]])
    p[k] <- ncol(Y[[k]])
    if(n_temp!=n){
      stop("Sample mismatch!")
    }
  }

  # critical value
  grandY <- NULL
  for(k in 1:K){
    grandY <- cbind(grandY,Y[[k]])
  }
  grandse2 <- NULL #1*sum(p)
  grandSf0 <- t(diag(Sf0))#1*r0
  grandSf <- NULL #1*sum(r)
  grandV_temp <- NULL #combination of loadings: first few long columns are joint loadings, subsequent diagonal blocks are individual loadings
  grandB_temp <- NULL #concatenate B
  grandse2 <- rep(1,p[1])*se2[1]
  grandSf <- t(diag(Sf[[1]]))
  grandV_temp <- V_ind[[1]]
  grandB_temp <- B[[1]]
  for(k in 2:K){
    grandse2 <- c(grandse2, (rep(1,p[k])*se2[k]))
    grandSf <- cbind(grandSf, t(diag(Sf[[k]])))
    grandV_temp <- pracma::blkdiag(grandV_temp, V_ind[[k]])
    grandB_temp <- cbind(grandB_temp, B[[k]])
  }
  grandse2_inv <- 1/grandse2 #1*sum(p)
  grandSf_inv <- 1/grandSf #1*sum(r)
  grandSf0_inv <- 1/grandSf0 #1*r0
  grandV <- cbind(V_joint, grandV_temp)# sum(p)*(r0+sum(r))
  grandB <- cbind(B0, grandB_temp) #q*(r0+sum(r))
  Delta1 <- pracma::arrayfun("*", t(grandV), grandse2_inv)%*%grandV #[r0+sum(r)]*[r0+sum(r)]
  Delta2_inv <- solve(diag(cbind(grandSf0_inv,grandSf_inv))+Delta1)
  temp <- grandV%*%tcrossprod(Delta2_inv, grandV) #sum(p)*sum(p)
  SigmaY_inv <- diag(grandse2_inv)-pracma::arrayfun("*", pracma::arrayfun("*", temp, grandse2_inv), t(grandse2_inv)) #sum(p)*sum(p), not diagonal because of common structure, diff from SupSVD
  temp <- pracma::arrayfun("*", grandV, sqrt(cbind(grandSf0,grandSf)))
  SigmaY <- diag(as.vector(grandse2))+tcrossprod(temp, temp)

  #loglikelihood terms
  term1 <- -n*sum(log(diag(chol(SigmaY)))) #log det
  temp <- grandY - X%*%tcrossprod(grandB, grandV)
  term2 <- -1/2*sum(diag((SigmaY_inv%*%crossprod(temp, temp))))

  #final output
  out <- term1+term2

  return(out)
}
