#glmnet::glmnet
#pracma::bsxfun
EstU_fun <- function(X,U,relation = 'linear',sparsity = 0){
  q <- ncol(X)
  n <- nrow(U)
  r <- ncol(U)
  fX <- matrix(0,n,r)
  if(relation == 'linear' && sparsity == 0){# linear non-sparse
    B <- solve(crossprod(X,X), t(X))%*%U # q*sum(r)
    fX <- X%*%B # naturally column centered b/c X are centered
  }
  else if(relation == 'linear' && sparsity == 1){# linear and sparse
    B <- matrix(0,q,r)
    for(i in 1:r){
      fit <- glmnet::glmnet(X,U[,i],lambda.min.ratio = 0, standardize = TRUE)
      BIC_score <- n*stats::deviance(fit)+log(n)*fit$df
      ind <- which(BIC_score == min(BIC_score))
      B[,i] <- fit$beta[,ind]
    }
    fX <- X%*%B # naturally column centered b/c X are centered
  }
  else if(relation == 'univ_kernel'){
    if(q != 1){
      stop("Cannot deal with multivariate covariates...")
    }
    for(i in 1:r){# for each column of U
      out <- stats::ksmooth(X,U[,i])
      fX[,i] <- out$y
    }
    # center each column of fX (b/c we assume X and Y_i are all column centered)
    fX <- pracma::bsxfun("-", fX, colMeans(fX))
  }
  else{
    stop("No such relation function available...")
  }
}
