#' Fit a supervised singular value decomposition (SupSVD) model
#'
#' This function fits the SupSVD model:
#' X=UV' + E, U=YB + F
#' where X is an observed primary data matrix (to be decomposed), U is a latent score
#' matrix, V is a loading matrix, E is measurement noise, Y is an observed
#' auxiliary supervision matrix, B is a coefficient matrix, and F is a
#' random effect matrix. \cr
#' It is a generalization of principal component analysis (PCA) or singular
#' value decomposition (SVD). It decomposes the primary data matrix X into low-rank
#' components, while taking into account potential supervision from any auxiliary
#' data Y measured on the same samples. \cr
#' \cr
#' See more details in 2016 JMVA paper "Supervised singular value decomposition
#' and its asymptotic properties" by Gen Li, Dan Yang, Andrew B Nobel and
#' Haipeng Shen.
#'
#' @param Y n*q (column centered) auxiliary data matrix, rows are samples and columns are stats::variables
#'  (must have linearly independent columns to avoid overfitting)
#' @param X n*p (column centered) primary data matrix, which we want to decompose.
#'  rows are samples (matched with Y) and columns are variables
#' @param r positive scalar, prespecified rank (r < min(n,p))
#'
#' @return list with components \item{B:}{q*r coefficient matrix of Y on the scores of X, maybe sparse if gamma=1}
#' \item{V:}{p*r loading matrix of X, with orthonormal columns}
#' \item{U:}{n*r score matrix of X, conditional expectation of random scores}
#' \item{se2:}{scalar, variance of measurement error in the primary data X}
#' \item{Sf:}{r*r diagonal covariance matrix, for random effects (see paper)}
#' Note: Essentially, U and V are the most important output for dimension
#' reduction purpose as in PCA or SVD.
#'
#'
#' @export
#'
#' @examples
#' r=2
#' Y <- matrix(rnorm(400,0,1),nrow=100)
#' B <- c(-1,1,-sqrt(3/2),-1)
#' B <- cbind(B,c(1,1,-1,sqrt(3/2)))
#' V <- matrix(rnorm(68*2),68,2)
#' Fmatrix <- matrix(MASS::mvrnorm(n=1*100,rep(0,2),matrix(c(9,0,0,4),2,2)),100,2)
#' E <- matrix(rnorm(100*68,0,3),100,68)
#' Yc <- scale(Y,center=TRUE,scale=FALSE)
#'
#'# Case 1 (supsvd) X = YBV^T+FV^T+E
#'X1 <- Y%*%tcrossprod(B,V)+tcrossprod(Fmatrix,V)+E
#'X1c <- scale(X1,center=TRUE,scale=FALSE)
#'SupPCA(Yc,X1c,r)
#'#  Case 2 (PCA) X = FV^T+E
#'X2 <- tcrossprod(Fmatrix,V)+E
#'X2c <-scale(X2,center=TRUE,scale=FALSE)
#'SupPCA(Yc,X2c,r)

#'# Case 3 (RRR) X = YBV^T+E
#'X3 <- Y%*%tcrossprod(B,V)+E
#'X3c <- scale(X3,center=TRUE,scale=FALSE)
#' SupPCA(Yc,X3c,r)
SupPCA <- function(Y, X, r){
  if(is.vector(Y)=="TRUE"){
    Y <- as.matrix(Y,ncol=1)
  }
  n <- nrow(X)
  p <- ncol(X)
  n1 <- nrow(Y)
  q <- ncol(Y)

  ## Pre-Check
  if (n != n1){
    stop("X does not match Y! exit...")
  } else if(r>min(n,p)){
    stop("Rank is too greedy! exit...")
  } else if(qr(Y)$rank != q){
    stop("Columns of Y are linearly dependent! exit...")
  } else if(max(abs(colMeans(X))) > 0.01 || max(abs(colMeans(Y))) > 0.01){
    stop("Columns of X and Y are not centered. exit...")
  }

  ## SVD start
  if(r == 1){
    S <- RSpectra::svds(X,r)
    U <- S$u
    D <- S$d
    D <- diag(as.matrix(D))
    V <- S$v
    U <- U%*%D
    E <- X - tcrossprod(U,V)
    se2 <- stats::var(as.vector(E))
    B <- tcrossprod(solve(crossprod(Y,Y)), Y) %*% U
    Sf <- as.matrix(diag(as.matrix(diag(as.matrix((1/n)*crossprod(U-Y%*%B, U-Y%*%B))))))

    temp1 <- X%*%V-Y%*%B  ## n*r
    temp2 <- X-tcrossprod(Y%*%B, V)  ##n*p
    trace1 <- psych::tr(crossprod(temp1,temp1))
    trace2 <- psych::tr(crossprod(temp2,temp2))
    logl <- (-n/2)*(log(det(as.matrix(Sf+se2*diag(r))))+(p-r)*log(se2))-(0.5/se2)*trace2 - 0.5*psych::tr(crossprod(temp1,temp1)%*%solve(Sf+se2*diag(r))) + (0.5/se2)*trace1
    rec <- c(logl)

    max_niter = 10^5
    convg_thres = 10^(-6)
    Ldiff = 1
    Pdiff= 1
    niter = 0

    while ((niter <= max_niter) && (Pdiff > convg_thres)){
      ## record last iter
      logl_old <- logl
      se2_old <- se2
      Sf_old <- Sf
      V_old <- V
      B_old <- B

      ## E step
      ## some critical values
      Sfinv <- solve(Sf)
      weight <- solve(diag(as.matrix(r))+se2*Sfinv)    ## r*r
      cond_Mean <- (se2*Y%*%B%*%Sfinv + X%*%V) %*% weight ## E(U|X), n*r
      cond_var <- Sf %*% (diag(as.matrix(r))-weight)
      cond_quad <- n*cond_var + crossprod(cond_Mean,cond_Mean)

      ## M step
      V <- crossprod(X, cond_Mean)%*%solve(cond_quad)
      se2 <- (psych::tr(X%*%(t(X)-2*tcrossprod(V, cond_Mean))) + n*psych::tr(crossprod(V,V)%*%cond_var) +
                psych::tr(tcrossprod(tcrossprod(cond_Mean,V)%*%V, cond_Mean)))/(n*p)
      B <- solve(crossprod(Y,Y)) %*% t(Y) %*% cond_Mean
      Sf <- (cond_quad + crossprod(Y%*%B,Y%*%B)- crossprod(Y%*%B,cond_Mean) -
               crossprod(cond_Mean, Y%*%B))/n

      ## S step
      newS <- RSpectra::svds(V%*%Sf%*%t(V), r)
      newV <- newS$u
      newSf <- newS$d
      newSf <- as.matrix(diag(as.matrix(newSf)))
      Sf <- newSf[1:r,1:r,drop = FALSE]
      B <- tcrossprod(B,V) %*% newV[,1:r,drop = FALSE]
      V <- newV[,1:r,drop = FALSE]

      ## log likelihood
      temp1 <- X%*%V - Y%*%B  ## n*r
      temp2 <- X-Y%*%B%*%t(V)  ## n*p
      logl <- (-n/2)*(log(det(as.matrix(Sf+se2*diag(r)))) + (p-r)*log(se2)) - (0.5/se2)*trace2 -
        0.5*psych::tr(crossprod(temp1,temp1)%*%solve(Sf+se2*diag(r))) + (0.5/se2)*trace1
      rec1 <- c(rec,logl)

      ## iteration termination
      Ldiff <- logl-logl_old  ## should be positive
      Pdiff <- norm(V-V_old, type = "F")^2
      niter <- niter+1
    }

    if (niter<max_niter){
      cat("EMS converges at precision",convg_thres, "after", niter, "iterations.")
    } else{
      cat("EMS NOT converges at precision",convg_thres, "after", niter, "iterations!!!")
    }

    ## re-order V, and correspondingly B and Sf, U (one simple remedy for the improper order of V)
    I <- order(timeSeries::colStdevs(X%*%V), decreasing = TRUE)
    V <- V[,I,drop = FALSE]
    B <- B[,I,drop = FALSE]
    Sf <- Sf[I,I, drop = FALSE]

    ## correct sign of V for identifiability
    ## also correct B and Sf
    signchange <- as.matrix(sign(V[1,,drop = FALSE]))
    V <- pracma::bsxfun("*", V, pracma::repmat(signchange,nrow(V),1))
    B <- pracma::bsxfun("*", B, pracma::repmat(signchange,nrow(B),1))
    Sf <- diag(signchange)%*%Sf%*%diag(signchange)

    ## output U
    Sfinv <- solve(Sf)
    weight <- solve(diag(as.matrix(r))+se2*Sfinv)  ## r*r
    U <- (se2*Y%*%B%*%Sfinv + X%*%V)%*%weight}
  else{
    S <- RSpectra::svds(X,r)
    U <- S$u
    D <- S$d
    D <- diag(D)
    V <- S$v
    U <- U%*%D
    E <- X - tcrossprod(U,V)
    se2 <- stats::var(as.vector(E))
    B <- tcrossprod(solve(crossprod(Y,Y)), Y) %*% U
    Sf <- diag(diag((1/n)*crossprod(U-Y%*%B, U-Y%*%B)))

    temp1 <- X%*%V-Y%*%B  ## n*r
    temp2 <- X-tcrossprod(Y%*%B, V)  ##n*p
    trace1 <- psych::tr(crossprod(temp1,temp1))
    trace2 <- psych::tr(crossprod(temp2,temp2))
    logl <- (-n/2)*(log(det(Sf+se2*diag(r)))+(p-r)*log(se2))-(0.5/se2)*trace2 -
      0.5*psych::tr(crossprod(temp1,temp1)%*%solve(Sf+se2*diag(r))) + (0.5/se2)*trace1
    rec <- c(logl)

    max_niter = 10^5
    convg_thres = 10^(-6)
    Ldiff = 1
    Pdiff= 1
    niter = 0

    while ((niter <= max_niter) && (Pdiff > convg_thres)){
      ## record last iter
      logl_old <- logl
      se2_old <- se2
      Sf_old <- Sf
      V_old <- V
      B_old <- B

      ## E step
      ## some critical values
      Sfinv <- solve(Sf)
      weight <- solve(diag(r)+se2*Sfinv)    ## r*r
      cond_Mean <- (se2*Y%*%B%*%Sfinv + X%*%V) %*% weight ## E(U|X), n*r
      cond_var <- Sf %*% (diag(r)-weight)
      cond_quad <- n*cond_var + crossprod(cond_Mean,cond_Mean)

      ## M step
      V <- crossprod(X, cond_Mean)%*%solve(cond_quad)
      se2 <- (psych::tr(X%*%(t(X)-2*tcrossprod(V, cond_Mean))) + n*psych::tr(crossprod(V,V)%*%cond_var) +
                psych::tr(tcrossprod(tcrossprod(cond_Mean,V)%*%V, cond_Mean)))/(n*p)
      B <- solve(crossprod(Y,Y)) %*% t(Y) %*% cond_Mean
      Sf <- (cond_quad + crossprod(Y%*%B,Y%*%B)- crossprod(Y%*%B,cond_Mean) -
               crossprod(cond_Mean, Y%*%B))/n

      ## S step
      newS <- RSpectra::svds(V%*%Sf%*%t(V), r)
      newV <- newS$u
      newSf <- newS$d
      newSf <- diag(newSf)
      Sf <- newSf[1:r,1:r,drop = FALSE]
      B <- tcrossprod(B,V) %*% newV[,1:r,drop = FALSE]
      V <- newV[,1:r,drop = FALSE]

      ## log likelihood
      temp1 <- X%*%V - Y%*%B  ## n*r
      temp2 <- X-Y%*%B%*%t(V)  ## n*p
      logl <- (-n/2)*(log(det(Sf+se2*diag(r))) + (p-r)*log(se2)) - (0.5/se2)*trace2 -
        0.5*psych::tr(crossprod(temp1,temp1)%*%solve(Sf+se2*diag(r))) + (0.5/se2)*trace1
      rec1 <- c(rec,logl)

      ## iteration termination
      Ldiff <- logl-logl_old  ## should be positive
      Pdiff <- norm(V-V_old, type = "F")^2
      niter <- niter+1
    }

    if (niter<max_niter){
      cat("EMS converges at precision",convg_thres, "after", niter, "iterations.")
    } else{
      cat("EMS NOT converges at precision",convg_thres, "after", niter, "iterations!!!")
    }

    ## re-order V, and correspondingly B and Sf, U (one simple remedy for the improper order of V)
    I <- order(timeSeries::colStdevs(X%*%V), decreasing = TRUE)
    V <- V[,I,drop=FALSE]
    B <- B[,I,drop=FALSE]
    Sf <- Sf[I,I,drop=FALSE]

    ## correct sign of V for identifiability
    ## also correct B and Sf
    signchange <- sign(V[1,])
    V <- pracma::bsxfun("*", V, pracma::repmat(signchange,nrow(V),1))
    B <- pracma::bsxfun("*", B, pracma::repmat(signchange,nrow(B),1))
    Sf <- diag(signchange) %*% Sf %*%diag(signchange)

    ## output U
    Sfinv <- solve(Sf)
    weight <- solve(diag(r)+se2*Sfinv)  ## r*r
    U <- (se2*Y%*%B%*%Sfinv + X%*%V)%*%weight
  }
  return(list(U=as.matrix(U), V=as.matrix(V), B=as.matrix(B), Sf=as.matrix(Sf), se2=se2))
}

