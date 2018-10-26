#Matrix::Matrix::chol()
#Matrix::Matrix::tMatrix::crossprod()
#Matrix::Matrix::crossprod()
#Matrix::Matrix::solve()
#psych::psych::tr()
#Mass::MASS::mvrnorm
#pracma::pracma::num2str()

# library(Matrix)
#' Using EM algorithm to fit the SupCP model
#'
#' @param Y n*q full column rank reponse matrix(necessarily n>=q)
#' @param X n*p1*...*pk design array
#' @param R fixed rankd of approximation, R<=min(n,p)
#' @param AnnealIters Annealing iterations (default = 100)
#' @param ParafacStart binary argument for wether to initialize with Parafac factorization (default = 0)
#' @param max_niter maximum number of iterations (default = 1000)
#' @param convg_thres convergence threshold for difference in log likelihood (default = 10^-3)
#' @param Sf_diag whether Sf is diagnal (default=1,diagnoal)

#'
#' @return list with components
#' \item{B:}{q*r coefficient matrix for U~Y}
#' \item{V}{list of length K-1. V[k] is a p*r coefficient matrix with columns of norm 1}
#' \item{U:}{Conditional expectation of U: n*r}
#'
#' \item{se2:}{scalar, var(E)}
#' \item{Sf:}{r*r diagonal matrix, cov(F)}
#' \item{rec:}{log likelihood for each iteration}
#' @export
#'
#' @examples
#' sigmaF <- diag(c(100,64,36,16,4))
#' # F matrix n*r
#' Fmatrix1 <- matrix(MASS::mvrnorm(n=500,rep(0,5),sigmaF),100,5)
#' U<-Fmatrix1
#' V1 <- matrix(stats::rnorm(10*5),10,5)
#' V2 <- matrix(stats::rnorm(10*5),10,5)
#' L <- list(U,V1,V2)
#' X <- TensProd(L)
#' Y <- matrix(stats::rnorm(100*10),100,10)
#' R <-3
#' SupParafacEM(Y,X,R)
#'
SupParafacEM <- function(Y,X,R,AnnealIters=100,ParafacStart=0,max_niter=1000,convg_thres=10^-3,Sf_diag=1){

  n1 <- nrow(Y)
  q <-ncol(Y)
  m <- dim(X)
  n <- m[1]
  L <- length(m) # number of modes
  K <- L-1
  p <- prod(m[2:L]) #p1*p2*...*pK
  V <- list()

  # Pre-check
  if(n!=n1){
    stop("X does not match Y! exit...")
  }else if(qr(Y)$rank!=q){
    stop("Columns of Y are linearly dependent! exit...")
  }

  Index <- 1:L
  IndexV <- 1:(L-1)

  ## initialize via parafac
  if(ParafacStart==1){
    Init <- Parafac(X,R)$U #still has randomness in initial value, but in another layer
    V <- Init[2:L]
  }else{
    for(l in 2:L){
      V[[l-1]] <- normc(matrix(stats::rnorm(m[l]*R),nrow=m[l],ncol=R,byrow=T))
    }
  }
  #V <- Vstart
  Vmat <- matrix(0,nrow=p,ncol=R)# a very long matrix(p can be very large)
  for(r in 1:R){
    Temp <- TensProd(V,r)
    Vmat[,r] <- array(Temp,c(prod(dim(Temp)),1))
  }
  Xmat <- array(aperm(X,c(2:L,1)),c(prod(dim(X))/n,n))
  U <- crossprod(Xmat,Vmat)
  E <- X-TensProd(lapply(rapply(list(U,V), enquote, how="unlist"), eval))
  se2 <- stats::var(c(E))
  B <- tcrossprod(solve(crossprod(Y,Y)),Y)%*%U
  if(Sf_diag==1){
    Sf <- diag(diag((1/n)*crossprod((U-Y%*%B),(U-Y%*%B)))) # R*R, diagonal
  }else{
    Sf <- (1/n)*crossprod(U-Y%*%B,U-Y%*%B)
  }

  ##Compute determinant exactly, using Sylvester's  determinant theorem
  ##https://en.wikipedia.org/wiki/Determinant#Properties_of_the_determinant
  #MatForDet = sqrt(Sf)*Vmat'*Vmat*sqrt(Sf)./se2+eye(R); %R X R

  MatForDet <- (Sf^.5)%*%crossprod(Vmat,Vmat)%*%(Sf^.5)/se2+diag(R) # R*R
  #uses woodbury identity and trace properties
  logdet_VarMat <- 2*sum(log(diag(chol(MatForDet))))+p*log(se2)

  ResidMat <- t(Xmat)-Y%*%tcrossprod(B,Vmat) #n*p

  if(Sf_diag==1){
    Sfinv <- diag(1/diag(Sf))
  }else{
    Sfinv <- solve(Sf)
  }

  Trace <- (1/se2)*sum(diag(tcrossprod(ResidMat,ResidMat)))-(1/se2^2)*sum(diag(t(Vmat)%*%crossprod(ResidMat,ResidMat)%*%Vmat%*%solve(Sfinv+(1/se2)*crossprod(Vmat,Vmat))))
  logl <- (-n/2)*(logdet_VarMat)-0.5*Trace
  rec <- logl

  niter <- 1
  Pdiff <- convg_thres+1
  while(niter<=max_niter && (abs(Pdiff)>convg_thres)){
    cat(sprintf('This is the %.d th iterations, the Pdiff= %.4g \n',niter,Pdiff))
    niter <- niter+1

    #recod last iter
    logl_old <- logl
    se2_old <- se2
    Sf_old <- Sf
    Vmat_old <- Vmat
    V_old <- V
    B_old <- B

    #E step

    if(Sf_diag==1){
      Sfinv <- diag(1/diag(Sf))
    }else{
      Sfinv <- solve(Sf)
    }

    weight <- solve(crossprod(Vmat,Vmat)+se2*Sfinv) #r*r
    cond_Mean <- as.matrix((se2*Y%*%B%*%Sfinv+crossprod(Xmat,Vmat))%*%weight) #E(U|X), n*r
    U <- cond_Mean
    cond_Var <- solve((1/se2)*crossprod(Vmat,Vmat)+Sfinv) # cov(U(i)|X) r*r
    # %Add noise to the conditional mean of U.
    # %Variance of noise is a decreasing percantage of the variance of the
    # %true conditional mean.

    #cat("solve:",pracma::num2str(det((1/se2)*crossprod(Vmat,Vmat)+Sfinv)))
    #cat("matrix",print((1/se2)*crossprod(Vmat,Vmat)+Sfinv))

    if(niter<AnnealIters){
      anneal <- (AnnealIters-niter)/AnnealIters
      U <- matrix(MASS::mvrnorm(prod(dim(cond_Mean)),colMeans(cond_Mean),anneal*diag(matrixStats::colVars(U))),dim(cond_Mean)[1],dim(cond_Mean)[2])
    }
    cond_Mean <- U
    cond_quad <- n*cond_Var+crossprod(U,U) #E(U'U|X),r*r

    #Estimate V's
    for(l in 2:L){
      newIndex <- c(Index[-match(l,Index)],match(l,Index))
      ResponseMat <- array(aperm(X,newIndex),c(prod(dim(X))/m[l],m[l]))
      PredMat <- matrix(0,prod(m[Index[-match(l,Index)]]),R)
      VParams <- matrix(0,prod(m[Index[-match(c(1,l),Index)]]),R)
      for(r in 1:R){
        Temp <- TensProd(lapply(rapply(list(U,V[IndexV[-(l-1)]]), enquote, how="unlist"), eval),r)
        PredMat[,r] <- array(Temp,c(prod(dim(Temp)),1))
        if(L==3){
          Temp <- V[[IndexV[-(l-1)]]][,r,drop=FALSE]
        }else{
          Temp <- TensProd(V[IndexV[-(l-1)]],r)
        }
        VParams[,r] <- array(Temp,c(prod(dim(Temp)),1))
      }
      V[[l-1]] <- crossprod(ResponseMat,PredMat)%*%solve(crossprod(VParams,VParams)*cond_quad)
      #print(det(crossprod(VParams,VParams)*cond_quad))
      # V{l-1}=normc(V{l-1})
    }

    #estimate B
    B <- solve(crossprod(Y,Y),t(Y))%*%U

    #estimate Sf
    for(r in 1:R){
      Temp <- TensProd(V,r)
      Vmat[,r] <- array(Temp,c(prod(dim(Temp)),1))
    }
    se2 <- (sum(diag(crossprod(Xmat,(Xmat-2*tcrossprod(Vmat,cond_Mean)))))+n*sum(diag(crossprod(Vmat,Vmat)%*%cond_Var))+sum(diag(tcrossprod(cond_Mean,Vmat)%*%Vmat%*%t(cond_Mean))))/(n*p)

    #estimate diagnoal entries of covariance
    if(Sf_diag==1){
      Sf <- diag(diag((cond_quad+crossprod(Y%*%B,Y%*%B)-crossprod(Y%*%B,cond_Mean)-crossprod(cond_Mean,(Y%*%B)))/n))
    }else{# estimate full covariance
      Sf <- (cond_quad+crossprod(Y%*%B,Y%*%B)-crossprod(Y%*%B,cond_Mean)-crossprod(cond_Mean,Y%*%B))/n
    }

    #scaling
    for(l in 2:L){
      V[[l-1]] <- normc(V[[l-1]])
    }
    VmatS <- Vmat
    for(r in 1:R){
      Temp <- TensProd(V,r)
      VmatS[,r] <- array(Temp,c(prod(dim(Temp)),1))
    }
    Bscaling <- matrix(rep(1,q),q,1)%*%sqrt(colSums(Vmat^2))
    B <- B*Bscaling
    Sfscaling <- matrix(as.matrix(sqrt(colSums(Vmat^2)))%o%as.matrix(sqrt(colSums(Vmat^2))),length(sqrt(colSums(Vmat^2))))
    Sf <- Sf*Sfscaling
    Vmat <- VmatS

    # Calc likelihood
    if(Sf_diag==1){
      Sfinv <- diag(1/diag(Sf))
    }else{
      Sfinv <- solve(Sf)
    }

    ResidMat <- t(Xmat)-Y%*%tcrossprod(B,Vmat) #n*p
    MatForDet <- (Sf^.5)*crossprod(Vmat,Vmat)%*%(Sf^.5)/se2+diag(R) #R*R
    logdet_VarMat <- 2*sum(log(diag(chol(MatForDet))))+p*log(se2)
    Trace <- (1/se2)*sum(diag(tcrossprod(ResidMat,ResidMat)))-(1/se2^2)*sum(diag(t(Vmat)%*%t(ResidMat)%*%ResidMat%*%Vmat%*%solve(Sfinv+(1/se2)*crossprod(Vmat,Vmat))))
    logl <- (-n/2)*(logdet_VarMat)-0.5*Trace
    rec <- cbind(rec,logl)

    #iteration termination
    Ldiff <- logl-logl_old #should be positive
    Pdiff <- Ldiff
  }

  if(niter<max_niter){
    cat("EM converges at precision",pracma::num2str(convg_thres,fmt=10),"after",pracma::num2str(niter),"iterations.")
  }else{
    cat("EM does not converge at precision", pracma::num2str(convg_thres,fmt=10),"after",pracma::num2str(max_niter),"iterations!!!")
  }

  #re-order parameters
  tmp <- sort(diag(Sf),index.return=TRUE,decreasing = TRUE)
  I <- tmp$ix
  for(k in 1:(L-1)){
    V[[k]] <- V[[k]][,I]
  }
  B <- B[,I]
  Sf <- Sf[I,I]
  U <- U[,I]
  return(list(B=B,V=V,U=U,se2=se2,Sf=Sf,rec=rec))
}
