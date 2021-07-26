#' Supervised Integrated Factor Analysis
#'
#' @param X n*q matrix, centered covariate data
#' @param Y 1*K cell array, each cell is a n*pi centered primary
#'                   data, should roughly contribute equally to joint struct
#' @param r0 scalar, prespecified rank of common structure
#' @param r 1*K vector, prespecified rank of specific structures
#' @param relation 'linear' (default), use linear function to model
#'                           the relation between U and covariates
#'                        'univ_kernel', use kernel methods for single covariates
#' @param sparsity 1 (default), when est B0 or B, use LASSO with BIC
#'                           to select the best tuning, suitable for high
#' @param max_niter default 1000, max number of iteration
#' @param convg_thres default 0.01, overall convergence threshold
#' @param type "A" or "Anp" or "B" or "Bnp", condition to use
#'
#' @return list with components
#' \item{B0:}{q*r0 matrix, coefficient for joint structure (may be sparse)}
#' \item{B:}{1*K cell array, each is a q*ri coefficient matrix (may be sparse)}
#' \item{V_joint:}{sum(p)*r0 matrix, stacked joint loadings, with orthonormal
#'                   columns.}
#'
#' \item{V_ind:}{1*K array, each is a pi*ri loading matrix, with
#'                   orthonormal columns}
#' \item{se2:}{1*K vector, noise variance for each phenotypic data set}
#' \item{Sf0:}{r0*r0 matrix, diagonal covariance matrix}
#' \item{Sf:}{1*K array, each is a ri*ri diagonal covariance matrix}
#' \item{EU:}{n*(r0+sum(ri)) matrix, conditional expectation of joint and individual scores}
#' @export
#'
#' @examples
#' \dontrun{
#' r0 <- 2
#' r <- c(3,3)
#' V <- matrix(stats::rnorm(10*2),10,2)
#' Fmatrix <- matrix(MASS::mvrnorm(n=1*500,rep(0,2),matrix(c(9,0,0,4),2,2)),500,2)
#' E <- matrix(stats::rnorm(500*10,0,3),500,10)
#' X <- tcrossprod(Fmatrix,V)+E
#' X <-scale(X,center=TRUE,scale=FALSE)
#' Y1 <- matrix(stats::rnorm(500*200,0,1),500,200)
#' Y2 <- matrix(stats::rnorm(500*200,0,1),500,200)
#' Y <- list(Y1,Y2)
#' SIFA(X,Y,r0,r,max_niter=200)
#' }
SIFA<- function(X,Y,r0,r,relation = 'linear', sparsity = 0,max_niter=1000, convg_thres=0.01, type ="Anp"){


  #Check dimension
  K <- length(Y)
  p <- rep(0,K)
  n <- nrow(X)
  q <- ncol(X)
  if(n <= q && sparsity==0 ){
    warning("Must use sparse estimation to avoid overfitting!!!")
    sparsity <- 1
  }
  for(k in 1:K){
    n_temp <- nrow(Y[[k]])
    p[k] <- ncol(Y[[k]])
    if(n_temp != n){
      stop("Sample mismatch!")
    }
  }
  if(sum((r0+r)>min(n,p))){
    stop("Too greedy on ranks!")
  }

  # initial value
  grandY <- NULL
  for(k in 1:K){
    grandY <- cbind(grandY,Y[[k]])
  }

  grandU<-NULL
  grandV<-NULL

  #### A_np
  if(type =="Anp"){
    S <- RSpectra::svds(grandY,r0)
    U_joint_ini <- S$u
    D_joint_ini <- diag(S$d,nrow=length(S$d))
    V_joint <- S$v
    U_joint <- U_joint_ini%*%D_joint_ini
    V_ind <- list()
    Sf <- list()

    loc3 <- 1
    loc4 <- sum(p[1:1])
    Vcurrent <- V_joint[loc3:loc4,]
    Ycurrent <- Y[[1]]-tcrossprod(U_joint,Vcurrent)
    Sk <- RSpectra::svds(Ycurrent,r[1])
    U_k_ini <- Sk$u
    D_k_ini <- diag(Sk$d,nrow=length(Sk$d))
    V_k <- Sk$v
    U_k <- U_k_ini%*%D_k_ini
    V_ind[[1]] <- V_k
    grandU <- cbind(grandU, U_k)
    grandV <- pracma::blkdiag(V_ind[[1]])

    for(k in 2:K){
      loc3 <- sum(p[1:(k-1)])+1
      loc4 <- sum(p[1:k])
      Vcurrent <- V_joint[loc3:loc4,]
      Ycurrent <- Y[[k]]-tcrossprod(U_joint,Vcurrent)
      Sk <- RSpectra::svds(Ycurrent,r[k])
      U_k_ini <- Sk$u
      D_k_ini <- diag(Sk$d,nrow=length(Sk$d))
      V_k <- Sk$v
      U_k <- U_k_ini%*%D_k_ini
      V_ind[[k]] <- V_k
      grandU <- cbind(grandU, U_k)
      grandV <- pracma::blkdiag(grandV,V_ind[[k]])
    }

    grandU <- cbind(U_joint, grandU)
    grandV <- cbind(V_joint, grandV)


    fX <- EstU_fun(X,grandU,relation,sparsity)

    se2 <- rep(0,K)
    loc1 <- r0+1
    loc2 <- loc1 + r[1]-1
    loc3 <- 1
    loc4 <- sum(p[1:1])
    U_k <- grandU[,loc1:loc2]
    fkX <-fX[,loc1:loc2]
    Sf[[1]] <- diag(matrixStats::colSds(U_k-fkX)^2)
    temp <- tcrossprod(grandU,grandV)
    se2[1] <- norm(Y[[1]]-temp[,loc3:loc4],"f")^2/(n*p[1])


    for(k in 2:K){
      loc1 <- r0+sum(r[1:(k-1)])+1
      loc2 <- loc1 + r[k]-1
      loc3 <- sum(p[1:(k-1)])+1
      loc4 <- sum(p[1:k])
      U_k <- grandU[,loc1:loc2]
      fkX <-fX[,loc1:loc2]
      Sf[[k]] <- diag(matrixStats::colSds(U_k-fkX)^2)
      temp <- tcrossprod(grandU,grandV)
      se2[k] <- norm(Y[[k]]-temp[,loc3:loc4],"f")^2/(n*p[k])
    }
    Sf0 <- diag(matrixStats::colSds(U_joint-fX[,1:r0])^2,nrow=length(matrixStats::colSds(U_joint-fX[,1:r0])^2))# a matrix

    # disp('Initial est done!')
    loglik <- loglikelihood1(Y,fX,V_joint,V_ind,se2,Sf0,Sf)
    recloglik <- loglik
    maxloglik <- loglik

    niter <- 0
    diff <- Inf
    # diff_max <- inf
    # tol_thres <- -10
    recdiff <- list()
    while(niter <= max_niter && abs(diff)>convg_thres){
      cat(sprintf('This is the %.d th iterations, the diff= %.4g \n',niter,diff))
      niter <- niter+1
      # record last iter
      grandV_old <- grandV

      # E stop
      # some critical values
      grandse2 <- NULL # 1*sum(p)
      grandSf0 <- t(diag(Sf0)) # 1*r0
      grandSf <- NULL # 1*sum(r)
      Delta_temp <- NULL # 1*sum(r)
      for(k in 1:K){
        grandse2 <- c(grandse2, rep(1,p[k])*(se2[k]))
        grandSf <- c(grandSf, t(diag(Sf[[k]])))
        Delta_temp <- c(Delta_temp, rep(1,r[k])*(1/se2[k]))
      }

      grandse2_inv <- 1/grandse2 #1*sum(p)
      grandSf_inv <- 1/grandSf #1*sum(r)
      grandSf0_inv <- 1/grandSf0 #1*r0
      Delta1 <- pracma::arrayfun("*", t(grandV), grandse2_inv)%*%grandV #[r0+sum(r)]*[r0+sum(r)]
      Delta2_inv <- solve(diag(c(grandSf0_inv, grandSf_inv))+Delta1)
      temp <- grandV%*%tcrossprod(Delta2_inv,grandV)
      SigmaY_inv <- diag(grandse2_inv)-pracma::arrayfun("*",pracma::arrayfun("*",temp,grandse2_inv),t(grandse2_inv))
      VSigmaYinvV <- Delta1-Delta1%*%Delta2_inv%*%Delta1 #sum(p)*sum(p)
      EU <- fX%*%(diag(sum(r)+r0)-pracma::arrayfun("*",VSigmaYinvV,c(grandSf0,grandSf)))+grandY%*%SigmaY_inv%*%pracma::arrayfun("*",grandV,c(grandSf0,grandSf))
      covU <- diag(c(grandSf0,grandSf))-pracma::arrayfun("*",pracma::arrayfun("*",VSigmaYinvV,c(grandSf0,grandSf)),t(c(grandSf0,grandSf)))

      # M step
      # est V
      # for iter 1 = 1:3 alternate between V_joint and V_ind
      loc1 <- r0+1
      loc2 <- loc1 + r[1]-1
      loc3 <- 1
      loc4 <- sum(p[1:1])
      EU0 <- EU[,(1:r0)]
      EUk <- EU[,(loc1:loc2)]
      EUkU0 <- n*covU[(loc1:loc2),(1:r0)]+crossprod(EUk,EU0) # r(k)*r0
      EU0U0 <- n*covU[(1:r0),(1:r0)]+crossprod(EU0, EU0) #r0*r0
      Vk <- V_ind[[1]] #p(k)*r(k)

      # V_0k
      V_0k <- (crossprod(Y[[1]],EU0)-Vk%*%EUkU0)%*%solve(EU0U0)
      V_joint[(loc3:loc4),] <- V_0k

      # V_k
      S <- RSpectra::svds(crossprod(EUk,Y[[1]])-tcrossprod(EUkU0,V_0k),r[1])
      tempL <- S$u
      tempR <- S$v
      V_ind[[1]] <- tcrossprod(tempR,tempL)
      for(k in 2:K){
        loc1 <- r0+sum(r[1:(k-1)])+1
        loc2 <- loc1 + r[k] -1
        loc3 <- sum(p[1:(k-1)])+1
        loc4 <- sum(p[1:k])

        # critical value
        EU0 <- EU[,(1:r0)]
        EUk <- EU[,(loc1:loc2)]
        EUkU0 <- n*covU[(loc1:loc2),(1:r0)]+crossprod(EUk,EU0) # r(k)*r0
        EU0U0 <- n*covU[(1:r0),(1:r0)]+crossprod(EU0, EU0) #r0*r0
        Vk <- V_ind[[k]] #p(k)*r(k)

        # V_0k
        V_0k <- (crossprod(Y[[k]],EU0)-Vk%*%EUkU0)%*%solve(EU0U0)
        V_joint[(loc3:loc4),] <- V_0k

        # V_k
        S <- RSpectra::svds(crossprod(EUk,Y[[k]])-tcrossprod(EUkU0,V_0k),r[k])
        tempL <- S$u
        tempR <- S$v
        V_ind[[k]] <- tcrossprod(tempR,tempL)
      }
      # end;

      # est se2
      loc1 <- r0+1
      loc2 <- loc1 + r[1]-1
      loc3 <- 1
      loc4 <- sum(p[1:1])
      covUcurrent <- covU[c(1:r0,loc1:loc2),c(1:r0,loc1:loc2)] #(r0 + r(k))*(r0+r(k))
      Ycurrent <- Y[[1]]
      EUcurrent <- cbind(EU[,(1:r0)], EU[,(loc1:loc2)])
      Vcurrent <- cbind(V_joint[(loc3:loc4),],V_ind[[1]])
      temp1 <- sum(diag(tcrossprod(Ycurrent,Ycurrent)))
      temp2 <- 2*sum(diag(tcrossprod(tcrossprod(EUcurrent,Vcurrent),Ycurrent)))
      temp3 <- n*sum(diag(crossprod(Vcurrent,Vcurrent)%*%covUcurrent))
      temp4 <- sum(diag(crossprod(EUcurrent,EUcurrent)%*%crossprod(Vcurrent,Vcurrent)))
      se2[1] <- (temp1-temp2+temp3+temp4)/(n*p[1])
      for(k in 2:K){
        loc1 <- r0+sum(r[1:(k-1)])+1
        loc2 <- loc1 + r[k]-1
        loc3 <- sum(p[1:(k-1)])+1
        loc4 <- sum(p[1:k])
        covUcurrent <- covU[c(1:r0,loc1:loc2),c(1:r0,loc1:loc2)] #(r0 + r(k))*(r0+r(k))
        Ycurrent <- Y[[k]]
        EUcurrent <- cbind(EU[,(1:r0)], EU[,(loc1:loc2)])
        Vcurrent <- cbind(V_joint[(loc3:loc4),],V_ind[[k]])
        temp1 <- sum(diag(tcrossprod(Ycurrent,Ycurrent)))
        temp2 <- 2*sum(diag(tcrossprod(tcrossprod(EUcurrent,Vcurrent),Ycurrent)))
        temp3 <- n*sum(diag(crossprod(Vcurrent,Vcurrent)%*%covUcurrent))
        temp4 <- sum(diag(crossprod(EUcurrent,EUcurrent)%*%crossprod(Vcurrent,Vcurrent)))
        se2[k] <- (temp1-temp2+temp3+temp4)/(n*p[k])
      }

      # est fX
      fX <- EstU_fun(X,EU,relation,sparsity)

      # est Sf0 and Sf
      f0X <- fX[,(1:r0)]
      EUcurrent <- EU[,(1:r0)]
      covUcurrent <- covU[(1:r0),(1:r0)]
      temp1 <- n*covUcurrent
      temp2 <- crossprod(EUcurrent,EUcurrent)
      temp3 <- crossprod(f0X,f0X)
      temp4 <- crossprod(f0X,EUcurrent)
      temp5 <- crossprod(EUcurrent,f0X)
      #Sf0 <- (temp1+temp2+temp3-temp4-temp5)/n #questionable!!!
      Sf0 <- diag(diag(temp1+temp2+temp3-temp4-temp5)/n,nrow = length(diag(temp1+temp2+temp3-temp4-temp5)/n)) #exactly follow the draft

      loc1 <- r0+1
      loc2 <- loc1+r[1]-1
      fkX <- fX[,(loc1:loc2)]
      EUcurrent <- EU[,(loc1:loc2)]
      covUcurrent <- covU[(loc1:loc2),(loc1:loc2)]
      temp1 <- n*covUcurrent
      temp2 <- crossprod(EUcurrent,EUcurrent)
      temp3 <- crossprod(fkX,fkX)
      temp4 <- crossprod(fkX,EUcurrent)
      temp5 <- crossprod(EUcurrent,fkX)
      Sf[[1]] <- diag(diag((temp1+temp2+temp3-temp4-temp5)/n))

      for(k in 2:K){
        loc1 <- r0+sum(r[1:(k-1)])+1
        loc2 <- loc1+r[k]-1
        fkX <- fX[,(loc1:loc2)]
        EUcurrent <- EU[,(loc1:loc2)]
        covUcurrent <- covU[(loc1:loc2),(loc1:loc2)]
        temp1 <- n*covUcurrent
        temp2 <- crossprod(EUcurrent,EUcurrent)
        temp3 <- crossprod(fkX,fkX)
        temp4 <- crossprod(fkX,EUcurrent)
        temp5 <- crossprod(EUcurrent,fkX)
        Sf[[k]] <- diag(diag((temp1+temp2+temp3-temp4-temp5)/n))
      }

      # Post Standardization for V0 (and Sf0 and B0)
      S1 <- RSpectra::svds(V_joint%*%tcrossprod(Sf0,V_joint),r0)
      Sf0 <- diag(S1$d,nrow=length(S1$d))
      V_joint_new <- S1$u
      fX[,(1:r0)] <- fX[,(1:r0)]%*%crossprod(V_joint,V_joint_new)
      V_joint <- V_joint_new

      # reorder columns of B, V, and rows/columns of Sf
      tempA <- sort(diag(Sf[[1]]), decreasing = T, index.return = T)
      temp <- tempA$x
      I <- tempA$ix
      Sf[[1]] <- diag(temp)
      loc1 <- r0+1
      loc2 <- loc1+r[1]-1
      temp <- fX[,(loc1:loc2)]
      fX[,(loc1:loc2)] <- temp[,I]
      V_ind[[1]] <- V_ind[[1]][,I]
      for(k in 2:K){
        tempA <- sort(diag(Sf[[k]]), decreasing = T, index.return = T)
        temp <- tempA$x
        I <- tempA$ix
        Sf[[k]] <- diag(temp)
        loc1 <- r0+sum(r[1:(k-1)])+1
        loc2 <- loc1+r[k]-1
        temp <- fX[,(loc1:loc2)]
        fX[,(loc1:loc2)] <- temp[,I]
        V_ind[[k]] <- V_ind[[k]][,I]
      }

      # get grandV
      grandV <- pracma::blkdiag(V_ind[[1]])
      for(k in 2:K){
        grandV <- pracma::blkdiag(grandV,V_ind[[k]])
      }
      grandV <- cbind(V_joint,grandV) # sum(p)*(r0+sum(r))

      # stopping rule
      diff <- PrinAngle(grandV, grandV_old)
      recdiff <- cbind(recdiff,diff)

      # draw loglik trend
      #loglik <- loglikelihood1(Y,fX,V_joint,V_ind,se2,Sf0,Sf)
      #diff_max <- loglik-maxloglik #insurance, avoid likelihood decrease
      #maxloglik <- max(maxloglik,loglik)

      # check
      #recloglik <- cbind(recloglik,loglik)
      #plots

    }


    # calculate EU
    grandSf0 <- t(diag(Sf0))
    grandse2 <- rep(1,p[1])*se2[1]
    grandSf <- t(diag(Sf[[1]]))
    grandV_temp <- pracma::blkdiag(V_ind[[1]])
    Delta_temp <- rep(1,r[1])*(1/se2[1])
    for(k in 2:K){
      grandse2 <- c(grandse2,rep(1,p[k])*se2[k])
      grandSf <- c(grandSf, t(diag(Sf[[k]])))
      grandV_temp <- pracma::blkdiag(grandV_temp, V_ind[[k]])
      Delta_temp <- c(Delta_temp, rep(1,r[k])*(1/se2[k]))
    }

    grandse2_inv <- 1/grandse2 # 1*sum(p)
    grandSf_inv <- 1/grandSf # 1*sum(r)
    grandSf0_inv <- 1/grandSf0 # 1*r0
    grandV <- cbind(V_joint, grandV_temp) # sum(p)*(r0+sum(r))
    Delta1 <- pracma::arrayfun("*", t(grandV), grandse2_inv)%*%grandV # [r0+sum(r)]*[r0+sum(r)]
    Delta2_inv <- solve(diag(c(grandSf0_inv,grandSf_inv))+Delta1)
    temp <- grandV%*%tcrossprod(Delta2_inv,grandV) # sum(p)*sum(p)
    SigmaY_inv <- diag(grandse2_inv)-pracma::arrayfun("*",pracma::arrayfun("*",temp,grandse2_inv),t(grandse2_inv))
    # sum(p)*sum(p), not diagonal because of common structure, diff from SupSVD
    VSigmaYinvV <- Delta1 - Delta1%*%Delta2_inv%*%Delta1 #(r0+sum(r))*(r0+sum(r))
    EU <- fX%*%(diag(sum(r)+r0)-pracma::arrayfun("*",VSigmaYinvV,c(grandSf0,grandSf)))+
      grandY%*%SigmaY_inv%*%pracma::arrayfun("*",grandV,c(grandSf0,grandSf))
    # n*(r0+sum(r)), conditional mean

    # Print convergence information
    if(niter < max_niter){
      cat("SIPCA_A1 converges after", pracma::num2str(niter), "iterations.")
    }
    else{
      cat("SIPCA_A1 NOT converge after", pracma::num2str(max_niter), "iterations!!! Final change in angle:", pracma::num2str(diff))
    }
    return(list(V_joint=V_joint, V_ind=V_ind, se2=se2, Sf0=Sf0, Sf=Sf, EU=EU, fX=fX))
  }

  #### B_np
  else if(type == "Bnp"){
    S1 <- RSpectra::svds(grandY,r0) #note: segments of V_joint_ini correspond to each data set may not be orthogonal as desired.
    U_joint <- S1$u
    D_joint <- diag(S1$d,nrow=length(S1$d))
    V_joint_ini <- S1$v
    U_joint <- U_joint%*%D_joint
    V_ind <- list()
    Sf <- list()

    loc3 <- 1
    loc4 <- sum(p[1:1])
    Ycurrent <- Y[[1]]-tcrossprod(U_joint,V_joint_ini[(loc3:loc4),])
    S2 <- RSpectra::svds(Ycurrent,r[1])
    U_k_ini <- S2$u
    D_k_ini <- diag(S2$d,nrow=length(S2$d))
    V_k_ini <- S2$v
    U_k <- U_k_ini%*%D_k_ini
    V_ind[[1]] <- V_k_ini
    grandU <- cbind(grandU,U_k)
    grandV <- pracma::blkdiag(V_ind[[1]])
    for(k in 2:K){
      loc3 <- sum(p[1:(k-1)])+1
      loc4 <- sum(p[1:k])
      Ycurrent <- Y[[k]]-tcrossprod(U_joint,V_joint_ini[(loc3:loc4),])
      S2 <- RSpectra::svds(Ycurrent,r[k])
      U_k_ini <- S2$u
      D_k_ini <- diag(S2$d,nrow=length(S2$d))
      V_k_ini <- S2$v
      U_k <- U_k_ini%*%D_k_ini
      V_ind[[k]] <- V_k_ini
      grandU <- cbind(grandU,U_k)
      grandV <- pracma::blkdiag(grandV,V_ind[[k]])
    }
    # postprocess V_joint_ini such that it follows our identifiability condition

    V_joint <- matrix(0,sum(p),r0)

    loc3 <- 1
    loc4 <- sum(p[1:1])
    V_joint_k <- GramSchmidt(as.matrix(V_joint_ini[loc3:loc4,]),V_ind[[1]])$Q1
    V_joint[loc3:loc4,] <- V_joint_k*(1/sqrt(K))
    for(k in 2:K){
      loc3 <- sum(p[1:(k-1)])+1
      loc4 <- sum(p[1:k])
      V_joint_k <- GramSchmidt(as.matrix(V_joint_ini[loc3:loc4,]),V_ind[[k]])$Q1
      V_joint[loc3:loc4,] <- V_joint_k*(1/sqrt(K))
    }
    grandV <- cbind(V_joint,grandV)
    grandU <- cbind(U_joint,grandU)
    fX <- EstU_fun(X,grandU,relation,sparsity)

    se2 <- rep(0,K)
    loc1 <- r0+1
    loc2 <- loc1+r[1]-1
    loc3 <- 1
    loc4 <- sum(p[1:1])
    U_k <- grandU[,loc1:loc2]
    fkX <- fX[,loc1:loc2]
    Sf[[1]] <- diag(apply((U_k-fkX)^2,2,stats::sd))
    temp <- tcrossprod(grandU, grandV)
    se2[1] <- (norm(Y[[1]]-temp[,loc3:loc4],type = "F")%*%norm(Y[[1]]-temp[,loc3:loc4],type = "F"))%*%solve(n*p[1])
    for(k in 2:K){
      loc1 <- r0+sum(r[1:(k-1)])+1
      loc2 <- loc1+r[k]-1
      loc3 <- sum(p[1:(k-1)])+1
      loc4 <- sum(p[1:k])
      U_k <- grandU[,loc1:loc2]
      fkX <- fX[,loc1:loc2]
      Sf[[k]] <- diag(apply((U_k-fkX)^2,2,stats::sd))
      temp <- tcrossprod(grandU, grandV)
      se2[k] <- (norm(Y[[k]]-temp[,loc3:loc4],type = "F")%*%norm(Y[[k]]-temp[,loc3:loc4],type = "F"))%*%solve(n*p[k])
    }
    Sf0 <- diag(apply((U_joint-fX[,1:r0])^2,2,stats::sd),nrow=length(apply((U_joint-fX[,1:r0])^2,2,stats::sd)))

    #disp('Initial est done!')
    loglik <- loglikelihood1(Y,fX,V_joint,V_ind,se2,Sf0,Sf)
    recloglik <- loglik
    maxloglik <- loglik

    niter <- 0
    diff = Inf
    recdiff <- NULL
    while(niter<=max_niter && abs(diff)>convg_thres){
      # && diff_max > tol_thres
      cat(sprintf('This is the %.d th iterations, the diff= %.4g \n',niter,diff))
      niter <- niter+1


      #record last iter
      grandV_old <- grandV

      # E step
      # some critical values
      grandse2 <- NULL #1*sum(p)
      grandSf <- t(diag(Sf0)) #1*(r0+sum(r))
      Delta_temp <- rep(1,r0)*sum(se2^(-1))/K
      for(k in 1:K){
        grandse2 <- c(grandse2,(rep(1,p[k])*se2[k]))
        grandSf <- c(grandSf,t(diag(Sf[[k]])))
        Delta_temp <- c(Delta_temp,rep(1,r[k])*(1/se2[k]))
      }
      grandse2_inv <- 1/grandse2 #1*sum(p)
      grandSf_inv <- 1/grandSf #1*(r0+sum(r))
      Delta1 <- Delta_temp
      # bsxfun("*",t(grandV), grandse2_inv)%*%grandV # 1*(r0+sum(r))
      Delta2 <- grandSf_inv+Delta1 #1*(r0+sum(r))
      temp <- pracma::arrayfun("*", grandV, (1/Delta2))%*%t(grandV) #sum(p)*sum(p)
      SigmaY_inv <- diag(grandse2_inv)-pracma::arrayfun("*",pracma::arrayfun("*", temp, grandse2_inv), t(grandse2_inv))#sum(p)*sum(p), not diagonal because of common structure, diff from SupSVD
      VSigmaYinvV <- Delta1 - Delta1*(1/Delta2)*Delta1 #1*(r0+sum(r))
      EU <- fX%*%diag(rep(1,sum(r)+r0)-VSigmaYinvV*grandSf)+grandY%*%SigmaY_inv%*%pracma::arrayfun("*", grandV, grandSf) #conditional mean
      covU <- grandSf - grandSf*VSigmaYinvV*grandSf#1*(r0+sum(r)), conditional variance (turns out to be a diagonal matrix)

      # M step
      # est V and se2
      loc1 <- r0+1
      loc2 <- loc1+r[1]-1
      loc3 <- 1
      loc4 <- sum(p[1:1])
      Ycurrent <- Y[[1]]
      EUcurrent <- cbind(EU[,1:r0], EU[,loc1:loc2])
      EUcurrent_star <- cbind((1/sqrt(K))*EU[,1:r0],EU[,loc1:loc2])

      # V
      S3 <- RSpectra::svds(crossprod(Ycurrent,EUcurrent_star),(r0+r[1]))
      tempL <- S3$u
      tempR <- S3$v
      Vcurrent_star <- tcrossprod(tempL,tempR) # should have orthonormal columns
      V_joint[loc3:loc4,] <- Vcurrent_star[,1:r0]*(1/sqrt(K))
      V_ind[[1]] <- Vcurrent_star[,(r0+1):(r0+r[1])]
      Vcurrent <- cbind(V_joint[loc3:loc4,],V_ind[[1]])
      # clear EUcurrent_star Vcurrent_star

      # se2 (directly from SupSVD formula)
      covUcurrent <- covU[c((1:r0),(loc1:loc2)),drop = "FLASE"]#1*(r0+r(k))
      temp1 <- sum(diag(tcrossprod(Ycurrent,Ycurrent)))
      temp2 <- 2*sum(diag(tcrossprod(EUcurrent,Vcurrent)%*%t(Ycurrent)))
      temp3 <- n*sum(diag(crossprod(Vcurrent,Vcurrent)%*%diag(covUcurrent)))
      temp4 <- sum(diag(crossprod(EUcurrent,EUcurrent)%*%crossprod(Vcurrent,Vcurrent)))
      se2[1] <- (temp1-temp2+temp3+temp4)/(n*p[1])
      for(k in 2:K){
        loc1 <- r0+sum(r[1:(k-1)])+1
        loc2 <- loc1+r[k]-1
        loc3 <- sum(p[1:(k-1)])+1
        loc4 <- sum(p[1:k])
        Ycurrent <- Y[[k]]
        EUcurrent <- cbind(EU[,1:r0], EU[,loc1:loc2])
        EUcurrent_star <- cbind((1/sqrt(K))*EU[,1:r0],EU[,loc1:loc2])

        # V
        S3 <- RSpectra::svds(crossprod(Ycurrent,EUcurrent_star),(r0+r[k]))
        tempL <- S3$u
        tempR <- S3$v
        Vcurrent_star <- tcrossprod(tempL,tempR) # should have orthonormal columns
        V_joint[loc3:loc4,] <- Vcurrent_star[,1:r0]*(1/sqrt(K))
        V_ind[[k]] <- Vcurrent_star[,(r0+1):(r0+r[k])]
        Vcurrent <- cbind(V_joint[loc3:loc4,],V_ind[[k]])
        # clear EUcurrent_star Vcurrent_star

        # se2 (directly from SupSVD formula)
        covUcurrent <- covU[c((1:r0),(loc1:loc2)),drop = "FLASE"]#1*(r0+r(k))
        temp1 <- sum(diag(tcrossprod(Ycurrent,Ycurrent)))
        temp2 <- 2*sum(diag(tcrossprod(EUcurrent,Vcurrent)%*%t(Ycurrent)))
        temp3 <- n*sum(diag(crossprod(Vcurrent,Vcurrent)%*%diag(covUcurrent)))
        temp4 <- sum(diag(crossprod(EUcurrent,EUcurrent)%*%crossprod(Vcurrent,Vcurrent)))
        se2[k] <- (temp1-temp2+temp3+temp4)/(n*p[k])
      }

      # est fX
      fX <- EstU_fun(X, EU,relation, sparsity)# estimate an n*(r0+sum(ri)) matrix

      # est Sf
      f0X <- fX[,1:r0]
      EUcurrent <- EU[,1:r0]
      covUcurrent <- covU[1:r0,drop = "FLASE"]
      temp1 <- n*diag(covUcurrent,nrow=length(covUcurrent))
      temp2 <- crossprod(EUcurrent,EUcurrent)
      temp3 <- crossprod(f0X,EUcurrent)
      temp4 <- crossprod(f0X,EUcurrent)
      temp5 <- crossprod(EUcurrent,f0X)
      Sf0 <- diag(diag(temp1+temp2+temp3-temp4-temp5)/n,nrow=length(diag(temp1+temp2+temp3-temp4-temp5)/n))

      loc1 <- r0+1
      loc2 <- loc1+r[1]-1
      fkX <- fX[,loc1:loc2]
      EUcurrent <- EU[,loc1:loc2]
      covUcurrent <- covU[loc1:loc2]
      temp1 <- n*diag(covUcurrent)
      temp2 <- crossprod(EUcurrent,EUcurrent)
      temp3 <- crossprod(fkX, fkX)
      temp4 <- crossprod(fkX,EUcurrent)
      temp5 <- crossprod(EUcurrent,fkX)
      Sf[[1]] <- diag(diag((temp1+temp2+temp3-temp4-temp5)/n))

      for(k in 2:K){
        loc1 <- r0+sum(r[1:(k-1)])+1
        loc2 <- loc1+r[k]-1
        fkX <- fX[,loc1:loc2]
        EUcurrent <- EU[,loc1:loc2]
        covUcurrent <- covU[loc1:loc2]
        temp1 <- n*diag(covUcurrent)
        temp2 <- crossprod(EUcurrent,EUcurrent)
        temp3 <- crossprod(fkX, fkX)
        temp4 <- crossprod(fkX,EUcurrent)
        temp5 <- crossprod(EUcurrent,fkX)
        Sf[[k]] <- diag(diag((temp1+temp2+temp3-temp4-temp5)/n))

        # reorder columns of B, V, and rows/columns of Sf
        tempA <- sort(diag(Sf[[k]]), decreasing = T, index.return = T)
        temp <- tempA$x
        I <- tempA$ix
        Sf[[k]] <- diag(temp)
        loc1 <- r0+sum(r[1:(k-1)])+1
        loc2 <- loc1+r[k]-1
        temp <- fX[,(loc1:loc2)]
        fX[,(loc1:loc2)] <- temp[,I]
        V_ind[[k]] <- V_ind[[k]][,I]
      }

      #calc grandV
      grandV <- NULL
      grandV <- pracma::blkdiag(V_ind[[1]])
      for(k in 2:K){
        grandV <- pracma::blkdiag(grandV, V_ind[[k]])
      }
      grandV <- cbind(V_joint,grandV) #sum(p)*(r0+sum(r))

      #stopping rule
      diff <- PrinAngle(grandV, grandV_old)
      recdiff <- cbind(recdiff,diff)

      #draw loglik trend
      loglik <- loglikelihood1(Y, fX, V_joint,V_ind,se2, Sf0,Sf)
      diff_max <- loglik-maxloglik #insurance, avoid likelihood decrease
      maxloglik <- max(maxloglik,loglik)

      #check
      recloglik <- cbind(recloglik, loglik)
      #par(mfrow=c(2,1))
      #plot(recloglik,type = "b",main = "likelihood")
      #plot(recdiff, type = "b", main = paste("angle=",pracma::num2str(diff)))
    }

    # calc EU
    grandse2 <- NULL#1*sum(p)
    grandSf <- t(diag(Sf0)) #1*(r0+sum(r))
    grandV_temp <- NULL #combination of loadings: first few long columns are joint loadings, subsequent diagonal blocks are individual loadings
    Delta_temp <- rep(1,r0)*sum(se2^(-1))/K

    grandse2 <- c(grandse2,(rep(1,p[1])*se2[1]))
    grandSf <- c(grandSf, t(diag(Sf[[1]])))
    grandV_temp <- pracma::blkdiag(V_ind[[1]])
    Delta_temp <- c(Delta_temp,rep(1,r[1])*(1/se2[1]))
    for(k in 2:K){
      grandse2 <- c(grandse2,(rep(1,p[k])*se2[k]))
      grandSf <- c(grandSf, t(diag(Sf[[k]])))
      grandV_temp <- pracma::blkdiag(grandV_temp,V_ind[[k]])
      Delta_temp <- c(Delta_temp,rep(1,r[k])*(1/se2[k]))
    }
    grandse2_inv <- 1/grandse2 #1*sum(p)
    grandSf_inv <- 1/grandSf #1*(r0+sum(r))
    grandV <- cbind(V_joint,grandV_temp) #sum(p)*(r0+sum(r))
    Delta1 <- Delta_temp #bsxfun("*", t(grandV), grandse2_inv)%*%grandV #1*(r0+sum(r))
    Delta2 <- grandSf_inv+Delta1 #1*(r0+sum(r))
    temp <- pracma::arrayfun("*", grandV, (1/Delta2))%*%t(grandV) #sum(p)*sum(p)
    SigmaY_inv <- diag(grandse2_inv)-pracma::arrayfun("*", pracma::arrayfun("*",temp,grandse2_inv), t(grandse2_inv)) #sum(p)*sum(p), not diagonal because of common structure, diff from SupSVD
    VSigmaYinV <- Delta1-Delta1*(1/Delta2)*Delta1 #1*(r0+sum(r))
    EU <- fX%*%diag(rep(1,sum(r)+r0)-VSigmaYinV*grandSf)+grandY%*%SigmaY_inv%*%pracma::arrayfun("*", grandV,grandSf) #conditional mean

    #Print convergence information
    if(niter<max_niter){
      cat("SIPCA_B1 converges after",pracma::num2str(niter),"iterations.")
    }else{
      cat("SIPCA_B1 NOT converge after", pracma::num2str(max_niter),"iterations!!! Final change in angle:",pracma::num2str(diff))
    }
    return(list(V_joint=V_joint, V_ind = V_ind, se2 = se2, Sf0 = Sf0, Sf = Sf, EU = EU, fX=fX))
  }


  ###AAA
  else if(type=="A"){
    S <- RSpectra::svds(grandY,r0)
    U_joint_ini <- S$u
    D_joint_ini <- diag(S$d,nrow=length(S$d))
    V_joint <- S$v
    # note: segments of V_joint_ini corresponding to each data set may not be orthogonal as desired
    U_joint <- U_joint_ini%*%D_joint_ini
    if(sparsity==1){
      B0 <- matrix(0,q,r0)
      for(i in 1:r0){
        fit <- glmnet::glmnet(X, U_joint[,i], alpha = 1, lambda.min.ratio = 0, standardize = TRUE)
        BIC_score <- n*stats::deviance(fit)+log(n)*fit$df
        ind <- which(BIC_score == min(BIC_score))
        B0[,i] <- fit$beta[,ind]
      }
    }else{# if sparsity = 0, no B-sparsity
      B0 <- solve(crossprod(X,X), t(X))%*%U_joint
    }
    Sf0 <- diag(apply((U_joint-X%*%B0)^2,2,stats::sd),nrow=length(apply((U_joint-X%*%B0)^2,2,stats::sd)))
    V_ind <- list()
    B <- list()
    for(k in 1:K){
      B[[k]] <- matrix(0,q,r[k])
    }
    Sf <- list()

    loc3 <- 1
    loc4 <- sum(p[1:1])
    Vcurrent <- V_joint[(loc3:loc4),]
    Ycurrent <- Y[[1]]-tcrossprod(U_joint,Vcurrent)
    S2 <- RSpectra::svds(Ycurrent,r[1])
    U_k_ini <- S2$u
    D_k_ini <- diag(S2$d,nrow=length(S2$d))
    V_k <- S2$v
    U_k <- U_k_ini%*%D_k_ini
    V_ind[[1]] <- V_k
    if (sparsity==1){
      for(i in 1:r[1]){
        fit <- glmnet::glmnet(X, U_k[,i], alpha = 1, lambda.min.ratio = 0, standardize = TRUE)
        BIC_score <- n*stats::deviance(fit)+log(n)*fit$df
        ind <- which(BIC_score == min(BIC_score))
        B[[1]][,i] <- fit$beta[,ind]
      }
    }else{# if sparsity = 0, no B-sparsity
      B[[1]] <- solve(crossprod(X,X),t(X))%*%U_k
    }
    Sf[[1]] <- diag(apply((U_k - X%*%B[[1]])^2,2,stats::sd))
    se2 <- rep(0,K)
    se2[1] <- norm(Y[[k]]-tcrossprod(U_k,V_k)-tcrossprod(U_joint,Vcurrent), type = "F")^2/(n*p[1])

    for(k in 2:K){
      loc3 <- sum(p[1:(k-1)])+1
      loc4 <- sum(p[1:k])
      Vcurrent <- V_joint[(loc3:loc4),]
      Ycurrent <- Y[[k]]-tcrossprod(U_joint,Vcurrent)
      S2 <- RSpectra::svds(Ycurrent,r[k])
      U_k_ini <- S2$u
      D_k_ini <- diag(S2$d,nrow=length(S2$d))
      V_k <- S2$v
      U_k <- U_k_ini%*%D_k_ini
      V_ind[[k]] <- V_k
      if (sparsity==1){
        for(i in 1:r[k]){
          fit <- glmnet::glmnet(X, U_k[,i], alpha = 1, lambda.min.ratio = 0, standardize = TRUE)
          BIC_score <- n*stats::deviance(fit)+log(n)*fit$df
          ind <- which(BIC_score == min(BIC_score))
          B[[k]][,i] <- fit$beta[,ind]
        }
      }else{# if sparsity = 0, no B-sparsity
        B[[k]] <- solve(crossprod(X,X),t(X))%*%U_k
      }
      Sf[[k]] <- diag(apply((U_k - X%*%B[[k]])^2,2,stats::sd))
      se2[k] <- norm(Y[[k]]-tcrossprod(U_k,V_k)-tcrossprod(U_joint,Vcurrent), type = "F")^2/(n*p[k])
    }
    grandV <- NULL
    grandV <- pracma::blkdiag(V_ind[[1]])
    for (k in 2:K){
      grandV <- pracma::blkdiag(grandV,V_ind[[k]])
    }
    grandV <- cbind(V_joint, grandV) # sum(p)*(r0+sum(r))

    #disp('Initial set done!')
    loglik <- loglikelihood(X,Y,B0,B,V_joint,V_ind,se2,Sf0,Sf)
    recloglik <- loglik
    maxloglik <- loglik

    niter <- 0
    diff <- Inf
    diff_max <- Inf
    tol_thres <- -10
    recdiff <- NULL
    while(niter<=max_niter && abs(diff)>convg_thres)#&& diff_max>tol_thres
    {
      cat(sprintf('This is the %.d th iterations, the diff= %.4g \n',niter,diff))
      niter <- niter + 1


      #record last iter
      #loglik_old <- loglik
      #V_joint_old <- V_joint
      grandV_old <-grandV

      # E step
      # some critical values
      grandse2 <- NULL #1*sum(p)
      grandSf0 <- t(diag(Sf0)) # 1*r0
      grandSf <- NULL #1*sum(r)
      #grandV_temp #combination of loadings: first few long columns are joint loadings, subsequent diagonal blocks are individual loadings
      grandB_temp <- NULL #concatenate B
      Delta_temp <- NULL #1*sum(r)
      for(k in 1:K){
        grandse2 <- c(grandse2,(rep(1,p[k])*se2[k]))
        grandSf <- c(grandSf,t(diag(Sf[[k]])))
        # grandV_temp = pracma::blkdiag(grandV_temp, V_ind[[k]])
        grandB_temp <- cbind(grandB_temp, B[[k]])
        Delta_temp <- c(Delta_temp, rep(1,r[k])*(1/se2[k]))
      }
      grandse2_inv <- 1/grandse2 #1*sum(p)
      grandSf_inv <- 1/grandSf # 1*sum(r)
      grandSf0_inv <- 1/grandSf0 #1*r0
      # grandV = cbind(V_joint, grandV_temp) #sum(p)*(r0+sum(r))
      grandB <- cbind(B0,grandB_temp) # q*(r0+sum(r))
      Delta1 <- pracma::arrayfun("*", t(grandV), grandse2_inv)%*%grandV #[r0+sum(r)]*[r0+sum(r)]

      Delta2_inv <- solve(diag(c(grandSf0_inv,grandSf_inv))+Delta1)
      temp <- grandV%*%tcrossprod(Delta2_inv,grandV)
      SigmaY_inv <- diag(grandse2_inv) - pracma::arrayfun("*", pracma::arrayfun("*",temp, grandse2_inv), t(grandse2_inv)) #sum(p)*sum(p), not diagonal because of common structure, diff from SupSVD
      VSigmaYinvV <- Delta1 - Delta1%*%Delta2_inv%*%Delta1 #(r0+sum(r))*(r0+sum(r))
      EU <- X%*%grandB%*%(diag(sum(r)+r0)-pracma::arrayfun("*", VSigmaYinvV,c(grandSf0,grandSf)))+grandY%*%SigmaY_inv%*%pracma::arrayfun("*", grandV, c(grandSf0, grandSf)) # n*(r0+sum(r)), conditional mean
      covU <- diag(c(grandSf0,grandSf))-pracma::arrayfun("*", pracma::arrayfun("*", VSigmaYinvV, c(grandSf0,grandSf)), t(c(grandSf0, grandSf))) #(r0+sum(r))*(r0+sum(r))

      # M step
      # V and se2
      # for iter = 1:3 # alternate between V_joint and V_ind

      loc1 <- r0 + +1
      loc2 <- loc1 + r[1]-1
      loc3 <- 1
      loc4 <- sum(p[1:1])

      # critical values
      EU0 <- EU[,(1:r0)]
      EUk <- EU[,(loc1:loc2)]
      EUkU0 <- n*covU[(loc1:loc2), (1:r0)]+crossprod(EUk, EU0) #r(k)*r0
      EU0U0 <- n*covU[(1:r0), (1:r0)]+crossprod(EU0,EU0) #r0*r0
      Vk <- V_ind[[1]]

      #V_0k
      V_0k <- (crossprod(Y[[1]], EU0)-Vk%*%EUkU0)%*%solve(EU0U0)
      V_joint[(loc3:loc4),] <- V_0k

      #V_k
      S3 <- RSpectra::svds(crossprod(EUk, Y[[1]])-tcrossprod(EUkU0, V_0k),r[1])
      V_ind[[1]] <- tcrossprod(S3$v, S3$u)

      for(k in 2:K){
        loc1 <- r0 + sum(r[1:(k-1)])+1
        loc2 <- loc1 + r[k]-1
        loc3 <- sum(p[1:(k-1)])+1
        loc4 <- sum(p[1:k])

        # critical values
        EU0 <- EU[,(1:r0)]
        EUk <- EU[,(loc1:loc2)]
        EUkU0 <- n*covU[(loc1:loc2), (1:r0)]+crossprod(EUk, EU0) #r(k)*r0
        EU0U0 <- n*covU[(1:r0), (1:r0)]+crossprod(EU0,EU0) #r0*r0
        Vk <- V_ind[[k]]

        #V_0k
        V_0k <- (crossprod(Y[[k]], EU0)-Vk%*%EUkU0)%*%solve(EU0U0)
        V_joint[(loc3:loc4),] <- V_0k

        #V_k
        S3 <- RSpectra::svds(crossprod(EUk, Y[[k]])-tcrossprod(EUkU0, V_0k),r[k])
        V_ind[[k]] <- tcrossprod(S3$v, S3$u)
      }

      # est se2
      loc1 <- r0+1
      loc2 <- loc1 + r[1]-1
      loc3 <- 1
      loc4 <- sum(p[1:1])
      covUcurrent <- covU[c(1:r0, loc1:loc2), c(1:r0, loc1:loc2)] #(r0+r(k))*(r0+r(k))
      Ycurrent <- Y[[1]]
      EUcurrent <- cbind(EU[,(1:r0)], EU[,(loc1:loc2)])
      Vcurrent <- cbind(V_joint[(loc3:loc4),], V_ind[[1]])
      temp1 <- sum(diag(tcrossprod(Ycurrent, Ycurrent)))
      temp2 <- 2*sum(diag(tcrossprod(tcrossprod(EUcurrent,Vcurrent),Ycurrent)))
      temp3 <- n*sum(diag(crossprod(Vcurrent,Vcurrent)%*%covUcurrent))
      temp4 <- sum(diag(crossprod(EUcurrent,EUcurrent)%*%crossprod(Vcurrent,Vcurrent)))
      se2[1] <- (temp1-temp2+temp3+temp4)/(n*p[1])

      for(k in 2:K){
        loc1 <- r0+sum(r[1:(k-1)])+1
        loc2 <- loc1 + r[k]-1
        loc3 <- sum(p[1:(k-1)])+1
        loc4 <- sum(p[1:k])
        covUcurrent <- covU[c(1:r0, loc1:loc2), c(1:r0, loc1:loc2)] #(r0+r(k))*(r0+r(k))
        Ycurrent <- Y[[k]]
        EUcurrent <- cbind(EU[,(1:r0)], EU[,(loc1:loc2)])
        Vcurrent <- cbind(V_joint[(loc3:loc4),], V_ind[[k]])
        temp1 <- sum(diag(tcrossprod(Ycurrent, Ycurrent)))
        temp2 <- 2*sum(diag(tcrossprod(tcrossprod(EUcurrent,Vcurrent),Ycurrent)))
        temp3 <- n*sum(diag(crossprod(Vcurrent,Vcurrent)%*%covUcurrent))
        temp4 <- sum(diag(crossprod(EUcurrent,EUcurrent)%*%crossprod(Vcurrent,Vcurrent)))
        se2[k] <- (temp1-temp2+temp3+temp4)/(n*p[k])
      }

      # B and Sf
      EUcurrent <- EU[,(1:r0)]
      if (sparsity==1){
        for(i in 1:r0){
          fit <- glmnet::glmnet(X,EUcurrent[,i],lambda.min.ratio = 0, standardize = TRUE)
          BIC_score <- n*stats::deviance(fit)+log(n)*fit$df
          ind <- which(BIC_score == min(BIC_score))
          B0[,i] <- fit$beta[,ind]
        }
      }else{
        B0 <- solve(crossprod(X,X),t(X))%*%EUcurrent
      }
      covUcurrent <- covU[(1:r0), (1:r0)]
      temp1 <- n*covUcurrent
      temp2 <- crossprod(EUcurrent,EUcurrent)
      temp3 <- crossprod(B0, crossprod(X,X))%*%B0
      temp4 <- t(B0)%*%crossprod(X,EUcurrent)
      temp5 <- crossprod(EUcurrent, X)%*%B0
      Sf0 <- (temp1+temp2+temp3-temp4-temp5)/n #not diagonal


      loc1 <- r0+1
      loc2 <- loc1+r[1]-1
      EUcurrent <- EU[,(loc1:loc2)]
      if(sparsity==1){
        for(i in 1:r[1]){
          fit <- glmnet::glmnet(X,EUcurrent[,i],lambda.min.ratio = 0, standardize = TRUE)
          BIC_score <- n*stats::deviance(fit)+log(n)*fit$df
          ind <- which(BIC_score == min(BIC_score))
          B[[1]][,i] <- fit$beta[,ind]
        }
      }else{ # no B-sparsity
        B[[1]] <- solve(crossprod(X,X),t(X))%*%EUcurrent
      }
      covUcurrent <- covU[(loc1:loc2), (loc1:loc2)]
      temp1 <- n*covUcurrent
      temp2 <- crossprod(EUcurrent,EUcurrent)
      temp3 <- crossprod(B[[1]], crossprod(X,X))%*%B[[1]]
      temp4 <- t(B[[1]])%*%crossprod(X,EUcurrent)
      temp5 <- crossprod(EUcurrent, X)%*%B[[1]]
      Sf[[1]] <- diag(diag((temp1+temp2+temp3-temp4-temp5)/n))

      for(k in 2:K){
        loc1 <- r0+sum(r[1:(k-1)])+1
        loc2 <- loc1+r[k]-1
        EUcurrent <- EU[,(loc1:loc2)]
        if(sparsity==1){
          for(i in 1:r[k]){
            fit <- glmnet::glmnet(X,EUcurrent[,i],lambda.min.ratio = 0, standardize = TRUE)
            BIC_score <- n*stats::deviance(fit)+log(n)*fit$df
            ind <- which(BIC_score == min(BIC_score))
            B[[k]][,i] <- fit$beta[,ind]
          }
        }else{ # no B-sparsity
          B[[k]] <- solve(crossprod(X,X),t(X))%*%EUcurrent
        }
        covUcurrent <- covU[(loc1:loc2), (loc1:loc2)]
        temp1 <- n*covUcurrent
        temp2 <- crossprod(EUcurrent,EUcurrent)
        temp3 <- crossprod(B[[k]], crossprod(X,X))%*%B[[k]]
        temp4 <- t(B[[k]])%*%crossprod(X,EUcurrent)
        temp5 <- crossprod(EUcurrent, X)%*%B[[k]]
        Sf[[k]] <- diag(diag((temp1+temp2+temp3-temp4-temp5)/n))
      }

      #Post Standardization for V0 (and Sf0 and B0)
      S4 <- RSpectra::svds(V_joint%*%tcrossprod(Sf0, V_joint),r0)
      V_joint_new <- S4$u
      Sf0 <- diag(S4$d,nrow=length(S4$d))
      V_joint_new <- pracma::arrayfun("*", V_joint_new, sign(V_joint_new[1,])) # make sure no sign alternation
      B0 <- tcrossprod(B0,V_joint)%*%V_joint_new #will sacrifice any sparsity in B0
      V_joint <- V_joint_new #reorder columns of B, V, and rows/columns of Sf
      for(k in 1:K){
        tempA <- sort(diag(Sf[[k]]), decreasing = T, index.return = T)
        temp <- tempA$x
        I <- tempA$ix
        Sf[[k]] <- diag(temp)
        B[[k]] <- B[[k]][,I]
        V_ind[[k]] <- V_ind[[k]][,I]
      }
      # get grandV
      grandV <- pracma::blkdiag(V_ind[[1]])
      for(k in 2:K){
        grandV <- pracma::blkdiag(grandV, V_ind[[k]])
      }
      grandV <- cbind(V_joint, grandV) #sum(p)*(r0+sum(r))

      #Stopping rule
      diff <- PrinAngle(grandV,grandV_old)
      recdiff <- cbind(recdiff,diff)

      #draw loglik trend
      loglik <- loglikelihood(X,Y,B0,B,V_joint,V_ind,se2,Sf0,Sf)
      diff_max <- loglik-maxloglik # insurance, avoid likelihood decrease
      maxloglik <- max(maxloglik,loglik)

    }

    #calculate EU
    grandse2 <- NULL #1*sum(p)
    grandSf0 <- t(diag(Sf0)) #1*r0, r0>1
    grandSf <- NULL#1*sum(r)
    grandV_temp <- NULL #combination of loadings: first few long columns are joint loadings, subsequent diagonal blocks are individual loadings.
    grandB_temp <- NULL #concatenate B
    Delta_temp <- NULL #1*sum(r)

    grandse2 <- (rep(1,p[1])*se2[1])
    grandSf <- t(diag(Sf[[1]]))
    grandV_temp <- pracma::blkdiag(V_ind[[1]])
    grandB_temp <- B[[1]]
    Delta_temp <- rep(1,r[1])*(1/se2[1])
    for(k in 2:K){
      grandse2 <- c(grandse2, (rep(1,p[k])*se2[k]))
      grandSf <- c(grandSf, t(diag(Sf[[k]])))
      grandV_temp <- pracma::blkdiag(grandV_temp,V_ind[[k]])
      grandB_temp <- cbind(grandB_temp,B[[k]])
      Delta_temp <- c(Delta_temp,rep(1,r[k])*(1/se2[k]))
    }
    grandse2_inv <- 1/grandse2 #1*sum(p)
    grandSf_inv <- 1/grandSf #1*sum(r)
    grandSf0_inv <- 1/grandSf0 #1*r0
    grandV <- cbind(V_joint, grandV_temp) #sum(p)*(r0+sum(r))
    grandB <- cbind(B0,grandB_temp) #q*(r0+sum(r))
    Delta1 <- pracma::arrayfun("*", t(grandV), grandse2_inv)%*%grandV # [r0+sum(r)]*[r0+sum(r)]
    Delta2_inv <- solve(diag(c(grandSf0_inv,grandSf_inv))+Delta1)
    temp <- grandV%*%tcrossprod(Delta2_inv,grandV) #sum(p)*sum(p)
    SigmaY_inv <- diag(grandse2_inv)-pracma::arrayfun("*",pracma::arrayfun("*",temp,grandse2_inv),t(grandse2_inv)) #sum(p)*sum(p), not diagonal because of common structure, diff from SupSVD
    VSigmaYinvV <- Delta1 - Delta1%*%Delta2_inv%*%Delta1 #(r0+sum(r))*(r0+sum(r))
    EU <- X%*%grandB%*%(diag(sum(r)+r0)-pracma::arrayfun("*",VSigmaYinvV,c(grandSf0,grandSf)))+grandY%*%SigmaY_inv%*%pracma::arrayfun("*",grandV,c(grandSf0,grandSf)) #n*(r0+sum(r)), conditional mean

    #Print convergence information
    if(niter < max_niter){
      cat("SIPCA_A converges after",pracma::num2str(niter), "iterations.")
    }else{
      cat("SIPCA_A NOT converge after", pracma::num2str(max_niter), "iterations!!! Final change in angle:",pracma::num2str(diff))
    }
    return(list(B0=B0,B=B,V_join=V_joint,V_ind=V_ind, se2=se2, Sf0=Sf0,Sf=Sf, EU=EU))
  }

  ###BBB
  else if(type=="B"){
    S <- RSpectra::svds(grandY,r0)
    U_joint_ini <- S$u
    D_joint_ini <- diag(S$d,nrow=length(S$d))
    V_joint_ini <- S$v
    # note: segments of V_joint_ini corresponding to each data set may not be orthogonal as desired
    U_joint_ini <- U_joint_ini%*%D_joint_ini
    B0 <- matrix(0,q,r0)
    if (sparsity ==1){
      for (i in 1:r0){
        # tic
        #Attention: lasso fcn alwayas center X,Y, and coefficients are calcualted based on centered X and Y. The function will return a separate column for intercept; if we center data outside the function, the intercept will be nearly 0.
        fit <- glmnet::glmnet(X, U_joint_ini[,i], alpha = 1, lambda.min.ratio = 0, standardize = TRUE)
        BIC_score <- n*stats::deviance(fit)+log(n)*fit$df
        ind <- which(BIC_score == min(BIC_score))
        B0[,i] <- fit$beta[,ind]
        #toc
      }
    }else{
      # if sparsity = 0, no B-sparsity
      B0 <- solve(crossprod(X,X), t(X))%*%U_joint_ini
    }
    Sf0 <- diag(apply((U_joint_ini-X%*%B0)^2,2,stats::sd),nrow=length(apply((U_joint_ini-X%*%B0)^2,2,stats::sd)))
    V_ind <- list()
    B <- list()
    for(k in 1:K){
      B[[k]] <- matrix(0,q,r[k])
    }
    Sf <- list()

    se2 <- rep(0,K)
    loc3 <- 1
    loc4 <- sum(p[1:1])
    Ycurrent <- Y[[1]]-tcrossprod(U_joint_ini,V_joint_ini[loc3:loc4,])
    S2 <- RSpectra::svds(Ycurrent,r[1])
    U_k_ini <- S2$u
    D_k_ini <- diag(S2$d,nrow=length(S2$d))
    V_k_ini <- S2$v
    U_k_ini <- U_k_ini%*%D_k_ini
    V_ind[[1]] <- V_k_ini
    if (sparsity==1){
      for(i in 1:r[1]){
        fit <- glmnet::glmnet(X, U_k_ini[,i], alpha = 1, lambda.min.ratio = 0, standardize = TRUE)
        BIC_score <- n*stats::deviance(fit)+log(n)*fit$df
        ind <- which(BIC_score == min(BIC_score))
        B[[1]][,i] <- fit$beta[,ind]
      }
    }else{
      # if sparsity = 0, no B-sparsity
      B[[1]] <- solve(crossprod(X,X),t(X))%*%U_k_ini
    }
    Sf[[1]] <- diag(apply((U_k_ini - X%*%B[[1]])^2,2,stats::sd))
    se2[1] <- norm(Ycurrent-tcrossprod(U_k_ini,V_k_ini), type = "F")^2/(n*p[1])
    for(k in 2:K){
      loc3 <- sum(p[1:(k-1)])+1
      loc4 <- sum(p[1:k])
      Ycurrent <- Y[[k]]-tcrossprod(U_joint_ini,V_joint_ini[loc3:loc4,])
      S2 <- RSpectra::svds(Ycurrent,r[k])
      U_k_ini <- S2$u
      D_k_ini <- diag(S2$d,nrow=length(S2$d))
      V_k_ini <- S2$v
      U_k_ini <- U_k_ini%*%D_k_ini
      V_ind[[k]] <- V_k_ini
      if (sparsity==1){
        for(i in 1:r[k]){
          fit <- glmnet::glmnet(X, U_k_ini[,i], alpha = 1, lambda.min.ratio = 0, standardize = TRUE)
          BIC_score <- n*stats::deviance(fit)+log(n)*fit$df
          ind <- which(BIC_score == min(BIC_score))
          B[[k]][,i] <- fit$beta[,ind]
        }
      }else{
        # if sparsity = 0, no B-sparsity
        B[[k]] <- solve(crossprod(X,X),t(X))%*%U_k_ini
      }
      Sf[[k]] <- diag(apply((U_k_ini - X%*%B[[k]])^2,2,stats::sd))
      se2[k] <- norm(Ycurrent-tcrossprod(U_k_ini,V_k_ini), type = "F")^2/(n*p[k])
    }

    #postprocess V_joint_ini such that it follows our identifiavility condition
    V_joint <- matrix(0,sum(p),r0)

    loc3 <- 1
    loc4 <- sum(p[1:1])
    V_joint_k <- GramSchmidt(as.matrix(V_joint_ini[loc3:loc4,]),V_ind[[1]])$Q1
    V_joint[loc3:loc4,] <- V_joint_k*(1/sqrt(K))
    for (k in 2:K){
      loc3 <- sum(p[1:(k-1)])+1
      loc4 <- sum(p[1:k])
      V_joint_k <- GramSchmidt(as.matrix(V_joint_ini[loc3:loc4,]),V_ind[[k]])$Q1
      V_joint[loc3:loc4,] <- V_joint_k*(1/sqrt(K))
    }
    grandV <- NULL

    grandV <- pracma::blkdiag(V_ind[[1]])
    for(k in 2:K){
      grandV <- pracma::blkdiag(grandV,V_ind[[k]])
    }
    grandV <- cbind(V_joint, grandV) # sum(p)*(r0+sum(r))


    #disp('Initial set done!')
    loglik <- loglikelihood(X,Y,B0,B,V_joint,V_ind,se2,Sf0,Sf)
    recloglik <- loglik
    maxloglik <- loglik

    niter <- 0
    diff <- Inf
    diff_max <- Inf
    tol_thres <- -10
    recdiff <- NULL

    recdiff_grandV <- NULL
    recdiff_V_joint <- NULL
    recdiff_V1 <- NULL
    recdiff_V2 <- NULL

    while(niter<=max_niter && abs(diff)>convg_thres)#&& diff_max>tol_thres
    {
      cat(sprintf('This is the %.d th iterations, the diff= %.4g \n',niter,diff))
      niter <- niter+1


      #record last iter
      #loglik_old <- loglik
      grandV_old <- grandV
      V1_old <- V_ind[[1]]
      V2_old <- V_ind[[2]]

      #E step
      grandse2 <- NULL #1*sum(p)
      grandSf <- t(diag(Sf0)) # 1*(r0+sum(r))
      # grandV_temp = NULL #combination of loadings: first few long columns are joint loadings, subsequent diagonal blocks are individual loadings
      grandB_temp <- NULL #concatenate B
      Delta_temp <- rep(1,r0)*sum(se2^(-1))/K
      for(k in 1:K){
        grandse2 <- c(grandse2, (rep(1,p[k])*se2[k]))
        grandSf <- c(grandSf,t(diag(Sf[[k]])))
        # grandV_temp <- pracma::blkdiag(grandV_temp, V_ind[k])
        grandB_temp <- cbind(grandB_temp, B[[k]])
        Delta_temp <- c(Delta_temp, rep(1,r[k])*(1/se2[k]))
      }
      grandse2_inv <- 1/grandse2 #1*sum(p)
      grandSf_inv <- 1/grandSf #1*(r0+sum(r))
      # grandV <- cbind(V_joint, grandV_temp) #sum(p)*(r0+sum(r))
      grandB <- cbind(B0, grandB_temp)#q*(r0+sum(r))
      Delta1 <- Delta_temp #bsxfun(@times, grandV', grandse2_inv)*grandV; % 1*(r0+sum(r))
      Delta2 <- grandSf_inv+Delta1 #1*(r0+sum(r))
      temp <- tcrossprod(pracma::arrayfun("*",grandV, (1/Delta2)),grandV) #sum(p)*sum(p)
      SigmaY_inv <- diag(grandse2_inv)-pracma::arrayfun("*",pracma::arrayfun("*",temp, grandse2_inv),t(grandse2_inv)) #sum(p)*sum(p), not diagonal because of common structure, diff from SupSVD
      VSigmaYinvV <- Delta1 - Delta1*(1/Delta2)*Delta1 #1*(r0+sum(r))
      EU <- X%*%pracma::arrayfun("*",grandB, (rep(1,sum(r)+r0)-VSigmaYinvV*grandSf))+grandY%*%SigmaY_inv%*%pracma::arrayfun("*",grandV,grandSf) # conditional mean
      covU <- grandSf-grandSf*VSigmaYinvV*grandSf #1*(r0+sum(r)), conditional variance (turns out to be a diagonal matrix)


      # M step
      # V and se2
      loc1 <- r0+1
      loc2 <- loc1+r[1]-1
      loc3 <- 1
      loc4 <- sum(p[1:1])
      Ycurrent <- Y[[1]]
      EUcurrent <- cbind(EU[,1:r0],EU[,loc1:loc2])
      EUcurrent_star <- cbind((1/sqrt(K))*EU[,1:r0],EU[,loc1:loc2])

      #V
      S3 <- RSpectra::svds(crossprod(Ycurrent,EUcurrent_star),(r0+r[1]))
      tempL <- S3$u
      tempR <- S3$v
      Vcurrent_star <- tcrossprod(tempL,tempR) # should have orthonormal columns
      V_joint[loc3:loc4,] <- Vcurrent_star[,1:r0,drop=F]*(1/sqrt(K))
      V_ind[[1]] <- Vcurrent_star[,(r0+1):(r0+r[1])]
      Vcurrent <- cbind(V_joint[loc3:loc4,],V_ind[[1]])
      # clear EUcurrent_star Vcurrent_star

      # se2 (directly from SupSVD formula)
      covUcurrent <- covU[c(1:r0,loc1:loc2)] #1*(r0+r[k])
      temp1 <- sum(diag(tcrossprod(Ycurrent,Ycurrent)))
      temp2 <- 2*sum(diag(tcrossprod(EUcurrent,Vcurrent)%*%t(Ycurrent)))
      temp3 <- n*sum(diag((crossprod(Vcurrent,Vcurrent))%*%diag(covUcurrent)))
      temp4 <- sum(diag((crossprod(EUcurrent,EUcurrent))%*%(crossprod(Vcurrent,Vcurrent))))
      se2[1] <- (temp1 - temp2+temp3+temp4)/(n*p[1])

      for (k in 2:K){
        loc1 <- r0+sum(r[1:(k-1)])+1
        loc2 <- loc1+r[k]-1
        loc3 <- sum(p[1:(k-1)])+1
        loc4 <- sum(p[1:k])
        Ycurrent <- Y[[k]]
        EUcurrent <- cbind(EU[,1:r0],EU[,loc1:loc2])
        EUcurrent_star <- cbind((1/sqrt(K))*EU[,1:r0],EU[,loc1:loc2])

        #V
        S3 <- RSpectra::svds(crossprod(Ycurrent,EUcurrent_star),(r0+r[k]))
        tempL <- S3$u
        tempR <- S3$v
        Vcurrent_star <- tcrossprod(tempL,tempR) # should have orthonormal columns
        V_joint[loc3:loc4,] <- Vcurrent_star[,1:r0,drop=F]*(1/sqrt(K))
        V_ind[[k]] <- Vcurrent_star[,(r0+1):(r0+r[k])]
        Vcurrent <- cbind(V_joint[loc3:loc4,],V_ind[[k]])
        # clear EUcurrent_star Vcurrent_star

        # se2 (directly from SupSVD formula)
        covUcurrent <- covU[c(1:r0,loc1:loc2)] #1*(r0+r[k])
        temp1 <- sum(diag(tcrossprod(Ycurrent,Ycurrent)))
        temp2 <- 2*sum(diag(tcrossprod(EUcurrent,Vcurrent)%*%t(Ycurrent)))
        temp3 <- n*sum(diag((crossprod(Vcurrent,Vcurrent))%*%diag(covUcurrent)))
        temp4 <- sum(diag((crossprod(EUcurrent,EUcurrent))%*%(crossprod(Vcurrent,Vcurrent))))
        se2[k] <- (temp1 - temp2+temp3+temp4)/(n*p[k])
      }

      # B and Sf
      EUcurrent <- EU[,1:r0]
      if(sparsity ==1){
        cat("Start estimating B0...")
        for(i in 1:r0){
          fit <- glmnet::glmnet(X,EUcurrent[,i],lambda.min.ratio = 0, standardize = TRUE,alpha=1)
          BIC_score <- n*stats::deviance(fit)+log(n)*fit$df
          ind <- which(BIC_score == min(BIC_score))
          B0[,i] <- fit$beta[,ind]
        }
        cat("Finish estimating B0!")
      }else{
        # if sparsity = 0, no B-sparsity
        B0 <- solve(crossprod(X,X),t(X))%*%EUcurrent
      }
      covUcurrent <- covU[1:r0]
      temp1 <- n*diag(covUcurrent,nrow=length(covUcurrent))
      temp2 <- crossprod(EUcurrent,EUcurrent)
      temp3 <- crossprod(B0,crossprod(X,X))%*%B0
      temp4 <- t(B0)%*%crossprod(X,EUcurrent)
      temp5 <- crossprod(EUcurrent,X)%*%B0
      Sf0 <- diag(diag(temp1+temp2+temp3-temp4-temp5)/n,nrow=length(diag(temp1+temp2+temp3-temp4-temp5)/n))

      loc1 <- r0+1
      loc2 <- loc1+r[1]-1
      EUcurrent <- EU[,loc1:loc2]
      if(sparsity == 1){
        #cat("start estimating B", pracma::num2str(k),"...")
        for(i in 1:r[1]){
          fit <- glmnet::glmnet(X,EUcurrent[,i],lambda.min.ratio = 0, standardize = TRUE,alpha=1)
          BIC_score <- n*stats::deviance(fit)+log(n)*fit$df
          ind <- which(BIC_score == min(BIC_score))
          B[[1]][,i] <- fit$beta[,ind]
        }
      }else{ # no B-sparsity
        B[[1]] <- solve(crossprod(X,X),t(X))%*%EUcurrent
      }
      covUcurrent <- covU[loc1:loc2]
      temp1 <- n*diag(covUcurrent)
      temp2 <- crossprod(EUcurrent,EUcurrent)
      temp3 <- crossprod(B[[1]],crossprod(X,X))%*%B[[1]]
      temp4 <- t(B[[1]])%*%crossprod(X,EUcurrent)
      temp5 <- crossprod(EUcurrent,X)%*%B[[1]]
      Sf[[1]] <- diag(diag(temp1+temp2+temp3-temp4-temp5)/n)
      for(k in 2:K){
        loc1 <- r0+sum(r[1:(k-1)])+1
        loc2 <- loc1+r[k]-1
        EUcurrent <- EU[,loc1:loc2]
        if(sparsity == 1){
          #cat("start estimating B", pracma::num2str(k),"...")
          for(i in 1:r[k]){
            fit <- glmnet::glmnet(X,EUcurrent[,i],lambda.min.ratio = 0, standardize = TRUE,alpha=1)
            BIC_score <- n*stats::deviance(fit)+log(n)*fit$df
            ind <- which(BIC_score == min(BIC_score))
            B[[k]][,i] <- fit$beta[,ind]
          }
        }else{ # no B-sparsity
          B[[k]] <- solve(crossprod(X,X),t(X))%*%EUcurrent
        }
        covUcurrent <- covU[loc1:loc2]
        temp1 <- n*diag(covUcurrent)
        temp2 <- crossprod(EUcurrent,EUcurrent)
        temp3 <- crossprod(B[[k]],crossprod(X,X))%*%B[[k]]
        temp4 <- t(B[[k]])%*%crossprod(X,EUcurrent)
        temp5 <- crossprod(EUcurrent,X)%*%B[[k]]
        Sf[[k]] <- diag(diag(temp1+temp2+temp3-temp4-temp5)/n)
      }

      #reorder columns of B,V,and rows/columns of Sf
      tempA <- sort(diag(Sf0), decreasing = T, index.return = T)
      temp <- tempA$x
      I <- tempA$ix
      Sf0 <- diag(temp)
      B0 <- B0[,I]
      V_joint <- V_joint[,I,drop=F]
      for(k in 1:K){
        tempB <- sort(diag(Sf[[k]]), decreasing = T, index.return = T)
        temp <- tempB$x
        I <- tempB$ix
        Sf[[k]] <- diag(temp)
        B[[k]]<-B[[k]][,I]
        V_ind[[k]] <- V_ind[[k]][,I]
      }

      #calc grandV
      grandV <- NULL
      grandV<-pracma::blkdiag(V_ind[[1]])
      for(k in 2:K){
        grandV<-pracma::blkdiag(grandV,V_ind[[k]])
      }
      grandV <- cbind(V_joint,grandV) #sum(p)*(r0+sum(r))

      #stopping rule
      diff <- PrinAngle(grandV,grandV_old)
      recdiff <- cbind(recdiff,diff)
      #diff

      #draw loglik trend
      loglik <- loglikelihood(X,Y,B0,B,V_joint,V_ind,se2,Sf0,Sf)
      diff_max<- loglik-maxloglik #insurance, avoid likelihood decrease
      maxloglik <- max(maxloglik,loglik)
    }

    #calc EU
    grandse2 <- NULL #1*sum(p)
    grandSf <- t(diag(Sf0)) #1*(r0+sum(r))
    grandV_temp <- NULL #combination of loadings: first few long columns are joint loading, subsequent diagonal blocks are individual loadings
    grandB_temp <- NULL #concatenate B
    Delta_temp <- rep(1,r0)*sum(se2^(-1))/K

    grandse2 <- c(grandse2,(rep(1,p[1])*se2[1]))
    grandSf <- c(grandSf,t(diag(Sf[[1]])))
    grandV_temp <- pracma::blkdiag(V_ind[[1]])
    grandB_temp <- cbind(grandB_temp,B[[1]])
    Delta_temp <- c(Delta_temp,rep(1,r[1])*(1/se2[1]))
    for(k in 2:K){
      grandse2 <- c(grandse2,(rep(1,p[k])*se2[k]))
      grandSf <- c(grandSf,t(diag(Sf[[k]])))
      grandV_temp <- pracma::blkdiag(grandV_temp,V_ind[[k]])
      grandB_temp <- cbind(grandB_temp,B[[k]])
      Delta_temp <- c(Delta_temp,rep(1,r[k])*(1/se2[k]))
    }
    grandse2_inv <- 1/grandse2 #1*sum(p)
    grandSf_inv <- 1/grandSf #1*(r0+sum(r))
    grandV <- cbind(V_joint, grandV_temp)#sum(p)*(r0+sum(r))
    grandB <- cbind(B0,grandB_temp) #q*(r0+sum(r))
    Delta1 <- Delta_temp #bsxfun("*", grandV',grandse2_inv)%*%grandV #1*(r0+sum(r))
    Delta2 <- grandSf_inv+Delta1 #1*(r0+sum(r))
    temp <- tcrossprod(pracma::arrayfun("*",grandV, (1/Delta2)),grandV) #sum(p)*sum(p)
    SigmaY_inv <- diag(grandse2_inv)-pracma::arrayfun("*",pracma::arrayfun("*",temp,grandse2_inv),t(grandse2_inv)) #sum(p)*sum(p), not diagonal because of common structure, diff from SupSVD
    VSigmaYinvV <- Delta1-Delta1*(1/Delta2)*Delta1 #1*(r0+sum(r))
    EU <- X%*%pracma::arrayfun("*",grandB,(rep(1,sum(r)+r0)-VSigmaYinvV*grandSf))+grandY%*%SigmaY_inv%*%pracma::arrayfun("*",grandV,grandSf) # conditional mean

    #Print convergence information
    if(niter<max_niter){
      cat("SIPCA_B converges after",pracma::num2str(niter),"iterations.")
    }else{
      cat("SIPCA_B NOT converge after",pracma::num2str(max_niter),"iterations!!! Final change in angle: ",pracma::num2str(diff))
    }
    return(list(B0=B0, B=B, V_joint=V_joint,V_ind=V_ind,se2=se2,Sf0=Sf0,Sf=Sf,EU=EU))
  }
}
