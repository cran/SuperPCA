#' Supervised Sparse and Functional Principal Component Analysis
#'
#' This function conducts supervised sparse and functional principal
#' component analysis by fitting the SupSVD model
#'       X=UV' + E
#'       U=YB + F
#' where X is an observed primary data matrix (to be decomposed), U is a latent score
#' matrix, V is a loading matrix, E is measurement noise, Y is an observed
#' auxiliary supervision matrix, B is a coefficient matrix, and F is a
#' random effect matrix. \cr
#' \cr
#' It decomposes the primary data matrix X into low-rank
#' components, while taking into account many different features: 1)
#' potential supervision from any auxiliary data Y measured on the same
#' samples; 2) potential smoothness for loading vectors V (for functional
#' data); 3) sparsity in supervision coefficients B and loadings V (for variable
#' selection). \cr
#' \cr
#' It is a very general dimension reduction method that subsumes
#' PCA, sparse PCA, functional PCA, supervised PCA, etc as special cases.
#' See more details in 2016 JCGS paper "Supervised sparse and
#' functional principal component analysis" by Gen Li, Haipeng Shen, and
#' Jianhua Z. Huang.
#'
#' @param Y n*q (column centered) auxiliary data matrix, rows are samples and columns are variables
#' @param X n*p (column centered) primary data matrix, which we want to decompose.
#'          rows are samples (matched with Y) and columns are variables
#' @param r positive scalar, prespecified rank (r should be smaller than n and p)
#' @param ind_lam 0 or 1 (default=1, sparse loading), sparsity index for loadings
#' @param ind_alp 0 or 1 (default=1, smooth loading), smoothness index for loadings
#' @param ind_gam 0 or 1 (default=1, sparse coefficient), sparsity index
#'                   for supervision coefficients. Note: if gamma is set to be 0,
#'                   Y must have q<n to avoid overfitting; if gamma is set to be 1,
#'                   then it can handle high dimensional supervision Y
#' @param ind_Omg p*p symmetric positive semi-definite matrix for
#'                   smoothness penalty (default is for evenly spaced data)
#'                   Note: only change this if you have unevenly spaced
#'                   functional data X
#' @param Omega ??
#' @param max_niter scalar (default=1E3), max number of overall iterations
#' @param convg_thres positive scalar (default=1E-6), overall convergence threshold
#' @param vmax_niter scalar (default=1E2), max number of iterations for estimating
#'                   each loading vector
#' @param vconvg_thres positive scalar (default=1E-4), convergence
#'                      threshold for the proximal gradient descent algorithm
#'                      for estimating each loading vector
#'
#' @return list with components
#' \item{B:}{q*r coefficient matrix of Y on the scores of X,maybe sparse if gamma=1}
#' \item{V:}{p*r loading matrix of X, each column has norm 1, but no strict orthogonality
#'           because of sparsity and smoothness. If lambda=1, V is sparse;
#'           if alpha=1, each column of V is smooth}
#' \item{U:}{n*r score matrix of X, conditional expectation of random scores,
#'           no strict orthogonality}
#'
#' \item{se2:}{scalar, variance of measurement error in the primary data X}
#' \item{Sf:}{r*r diagonal covariance matrix, for random effects (see paper)}
#' Note: Essentially, U and V are the most important output for dimension
#' reduction purpose as in PCA or SVD.
#' @export
#'
#' @examples
#' \dontrun{
#' library(spls)
#' data(yeast)
#' r <- 4
#' ydata <- as.data.frame(yeast[1])
#' xdata <- as.data.frame(yeast[2])
#' yc <- scale(ydata,center = TRUE,scale=FALSE)
#' xc <- scale(xdata,center=TRUE,scale=FALSE)
#' SupSFPCA(yc,xc,r)
#' }
SupSFPCA <- function(Y,X,r,ind_lam = 1,
                     ind_alp = 1,
                     ind_gam = 1,
                     ind_Omg = 1,
                     Omega = 0,
                     max_niter = 10^3,
                     convg_thres = 10^-6,
                     vmax_niter = 10^2,
                     vconvg_thres = 10^-4){

  # pre-check
  if(max(abs(colMeans(X)))>0.01 || max(abs(colMeans(Y)))>0.01){
    stop("Columns of X and Y are not centered. exit...")}


  n <- nrow(X)
  p <- ncol(X)
  n1 <- nrow(Y)
  q <- ncol(Y)
  if(qr(Y)$rank < q){
    stop("Do not run this code! Change initial and BIC df...")
    stop("gamma cannot be set to zero!")
  }

  # Pre-check
  if (n!=n1){
    stop("X does not match Y! exit...")
  }else if(qr(Y)$rank != q){
    stop("Columns of Y are linearly dependent! exit...")
  }else if(r>n || r>p){
    stop("Too greedy on ranks! exit...")
  }


  # set Omega
  if (ind_Omg == 1){
    Q1 <- Matrix::sparseMatrix(1:p, 1:p, x=-2, dims = c(p,p))
    Q2 <- Matrix::sparseMatrix(1:(p-1), 2:p, x = 1, dims = c(p,p))
    Q3 <- Matrix::sparseMatrix(2:p, 1:(p-1), x = 1, dims = c(p,p))
    Q <- Q1+Q2+Q3
    Q <- Q[,2:(p-1)]
    R1 <- Matrix::sparseMatrix(1:(p-2), 1:(p-2), x=2/3, dims = c(p-2,p-2))
    R2 <- Matrix::sparseMatrix(1:(p-1-2), 2:(p-2), x = 1/6, dims = c(p-2,p-2))
    R3 <- Matrix::sparseMatrix(2:(p-2), 1:(p-1-2), x = 1/6, dims = c(p-2,p-2))
    R <- R1+R2+R3
    Omega <- Q%*%solve(R)%*%t(Q) #p*p
  }
  Omega <- as.matrix(Omega)
  oeig <- max(eigen(Omega)$value)

  if(r == 1){
    # initial est
    S <- RSpectra::svds(X,r)
    U <- S$u
    D <- S$d
    D <-diag(as.matrix(D))
    V <- S$v
    U <-U%*%D
    E <- X - tcrossprod(U,V)
    se2 <- stats::var(as.vector(E))
    B <- tcrossprod(solve(crossprod(Y,Y)),Y)%*%U
    Sf <- diag(as.matrix(diag((1/n)*crossprod(U-Y%*%B, U-Y%*%B)))) #initial Sf clear E D

    diff <- 1 # EM criterion
    niter <- 0 # number of EM iter
    #max_niter <- 18
    while (niter <= max_niter && diff > convg_thres){ # record last iter
      if (pracma::mod(niter,2)==0){
        cat(sprintf('This is the %.d th iterations, the diff= %.4g \n',niter,diff))
      }

      se2_old <- se2
      Sf_old <- Sf
      V_old <- V
      B_old <- B

      # E step
      # some critical values
      Sfinv <- solve(Sf)
      weight <- solve(diag(r) + se2*Sfinv) #r*r
      cond_Mean <- (se2*Y%*%B%*%Sfinv + X%*%V) %*% weight # E(U|X), n*r
      cond_Var = Sf%*%(diag(r) - weight) #cov(U(i)|X), r*r
      cond_quad = n*cond_Var + crossprod(cond_Mean,cond_Mean) # E(U'U|X), r*r

      # M step
      # estimate B and Sf
      if(ind_gam != 0){
        for(k in 1:r){
          #Attention: lasso fcn always center X,Y, and return a
          # separate column of intercept; by default, it standardize
          # each column of X
          # Therefore, when using lasso, we always center the columns of
          # X and Y to aviod the trouble of intercept.
          fit <- glmnet::glmnet(Y, cond_Mean[,k], alpha = 1, lambda.min.ratio = 0, standardize = FALSE)
          BIC_score <- n*stats::deviance(fit)+log(n)*fit$df
          ind <- which(BIC_score == min(BIC_score))
          #lambda11 <- fit$lambda[ind]
          B[,k] <- fit$beta[,ind]
        }
      }
      # if gamma = 0, no B-sparsity
      else{
        B <- tcrossprod(solve(crossprod(Y,Y)),Y)%*%cond_Mean
      }

      # estimate Sf
      Sf <- diag(as.matrix(diag((cond_quad + crossprod(Y%*%B, Y%*%B) - t(Y%*%B)%*%cond_Mean-t(cond_Mean)%*%(Y%*%B))/n))) #r*r

      # estimate V
      V <- NULL
      for (k in 1:r){ # kth layer, some critical values
        theta <- t(X)%*%cond_Mean[,k] - (V_old%*%cond_quad[,k]-V_old[,k]*cond_quad[k,k]) #p*1
        c <- cond_quad[k,k] # E(u_k'u_k|X), 1*1

        # select smooth tuning (LOOCV w/o refitting)
        if(ind_alp != 0){
          alphavec <- seq(0,10,0.1) # smooth tuning range
          cv_score <- rep(0,length(alphavec))
          for (ialp in 1:length(alphavec)){
            alpha <- alphavec[ialp]
            hat <- solve(diag(p)+alpha*Omega)
            vest <- hat%*%theta/c
            cv_score[ialp] <- (1/p)*colSums((((theta/c)-vest)/(1-diag(hat)))^2)
          }

          I <- which.min(cv_score)
          optalp <- alphavec[I] # Optimal alpha for this iteration
        }
        else{ # no V-smoothness
          optalp <- 0
        }

        # specify sparsity tuning (for gaussian error)
        if(ind_lam != 0){
          optlam <- sqrt(2*log(p)*se2_old/c)
        }
        else{ # no V sparsity
          optlam <- 0
        }

        L <- 1 + optalp*oeig
        vk_old <- V_old[,k] # initial value for v_k is from last EM iteration
        vdiff <- 1
        vniter <- 0
        while(vniter <= vmax_niter && vdiff > vconvg_thres){
          # iteration for estimating v_k
          df <- -theta/c + (diag(p) + optalp*Omega)%*%vk_old
          vk <- soft_thr(vk_old - (1/L)*df, optlam/L)

          # set norm = 1
          if(norm(vk,type = "F") == 0){
            warning('zero vector v!')
          }
          else{
            vk <- vk/norm(vk,type = "F")
          }
          vdiff <- norm(matrix(vk-vk_old),type = "F")^2
          vk_old <- vk
          vniter <- vniter+1
        }
        V <- cbind(V,matrix(vk))
      }

      # Estimate se2
      se2 <- (psych::tr(X%*%(t(X) - 2*tcrossprod(V, cond_Mean)))
              + n*psych::tr(crossprod(V,V)%*%cond_Var)
              + psych::tr(tcrossprod(tcrossprod(cond_Mean, V)%*%V, cond_Mean))) / (n*p)

      # stopping rule
      diff <- norm(V-V_old, type = "F")^2
      niter <- niter+1
    }

    # Print convergence information
    if (niter < max_niter){
      cat("SupFPCA converges at precision", pracma::num2str(convg_thres,fmt=10), "after", pracma::num2str(niter), "iterations.")
    }else{
      cat("SupFPCA NOT converge at precision", pracma::num2str(convg_thres,fmt=10), "after", pracma::num2str(max_niter), "iterations!!!")
    }

    # reorder V and others
    I <- order(diag(as.matrix(t(V)%*%t(X)%*%X%*%V)), decreasing = TRUE)
    V <- V[,I]
    B <- B[,I]
    Sf <- as.matrix(Sf)[I,I]

    # output U
    Shinv <- solve(Sf)
    weight <- solve(diag(r) + se2*Sfinv) # r*r
    U <- (se2*Y%*%B%*%Sfinv + X%*%V)%*%weight
  }

  else{
    # initial est
    S <- RSpectra::svds(X,r)
    U <- S$u
    D <- S$d
    D <-diag(D)
    V <- S$v
    U <-U%*%D
    E <- X - tcrossprod(U,V)
    se2 <- stats::var(as.vector(E))
    B <- tcrossprod(solve(crossprod(Y,Y)),Y)%*%U
    Sf <- diag(diag((1/n)*crossprod(U-Y%*%B, U-Y%*%B))) #initial Sf clear E D

    diff <- 1 # EM criterion
    niter <- 0 # number of EM iter
    #max_niter <- 18
    while (niter <= max_niter && diff > convg_thres){ # record last iter
      if (pracma::mod(niter,2)==0){
        cat(sprintf('This is the %.d th iterations, the diff= %.4g \n',niter,diff))
      }

      se2_old <- se2
      Sf_old <- Sf
      V_old <- V
      B_old <- B

      # E step
      # some critical values
      Sfinv <- solve(Sf)
      weight <- solve(diag(r) + se2*Sfinv) #r*r
      cond_Mean <- (se2*Y%*%B%*%Sfinv + X%*%V) %*% weight # E(U|X), n*r
      cond_Var = Sf%*%(diag(r) - weight) #cov(U(i)|X), r*r
      cond_quad = n*cond_Var + crossprod(cond_Mean,cond_Mean) # E(U'U|X), r*r

      # M step
      # estimate B and Sf
      if(ind_gam != 0){
        for(k in 1:r){
          #Attention: lasso fcn always center X,Y, and return a
          # separate column of intercept; by default, it standardize
          # each column of X
          # Therefore, when using lasso, we always center the columns of
          # X and Y to aviod the trouble of intercept.
          fit <- glmnet::glmnet(Y, cond_Mean[,k], alpha = 1, lambda.min.ratio = 0, standardize = FALSE)
          BIC_score <- n*stats::deviance(fit)+log(n)*fit$df
          ind <- which(BIC_score == min(BIC_score))
          #lambda11 <- fit$lambda[ind]
          B[,k] <- fit$beta[,ind]
        }
      }
      # if gamma = 0, no B-sparsity
      else{
        B <- tcrossprod(solve(crossprod(Y,Y)),Y)%*%cond_Mean
      }

      # estimate Sf
      Sf <- diag(diag((cond_quad + crossprod(Y%*%B, Y%*%B) - t(Y%*%B)%*%cond_Mean-t(cond_Mean)%*%(Y%*%B))/n)) #r*r

      # estimate V
      V <- NULL
      for (k in 1:r){ # kth layer, some critical values
        theta <- t(X)%*%cond_Mean[,k] - (V_old%*%cond_quad[,k]-V_old[,k]*cond_quad[k,k]) #p*1
        c <- cond_quad[k,k] # E(u_k'u_k|X), 1*1

        # select smooth tuning (LOOCV w/o refitting)
        if(ind_alp != 0){
          alphavec <- seq(0,10,0.1) # smooth tuning range
          cv_score <- rep(0,length(alphavec))
          for (ialp in 1:length(alphavec)){
            alpha <- alphavec[ialp]
            hat <- solve(diag(p)+alpha*Omega)
            vest <- hat%*%theta/c
            cv_score[ialp] <- (1/p)*colSums((((theta/c)-vest)/(1-diag(hat)))^2)
          }

          I <- which.min(cv_score)
          optalp <- alphavec[I] # Optimal alpha for this iteration
        }
        else{ # no V-smoothness
          optalp <- 0
        }

        # specify sparsity tuning (for gaussian error)
        if(ind_lam != 0){
          optlam <- sqrt(2*log(p)*se2_old/c)
        }
        else{ # no V sparsity
          optlam <- 0
        }

        L <- 1 + optalp*oeig
        vk_old <- V_old[,k] # initial value for v_k is from last EM iteration
        vdiff <- 1
        vniter <- 0
        while(vniter <= vmax_niter && vdiff > vconvg_thres){
          # iteration for estimating v_k
          df <- -theta/c + (diag(p) + optalp*Omega)%*%vk_old
          vk <- soft_thr(vk_old - (1/L)*df, optlam/L)

          # set norm = 1
          if(norm(vk,type = "F") == 0){
            warning('zero vector v!')
          }
          else{
            vk <- vk/norm(vk,type = "F")
          }
          vdiff <- norm(matrix(vk-vk_old),type = "F")^2
          vk_old <- vk
          vniter <- vniter+1
        }

        V <- cbind(V,matrix(vk))
      }


      # Estimate se2
      se2 <- (psych::tr(X%*%(t(X) - 2*tcrossprod(V, cond_Mean)))
              + n*psych::tr(crossprod(V,V)%*%cond_Var)
              + psych::tr(tcrossprod(tcrossprod(cond_Mean, V)%*%V, cond_Mean))) / (n*p)

      # stopping rule
      diff <- norm(V-V_old, type = "F")^2
      niter <- niter+1
    }

    # Print convergence information
    if (niter < max_niter){
      cat("SupFPCA converges at precision", pracma::num2str(convg_thres,fmt=10), "after", pracma::num2str(niter), "iterations.")
    }else{
      cat("SupFPCA NOT converge at precision", pracma::num2str(convg_thres,fmt=10), "after", pracma::num2str(max_niter), "iterations!!!")
    }

    # reorder V and others
    I <- order(diag(t(V)%*%t(X)%*%X%*%V), decreasing = TRUE)
    V <- V[,I]
    B <- B[,I]
    Sf <- Sf[I,I]

    # output U
    Shinv <- solve(Sf)
    weight <- solve(diag(r) + se2*Sfinv) # r*r
    U <- (se2*Y%*%B%*%Sfinv + X%*%V)%*%weight
  }
  return(list(U=as.matrix(U), V=as.matrix(V), B=as.matrix(B), Sf=as.matrix(Sf), se2=se2))
}

soft_thr<-function(In, lambda){
  out <- sign(In)*pmax(abs(In) - lambda, 0)
  return(out)
}
