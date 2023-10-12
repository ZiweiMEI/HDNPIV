library(splines)
library(splines2)
library(glmnet)
library(MASS) 

normalize <- function(X, base = NULL){
  X <- as.matrix(X)
  if (is.null(base)){
    base = X
  }
  base = as.matrix(base)
  
  p = dim(X)[2]
    for (ip in 1:p){
      Xcol = X[,ip]
      base.col = base[,ip]
      if (var(Xcol) > 0){
        X[,ip] = (Xcol - min(base)) / (max(base) - min(base))
      }
    }
  return(X)
}

SoftThreshold <- function( x, lambda ) {
  #
  # Standard soft thresholding
  #
  if (x>lambda){
    return (x-lambda);}
  else {
    if (x< (-lambda)){
      return (x+lambda);}
    else {
      return (0); }
  }
}


InverseLinftyOneRow <- function ( sigma, i, mu, maxiter=50, threshold=1e-2 ) {
  p <- nrow(sigma);
  rho <- max(abs(sigma[i,-i])) / sigma[i,i];
  mu0 <- rho/(1+rho);
  beta <- rep(0,p);
  
  if (mu >= mu0){
    beta[i] <- (1-mu0)/sigma[i,i];
    returnlist <- list("optsol" = beta, "iter" = 0);
    return(returnlist);
  }
  
  diff.norm2 <- 1;
  last.norm2 <- 1;
  iter <- 1;
  iter.old <- 1;
  beta[i] <- (1-mu0)/sigma[i,i];
  beta.old <- beta;
  sigma.tilde <- sigma;
  diag(sigma.tilde) <- 0;
  vs <- -sigma.tilde%*%beta;
  
  while ((iter <= maxiter) && (diff.norm2 >= threshold*last.norm2)){
    
    for (j in 1:p){
      oldval <- beta[j];
      v <- vs[j];
      if (j==i)
        v <- v+1;
      beta[j] <- SoftThreshold(v,mu)/sigma[j,j];
      
      if (oldval != beta[j]){
        vs <- vs + (oldval-beta[j])*sigma.tilde[,j];
      }
    }
    
    iter <- iter + 1;
    if (iter==2*iter.old){
      d <- beta - beta.old;
      diff.norm2 <- sqrt(sum(d*d));
      last.norm2 <-sqrt(sum(beta*beta));
      iter.old <- iter;
      beta.old <- beta;
      if (iter>10)
        vs <- -sigma.tilde%*%beta;
    }
  }
  
  returnlist <- list("optsol" = beta, "iter" = iter)
  return(returnlist)
}


InverseLinfty <- function(sigma, n, mm, resol=1.5, mu=NULL, maxiter=50, threshold=1e-2, verbose = TRUE) {
  isgiven <- 1;
  if (is.null(mu)){
    isgiven <- 0;
  }
  
  p <- nrow(sigma);
  M <- matrix(0, mm, p);
  xperc = 0;
  xp = round(p/10);
  mus <- numeric(mm)
  for (i in 1:mm) {
    if ((i %% xp)==0){
      xperc = xperc+10;
      if (verbose) {
        print(paste(xperc,"% done",sep="")); }
    }
    if (isgiven==0){
      # mu <- (1/sqrt(n)) * qnorm(1-(0.1/(p^2)));
      mu <- 0.01
    }
    mu.stop <- 0;
    try.no <- 1;
    incr <- 0;
    while ((mu.stop != 1)&&(try.no<10)){
      last.beta <- beta
      output <- InverseLinftyOneRow(sigma, i, mu, maxiter=maxiter, threshold=threshold)
      beta <- output$optsol
      
      iter <- output$iter
      if (isgiven==1){
        mu.stop <- 1
      }
      else{
        if (try.no==1){
          if (iter == (maxiter+1)){
            incr <- 1;
            mu <- mu*resol;
          } else {
            incr <- 0;
            mu <- mu/resol;
          }
        }
        if (try.no > 1){
          if ((incr == 1)&&(iter == (maxiter+1))){
            mu <- mu*resol;
          }
          if ((incr == 1)&&(iter < (maxiter+1))){
            mu.stop <- 1;
          }
          if ((incr == 0)&&(iter < (maxiter+1))){
            mu <- mu/resol;
          }
          if ((incr == 0)&&(iter == (maxiter+1))){
            mu <- mu*resol;
            beta <- last.beta;
            mu.stop <- 1;
          }
        }
      }
      try.no <- try.no+1
    }
    mus[i] <- mu
     
    M[i,] <- beta;
  }
  return(list(Theta.hat = M, mus = mus))
}


GenerateSplineMat <- function(X,degree = 3, M = NULL, gamma = 13/84, intercept = FALSE){
  X = as.matrix(X)
  n = dim(X)[1]
  p = dim(X)[2]
  if (is.null(M)){
    M  = max(ceiling( (2*n)^(gamma)))  
  }
  SplineMat = NULL 
  SplineMat.deriv = NULL  
  Spline2.list = list()
  Spline2deriv.list = list()
  num.list = 1
  for (ii in 1:p){
    X.ii = (X[,ii])
    len = (max(X[,ii])-min(X[,ii]))
    h = len/M
    
    Spline.mat <- bs(X.ii, degree = degree, df = M, intercept = intercept)
    derivative.mat <- dbs(X.ii, degree = degree, df = M, intercept = intercept)
    SplineMat = cbind(SplineMat, Spline.mat ) 
    SplineMat.deriv = cbind(SplineMat.deriv,  derivative.mat )
    Spline2.list[[num.list]] <- Spline.mat
    Spline2deriv.list[[num.list]] <- derivative.mat
    num.list = num.list + 1
    
  }
  
  return(list(M = M, SplineMat = SplineMat, SplineMat.deriv = SplineMat.deriv,
              Spline2.list = Spline2.list, Spline2deriv.list = Spline2deriv.list ))
}




 





EstimateEqD <- function(D,Z,X,search_lambda = Lasso_tuning, intercept = TRUE, seed = 2023, M = NULL, Mmin = M, Mmax = M){
  old <- .Random.seed
  on.exit( { .Random.seed <<- old } )
  set.seed(seed)
  D = as.matrix(D)
  Z = as.matrix(Z)
  X = as.matrix(X)
  n = dim(D)[1]
  pz = dim(Z)[2]
  px = dim(X)[2]
  
  if (is.null(M)){
    K = GenerateSplineMat(Z[,1], M = M)$SplineMat
    M = dim(K)[2]
  }
  
  fold.n <- floor(n/10)
  fold.id <- c(rep(1:10,floor(n/10)) , (1:10)[1:(n-10*floor(n/10))] )
  
  fold.id <- sample(fold.id, n)
  
  cvm.vec <- rep(NA,M - Mmin + 1)
  
  for (mm in Mmin:(M)){
    
    K = NULL 
    for (ell in 1:pz){
      K = cbind(K, GenerateSplineMat(Z[,ell], M = mm)$SplineMat) 
    }
    
    out.Lasso = cv.glmnet(cbind(X,K),D,intercept = intercept, foldid = fold.id, penalty.factor = c(rep(1,px),rep(0,dim(K)[2])))
    cvm.vec[mm - Mmin + 1] <- min(out.Lasso$cvm)
  }
  
  M = which.min(cvm.vec) + Mmin - 1
  K = NULL 
  for (ell in 1:pz){
    K = cbind(K, GenerateSplineMat(Z[,ell], M = M)$SplineMat) 
  }
  
  out.Lasso = cv.glmnet(cbind(X,K),D,intercept = intercept, foldid = fold.id, penalty.factor = c(rep(1,px),rep(0,dim(K)[2])))
  
  if (is.null(search_lambda)){
    Lasso.coef = as.matrix(as.vector(coef(out.Lasso, s = out.Lasso$lambda.min)),ncol = 1)
  }else if (search_lambda == "1se"){
    Lasso.coef = as.matrix(as.vector(coef(out.Lasso, s = out.Lasso$lambda.1se)),ncol = 1)
  }

  
 
  intercept.hat = Lasso.coef[1,]
  phi.hat = Lasso.coef[2:(px+1),]
  
  kappa.hat = Lasso.coef[-(1:(px+1)),]
  
  residual = D - intercept.hat - K %*% kappa.hat - X %*% phi.hat
  # print(residual)
  return( list(M = M, v.hat = residual, kappa.hat = kappa.hat, phi.hat = phi.hat, Spline.Z.Mat = K, intercept = intercept.hat ) )
}



EstimateEqY <- function(Y,D,Z,X,v.hat = NULL, H = NULL, degree = 3, search_lambda = Lasso_tuning, intercept = TRUE, 
                        seed = 2023, M = NULL, Mmin = M, Mmax = M){
  old <- .Random.seed
  on.exit( { .Random.seed <<- old } )
  set.seed(seed)
  # Y = as.matrix(Y)
  D = as.matrix(D)
  Z = as.matrix(Z)
  X = as.matrix(X)
  n = dim(D)[1]
  pz = dim(Z)[2]
  px = dim(X)[2]
  
  if (is.null(M)){
    Mat.D.list = GenerateSplineMat(D, M = M) 
    B = Mat.D.list$SplineMat
    M = dim(B)[2]
  }
  
  # print("M")
  # print(M)
  # 
  fold.n <- floor(n/10)
  fold.id <- c(rep(1:10,floor(n/10)) , (1:10)[1:(n-10*floor(n/10))] )
  
  fold.id <- sample(fold.id, n)
  
  cvm.vec <- rep(NA, Mmax - Mmin + 1)
  
  for (mm in Mmin:Mmax){
    
    if (is.null(v.hat)){
      out.EqD = EstimateEqD(D,Z,X, M = M,search_lambda = search_lambda)
      v.hat <- out.EqD$v.hat 
    }
    
    v.hat = matrix(v.hat) 
    
    Mat.D.list = GenerateSplineMat(D, M = mm) 
    Mat.v.list = GenerateSplineMat(v.hat, M = mm)
    B = Mat.D.list$SplineMat
    Bprime = Mat.D.list$Spline2deriv.list[[1]]
    if (is.null(H)){
      H = Mat.v.list$SplineMat
    }
    
    W = cbind(B,H) 
    
    out.Lasso = cv.glmnet(cbind(X,W), Y , intercept = intercept, foldid =  fold.id, penalty.factor = c(rep(1,px),rep(0,dim(W)[2])))
    
    cvm.vec[mm - Mmin + 1] <- min(out.Lasso$cvm)
  }
  
  M = which.min(cvm.vec) + Mmin - 1
   
  if (is.null(v.hat)){
    out.EqD = EstimateEqD(D,Z,X, M = M)
    v.hat <- out.EqD$v.hat 
  }
  
  v.hat = matrix(v.hat) 
  
  Mat.D.list = GenerateSplineMat(D, M = M) 
  Mat.v.list = GenerateSplineMat(v.hat, M = M)
  B = Mat.D.list$SplineMat
  # Bprime = Mat.D.list$Spline2deriv.list[[1]]
  if (is.null(H)){
    H = Mat.v.list$SplineMat
  }
  
  W = cbind(B,H) 
  
  out.Lasso = cv.glmnet(cbind(X,W), Y , intercept = intercept, foldid =  fold.id, penalty.factor = c(rep(1,px),rep(0,dim(W)[2])))
  
  if (is.null(search_lambda)){
    Lasso.coef = as.matrix(as.vector(coef(out.Lasso, s = out.Lasso$lambda.min)),ncol = 1)
  }else if (search_lambda == "1se"){
    Lasso.coef = as.matrix(as.vector(coef(out.Lasso, s = out.Lasso$lambda.1se)),ncol = 1)
  }
  
  intercept.hat = Lasso.coef[1,]
  
  theta.hat = Lasso.coef[2:(px+1),]
  
  
  omega.hat = Lasso.coef[-(1:(px+1)),]
  
  beta.hat = omega.hat[1:dim(B)[2]]
  eta.hat = omega.hat[-(1:dim(B)[2])]
  
  residual = Y - intercept.hat - X %*% theta.hat - H %*% eta.hat - B %*% beta.hat 
  
  return( list(M = M, B = B, H = H, Bprime = Bprime, beta.hat = beta.hat, eta.hat = eta.hat, theta.hat = theta.hat, v.hat = v.hat,
               Mat.D.list = Mat.D.list, Mat.v.list = Mat.v.list, residual = residual, intercept = intercept.hat ) )
}

