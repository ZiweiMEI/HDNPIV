library(Rmosek)
library(Jmisc)
library(matrixStats)

find.lambda <- function(B,H,X.mat,K.hat,B0,pz){
  
 
  B0 = as.matrix(B0)
  
  n = dim(X.mat)[1]
  px = dim(X.mat)[2]
  MB = dim(B)[2]
  MH = dim(H)[2]
  MK = dim(K.hat)[2]
  
  B0.L2 = sqrt(sum(B0^2))
  
  F.hat <- cbind(B,H,X.mat,K.hat)
  pF = ncol(F.hat)
  
  prob <- list(sense="min")
  ### specify sup_norm constraints
  prob$c <-c(rep(0,pF),1) 
  
  zero.temp = matrix(0,nrow = 2*px, ncol = 1+1+pz)
  
  AAA <- n^(-1) * t(F.hat) %*% F.hat

  prob$bc <- rbind(blc=c( c(B0,rep(0,0)),rep(-Inf,pF+0)),
                  buc=c(rep(Inf,pF+0),c(B0,rep(0,0)) ))
  
  prob$A <- Matrix( rbind(cbind(AAA, 1),
                          cbind(AAA, -1)), sparse=TRUE)
  
  
  prob$bx <- rbind(blx=c(rep(-Inf,pF), 0),
                   bux=c(rep(Inf,pF), Inf ))
  
  
  # Return the solution
  sol <- try(mosek(prob, list(verbose = 0)), silent=TRUE)
  solution.vec <- sol$sol$itr$xx 
  
  opt.obj <- solution.vec[pF+1] 
  lambda.found =  max(opt.obj)  / B0.L2 
  return(lambda.found)
}


proj.debias <- function(B,H,X.mat,K.hat,B0,pz,lambda = NULL,lambda.X = NULL,kappa0 = 1.2){
  
  B0 = as.matrix(B0)
  n = dim(X.mat)[1]
  px = dim(X.mat)[2]
  MB = dim(B)[2]
  MH = dim(H)[2]
  MK = dim(K.hat)[2]
  
  if (is.null(lambda)){
    # if (n <= px){
    #   lambda = find.lambda(B,H,X.mat,K.hat,B0,pz)
    # }else{
    #   n.tuning = sample(1:n,px)
    #   lambda = (px/n)^(3/14) * find.lambda(B[n.tuning,],H[n.tuning,],X.mat[n.tuning,],K.hat[n.tuning,],B0,pz)
    # }
    lambda = find.lambda(B,H,X.mat,K.hat,B0,pz)
    print(c("lambda",lambda))
  }
  
  
  B0.L2 = sqrt(sum(B0^2))
  
  F.hat <- cbind(B,H,X.mat,K.hat)
 
  pF = ncol(F.hat)
  xi0 = c(B0, rep(0,pF-length(B0)))
  
  optimal = FALSE
  while(!optimal){
    # print(optimal)
    if (lambda <= 0){
      lambda = 0.01
      } 
    lambda = kappa0 * lambda 
    prob <- list(sense="min")
   
   ### specify sup_norm constraints
    prob$c <-c(rep(0,pF),1) 
    # zero.temp = matrix(0,nrow = px, ncol = 1+1)
    prob$A <- Matrix( cbind(rbind(  n^(-1)*t(F.hat)%*%F.hat,   n^(-1/2) * F.hat ),0), sparse=TRUE)
    prob$bc <-rbind(blc=c(rep(-lambda,pF) + B0, rep(-MB^(1.5) * sqrt(log(pF)) * kappa0 , n )), 
                     buc=c(rep(lambda,pF) + B0, rep(MB^(1.5) * sqrt(log(pF)) * kappa0 , n)) )
    prob$bx <- rbind(blx=c(rep(-Inf,pF),0),
                     bux=c(rep(Inf,pF),Inf))
    
    
    prob$F <- Matrix( rbind(
      c(rep(0,pF),1),
      cbind(n^(-1/2)*F.hat, matrix(0,nrow = n, ncol = 1) )
    ), sparse = TRUE)
    prob$g <-  rep(0,nrow( prob$F ))
    # Solve the problem
    prob$cones <- matrix(list(), nrow=3, ncol=1)
    rownames(prob$cones) <- c("type","dim","conepar")
    
    prob$cones[,1] <- list("QUAD", n+1, NULL)
    # sol <- mosek(prob)
    sol <- try(mosek(prob, list(verbose=0)), silent=TRUE)

    optimal <- identical(sol$response$code, 0)
   }
  
  opt.sol <- sol$sol$itr$xx[1:pF]
  # print(max(abs(opt.sol)))
  opt.obj <- sol$sol$itr$xx[pF+1]
  
  return(list(m.hat = opt.sol, opt.obj = opt.obj, lambda = lambda))
 
}

Inference.gfun <- function(Y,D,Z,X,gamma = 13/84, M = NULL, split = FALSE, 
                           sample.1.index = NULL, d.seq = NULL, Lasso_tuning = NULL,
                           gfun = FALSE, sig.level = 0.05, NBB = 1e4, kappa0 = 1.2){
 
  D = as.matrix(D)
  Z = as.matrix(Z)
  X = as.matrix(X)
  n = dim(D)[1]
  pz = dim(Z)[2]
  px = dim(X)[2]
  
  if (!split){
    sample.1 = 1:n 
    sample.2 = 1:n 
  }
  
  if (split & is.null(sample.1.index)){
    sample.1 = 1:max(floor(n/2),100)
    sample.2 = (max(floor(n/2),100)+1):n 
  }
  
  if (split & !is.null(sample.1.index)){
    sample.1 = sample.1.index
    sample.2 = (1:n)[-sample.1.index]
  }
   
  gfun.init = rep(NA,length(d.seq))
  gderiv.init = rep(NA,length(d.seq))
  gfun.est = rep(NA,length(d.seq))
  gderiv.est = rep(NA,length(d.seq))
  gfun.sd = rep(NA,length(d.seq))
  gderiv.sd = rep(NA,length(d.seq))
  # gderiv.sd.uniform = rep(NA,length(d.seq))
  
  n.a = length(sample.1)
  n.b = length(sample.2)
  
  print(c("na",n.a))
  print(c("nb",n.b))
  
  if (!is.null(gamma)){
      print("M")
      
      M = max( ceiling((2*n)^gamma),5)  
      
      print(M)
      print(n)
  }
  
  # print(M)
  
  Y.a = as.matrix(Y[sample.1,])
  D.a = as.matrix(D[sample.1,])
  Z.a = as.matrix(Z[sample.1,])
  X.a = as.matrix(X[sample.1,])
  
  Y.b = as.matrix(Y[sample.2,])
  D.b = as.matrix(D[sample.2,])
  Z.b = as.matrix(Z[sample.2,])
  X.b = as.matrix(X[sample.2,])
  
  if (is.null(d.seq)){
    d.seq = quantile(D,seq(.05,.95,.05))
    }
  
   
  
  out.EqD.a = EstimateEqD(D.a,Z.a,X.a, M = M, search_lambda = Lasso_tuning)
  K.b = NULL 
  for (ell in 1:pz){
    # print(pz)
    # print(dim(Z.b))
    K.b = cbind(K.b, GenerateSplineMat(Z.b[,ell], M = length( out.EqD.a$kappa.hat)/pz )$SplineMat) 
  }
  
  
 
  v.hat.b = D.b - K.b %*% out.EqD.a$kappa.hat - X.b %*% out.EqD.a$phi.hat - out.EqD.a$intercept
 
  out.EqY.a = EstimateEqY(Y.a,D.a,Z.a,X.a, M = M, search_lambda = Lasso_tuning)
  # Hprime.a = out.EqY.a$Mat.v.list$SplineMat.deriv
  eta.hat.a =  out.EqY.a$eta.hat
  
  # K.b = out.EqD.b$Spline.Z.Mat
  out.EqY.b = EstimateEqY(Y.b,D.b,Z.b,X.b, M = M, search_lambda = Lasso_tuning, v.hat = v.hat.b )
  
  beta.hat.b <- out.EqY.b$beta.hat
   
  
  out.EqY = EstimateEqY(Y,D,Z,X, M = M, search_lambda = Lasso_tuning)
  
  beta.hat <- out.EqY$beta.hat
  eta.hat <- out.EqY$eta.hat
  theta.hat <- out.EqY$theta.hat
  
  B.b = out.EqY.b$B
  Bprime.b = GenerateSplineMat(D.b, M = dim(B.b)[2] )$Spline2deriv.list[[1]]
  Bprime = GenerateSplineMat(D, M = length(beta.hat) )$Spline2deriv.list[[1]]
  B.b = GenerateSplineMat(D.b, M = dim(B.b)[2] )$Spline2.list[[1]]
   
  H.b =  out.EqY.b$H
  Hprime.b = GenerateSplineMat(v.hat.b, M = length(eta.hat.a))$Spline2deriv.list[[1]]
  
  
  q.deriv.hat =  Hprime.b %*% eta.hat.a   
  
  X.hat.b = X.b * as.vector(q.deriv.hat)
  X.mat.b <- cbind(X.b,X.hat.b)
  K.hat.b = K.b * as.vector(q.deriv.hat)
   
  # print(B.b)
  F.hat.b = cbind(B.b,H.b,X.mat.b,K.hat.b)
  pF = dim(F.hat.b)[2]
  
  Sigma.hat = t( demean(F.hat.b) ) %*% demean(F.hat.b) / n.b  
 

  
  e.hat <- out.EqY$residual 
  var.hat = mean(e.hat^2)
  
  ehat.b = demean(e.hat[sample.2] + out.EqY$intercept)
  
  print(c(mean( e.hat ),out.EqY$intercept,var.hat))
  
  BB = GenerateSplineMat(D, M = length(beta.hat) )$Spline2.list[[1]]
  Bd0 = t(predict(BB,d.seq)) 
  gfun.init <- t(Bd0) %*% beta.hat
 
  H.bootstrap <- matrix(NA,BB,length(d.seq))
  
  proj.mat = matrix(0,pF,length(beta.hat))
  for (jj in 1:length(beta.hat)){ 
    ejj = rep(0,pF)
    ejj[jj] = 1
    
    debias.fit <- proj.debias( demean(B.b), demean(H.b), demean(X.mat.b), demean(K.hat.b), B0 = ejj, pz, kappa0 = kappa0)
     
    Omega.jj <- debias.fit$m.hat
    
    obj <- debias.fit$opt.obj
    # print("obj")
    # print(c(obj^2,t(Omega.jj) %*% Sigma.hat %*% Omega.jj, max(abs(n.b^(-1/2) * M^(-1.5) * F.hat.b%*%Omega.jj)) ))
    
    proj.mat[,jj] <- c(Omega.jj)
  }
 
  
   

  Bprime.d0 = t(predict(Bprime,d.seq))
 
  m.hat <- proj.mat  %*%  Bprime.d0
  
  gderiv.init <- t(Bprime.d0) %*% beta.hat
  
  gderiv.est <- gderiv.init + 1/n.b * t(m.hat) %*% t(demean(F.hat.b)) %*% ehat.b
  gderiv.sd = sqrt(diag(var.hat * t(m.hat) %*% Sigma.hat %*% m.hat / n.b))
  
  gfun.est <- gfun.init + 1/n.b * t(Bd0) %*% t(proj.mat) %*% t(demean(F.hat.b)) %*% ehat.b 
  
  H.bootstrap <- 1/sqrt(n.b) * t(m.hat) %*% t(demean(F.hat.b)) %*% matrix(rnorm(n.b*NBB ), nrow = n.b) 
  
  H.bootstrap = H.bootstrap * sqrt( diag(1 / (t(m.hat) %*% Sigma.hat %*% m.hat )))
   
  H.quantile <- quantile(colMaxs(abs(H.bootstrap)), 1-sig.level)
  
  return(list(gfun.init=gfun.init,gfun.est=gfun.est,gderiv.init=gderiv.init,gderiv.est=gderiv.est, gderiv.sd=gderiv.sd, 
                d.seq = d.seq,
                H.quantile = H.quantile,
                Lasso_tuning = Lasso_tuning)) 
  
  
}