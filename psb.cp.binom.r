psb.cp.binom <-  function(nGa,kGa,X,R,Rz,draws=5000,burn=1000,Thin=5){
  #v2 only models variance in sampled mode
 # Thin=5
  K=30; #cut-point for stick-breaking process
  N=rep(0,K)
  Pi = rep(0,K)
  pp = rbeta(K,1,1)
  Probs = array(dim=c(Mi,Mt,K))
  Mi <- dim(nGa)[1]
  Mt <- 4
  D = dim(X)[2]
  lambda=0.001
  
  #initialize tensor bases
  U = list()
  U[[1]] = matrix(rnorm(D*R),ncol=R)
  VV = list()
  Q=c(4,K-1)
  MM=2
  for(mm in 1:MM) VV[[mm]] = matrix(rnorm(Q[mm]*R),ncol=R)
  Vmat = matrix(nrow=prod(Q),ncol=R)
  for(r in 1:R){
    Vr = lapply(VV, function(x) x[,r])
    Vmat[,r] = as.vector(array(apply(expand.grid(Vr), 1, prod), dim=Q))
  }
  
  alpha <- 0
  save_counter <- 1
  V <- array(dim=c(Mi,Mt,K))
  Z = array(rnorm(M*(K-1)), dim=c(Mi,Mt,K-1))
  Z_B=Z
  Int = Z
  D = dim(X)[2]
  B=array(dim=c(D,4,K-1))
  C=array(dim=c(Mi,Mt))
  
  #random effects
  z = array(dim=c(Mi,Mt,K-1))
  ZZ <-list()
  for(Dim in 1:3) ZZ[[Dim]] = matrix(rnorm(dim(z)[Dim]*Rz),ncol=Rz)
  ZZarray <- array(rep(0,prod(dim(z))),dim=dim(z))
  sigmaZZ <- rep(1,Rz)
  
  Probs = array(dim=c(Mi,Mt,K))
  B=array(dim=c(D,4,K-1))
  C=array(dim=c(Mi,Mt))
  pp = rbeta(K,1,1)
  
  numsave <-floor(draws/Thin)
  Vpreds <- array(dim=c(D,R,numsave))
  Vregions <- array(dim=c(Mt,R,numsave))
  Vprobs = array(dim=c(K-1,R,numsave))
  ppVec = array(dim=c(K,numsave))
  IntVec=array(dim=c(Mt,K-1,numsave))
  VregionsZ = array(dim=c(Mt,R,numsave))
  VprobsZ = array(dim=c(K-1,R,numsave))
  sdZZ <- rep(0,numsave)
  sdZB <- rep(0,numsave)
  sigmaVec <- rep(0,numsave)
  Iters=draws+burn
  for(jj in 1:Iters){
    ##Update weights using stick-breaking construction
    V[,,1:(K-1)] <- pnorm(Z)	
    V[,,K]= 1
    Pi=V
    Pi[,,2] =  V[,,2]*(1-V[,,1])
    for(k in 3:K){ 
      Pi[,,k] = V[,,k]*apply(1-V[,,1:(k-1)],FUN=prod,MARGIN=c(1,2))
    }
    ####Find data probs; estimate cluster allocations C
    for(k in 1:K){
      Probs[,,k] = Pi[,,k]*dbinom(kGa,nGa,pp[k])
    }
    for(i in 1:Mi){ for(j in 1:Mt){
      C[i,j] = c(1:K)[rmultinom(1,1,Probs[i,j,])==1]
    }}
    for(k in 1:K) pp[k] = rbeta(1,1+sum(kGa[C==k]), 1+sum(nGa[C==k])-sum(kGa[C==k]))

    z = array(dim=c(Mi,Mt,K-1))
    for(m in 1:Mi){for(j in 1:Mt){
      if(C[m,j]>1) z[m,j,1:(C[m,j]-1)] = rtnorm(C[m,j]-1,mean=Z[m,j,1:(C[m,j]-1)],upper=0) 
      if(C[m,j] < K) {z[m,j,C[m,j]] =rtnorm(1,mean=Z[m,j,C[m,j]],lower=0)}
      if(C[m,j] < K-1) {z[m,j,((C[m,j]+1):(K-1))] =rnorm(K-1-C[m,j],mean=Z[m,j,(C[m,j]+1):(K-1)])}
    }}
    
    ##Update Ints
    z_diff <- z-Z_B-ZZarray
    for(k in 1:(K-1)){ for(j in 1:Mt){
      Int[,j,k] = rep(rnorm(1,(alpha+sum(z_diff[,j,k]))/(1+Mi),sqrt(1/(1+Mi))),Mi)}}
    alpha <- rnorm(1,sum(Int[1,,])/((K-1)*Mt+1),sqrt(1/((K-1)*Mt+1)))
    
    ###update B
    Yvec = as.vector(z-ZZarray-Int)
    Ymats = list()
    for(mm in 1:MM){
      perm = c(mm+1,c(1:(MM+1))[-(mm+1)])
      Y_perm = aperm(z-ZZarray-Int,perm) 
      Ymats[[mm]] = array(Y_perm,dim=c(Q[mm],Mi*prod(Q[-mm])))
    }
    ##update U
    X_red = matrix(nrow=Mi,ncol=D*R)
    for(i in 1:Mi){X_red[i,] = as.vector(rep(X[i,],R))} 
    CC = matrix(nrow=dim(X_red)[1]*prod(Q),ncol=dim(X_red)[2])
    for(r in 1:R){
      index = ((r-1)*D+1):(r*D)
      CC[,index] = kronecker(Vmat[,r],X_red[,index])
    }
    lambdaMat = kronecker(t(Vmat)%*%Vmat,lambda*diag(D))
    CP <- crossprod(X_red)*kronecker(t(Vmat)%*%Vmat, matrix(rep(1,D^2),nrow=D))
    regMat = solve(CP+lambdaMat) 
    UlvecMean = regMat%*%(t(CC)%*%Yvec)
    UlvecVar = regMat
    Ulvec = mvrnorm(1,UlvecMean,UlvecVar)
    #standardize:
    #     Ulvec <- Ulvec/norm(Ulvec,type='F')
    U = matrix(Ulvec,nrow=D)
    
    for(mm in 1:MM){
      Bmat = U
      BRmat = matrix(nrow = prod(c(D,Q[-mm])),ncol=R)
      DDD = matrix(nrow=Mi*prod(Q[-mm]),ncol=R)
      for(r in 1:R){
        Vecs <- lapply(c(list(U),VV[-mm]), function(x) x[,r])
        Br = array(apply(expand.grid(Vecs), 1, prod), dim=c(D,prod(Q[-mm])))
        BRmat[,r] = as.vector(Br)
        DDD[,r] = as.vector(X %*% Br)
      }
      regMat = solve(t(DDD)%*%DDD+lambda*t(BRmat)%*%BRmat)
      VV[[mm]] = Ymats[[mm]]%*%DDD%*%regMat 
      VV[[mm]] = VV[[mm]] + mvrnorm(Q[mm],rep(0,R),regMat)
    #  VV[[mm]] = scale(VV[[mm]],center=FALSE)
      for(r in 1:R){
        sc <- sqrt(mean(VV[[mm]][,r]^2))
        VV[[mm]][,r] = VV[[mm]][,r]/sc
        U[,r] = U[,r]*sc
        if(mm==1 && mean(VV[[mm]][,r])<0){
          VV[[mm]][,r]<- -VV[[mm]][,r]
          U[,r] <- -U[,r] 
        }
        Vr = lapply(VV, function(x) x[,r])
        Vmat[,r] = as.vector(array(apply(expand.grid(Vr), 1, prod), dim=Q))
      }
    }  
    B = array(Bmat%*%t(Vmat),dim = c(D,Q))
    for(j in 1:4){for(k in 1:(K-1)){Z_B[,j,k] = t(X%*%B[1:D,j,k])}} 
    
    ##Now update indiv effects ZZ   
    zmats = list()
    for(Dim in 1:3){
      perm = c(Dim,c(1:3)[-(Dim)])
      z_perm = aperm(z-Z_B-Int,perm) 
      zmats[[Dim]] = array(z_perm,dim=c(dim(z)[Dim],prod(dim(z)[-Dim])))
      BRmat = matrix(nrow = prod(dim(z)[-Dim]),ncol=Rz)
      DDD = matrix(nrow=prod(dim(z)[-Dim]),ncol=Rz)
      for(r in 1:Rz){
        Vecs <- lapply(c(ZZ[-Dim]), function(x) x[,r])
        DDD[,r] = apply(expand.grid(Vecs), 1, prod)
      }
      if(Dim==1&Rz>1) Reg = diag(1/sigmaZZ)
      if(Dim>1&Rz>1) Reg = diag(Rz)
      if(Dim==1&Rz==1) Reg = 1/sigmaZZ
      if(Dim>1&Rz==1) Reg = 1
      regMat = solve(t(DDD)%*%DDD+Reg)
      ZZ[[Dim]] = zmats[[Dim]]%*%DDD%*%regMat 
      ZZ[[Dim]] = ZZ[[Dim]] + mvrnorm(dim(z)[Dim],rep(0,Rz),regMat)
    }
    #fix to add at each R
    for(r in 1:Rz){
      Zr = lapply(ZZ, function(x) x[,r])
      ZZarray = array(apply(expand.grid(Zr), 1, prod), dim=dim(z))
    }
    Z=Z_B+ZZarray+Int
    
  #  IG[1,1] prior
    sigmaZZ <- 1/rgamma(Rz,4+Mi/2,0.2+(Mi/2)*colMeans(ZZ[[1]]^2))
    if(jj>burn & jj%%Thin==0){
      Vpreds[,,save_counter] = U
      Vregions[,,save_counter] = VV[[1]]
      Vprobs[,,save_counter] = VV[[2]]
      ppVec[,save_counter] = pp
      IntVec[,,save_counter]=Int[1,,]
      VregionsZ[,,save_counter] = ZZ[[2]]
      VprobsZ[,,save_counter] = ZZ[[3]]
      sigmaVec[save_counter] = sigmaZZ
      sdZZ[save_counter] <- sd(ZZarray[,,1:10])
      sdZB[save_counter] <- sd(Z_B[,,1:10])
      save_counter <- save_counter+1
    }
#    print(c(sd(ZZarray[,,1:10]),sd(Z_B[,,1:10])))
 #   print(sigmaZZ)
#    print(jj)
  }
  
  return(list(Vpreds=Vpreds,Vregions=Vregions,Vprobs=Vprobs,ppVec=ppVec,IntVec=IntVec,VregionsZ=VregionsZ,VprobsZ=VprobsZ,sigmaVec=sigmaVec,sdZZ=sdZZ,sdZB=sdZB))
}
