# stop()

# Preamble ----------------------------------------------------------------

library(sna)
library(igraph)
library(intergraph)
library(Matrix)
library(rARPACK)
library(Rcpp)
library(RcppArmadillo)
library(Hmisc)
library(mvtnorm)
library(plyr)
library(dplyr)
library(Ternary)
library(numDeriv)
library(RColorBrewer)
library(dfoptim)

# sourceCpp("H:/edgeClustering/eClusteringCode/eClusteringCPPCode_v5f.cpp")

options(stringsAsFactors = FALSE)


# spectral clustering -----------------------------------------------------

### Helper function for spectral clustering
.eAdjacency = function(A){
  if(is_directed(A))A <- as.undirected(A,"collapse")
  n = vcount(A)
  M = ecount(A)
  allEnds = ends(A,1:M,names=F)
  
  # incMat = as(as.matrix(asNetwork(A),matrix.type="incidence"),"dgCMatrix")
  incMat = Matrix(0,n,M)
  incMat[cbind(c(t(allEnds)),
               rep(1:M,each=2))] = 1
  EMat = crossprod(incMat)
  diag(EMat)=0
  # EMat = Matrix(0,M,M)
  # for(i in 1:n){ #SLOW.....
  #   temp = which(allEnds[,1]==i | allEnds[,2]==i)
  #   EMat[temp,temp] = 1L
  # }
  # diag(EMat) = 0
  return(as(EMat,"dgCMatrix"))
}



spec2dir = function(A,eCl){
  A1 <- as.undirected(A,"collapse")
  A1 <- set_edge_attr(A1,"cluster",value=eCl)
  
  A0.Mat = as_adjacency_matrix(A,"both",sparse=TRUE)
  A1.Mat = as_adjacency_matrix(A1,"both",attr="cluster",sparse=TRUE)
  A0.Mat <- A0.Mat*A1.Mat
  
  ELOrig = as.data.frame(as_edgelist(A,F))
  EL = data.frame(Source = A0.Mat@i+1,
                  Target = rep(seq_along(diff(A0.Mat@p)),diff(A0.Mat@p)),
                  cluster = A0.Mat@x)
  colnames(ELOrig) = c("Source","Target")
  ELNew = dplyr::left_join(ELOrig,EL,by = c("Source","Target"))
  
  return(ELNew)
}


getMMSpec = function(A,eCl,weights=NULL,Jeffreys=FALSE){#weights = E(A)$weight
  if(is_directed(A))A <- as.undirected(A,"collapse")
  n = vcount(A)
  M = ecount(A)
  allEnds = ends(A,1:M,names=F)
  if(is.null(weights)){
    ve = cbind(c(allEnds),rep(eCl,2))
    tab = table(ve[,1],ve[,2])
    tab = tab[order(as.numeric(rownames(tab))),]
    if(Jeffreys) tab = tab + 0.5
    rs = rowSums(tab)
    tab = cbind(vegetarian::normalize.rows(tab),rs)
    return(tab)
  }else{
    ve = data.frame(node=c(allEnds),
                    cluster=rep(eCl,2),
                    weights=rep(weights,2))
    tab = aggregate(weights~node+cluster,data=ve,FUN=sum)
    ret = matrix(0,n,max(tab$cluster))
    ret[as.matrix(tab[,1:2])] = tab[,3]
    return(ret)
  }
}


# EM algo -----------------------------------------------------------------

#' Edge clustering
#' 
#' 'eClustGEM' implements the generalized EM algorithm for clustering edges of a network using a latent space approach
#' 
#' @param A igraph object
#' @param K integer.  Number of clusters
#' @param p integer. Dimension of latent space
#' @param maxIter integer. Maximum number of iterations of the GEM algorithm
#' @param eps numeric. Convergence threshold for relative change of log posterior
#' @param CGSteps integer.  Maximum number of conjugate gradient steps during each iteration
#' @param QNSteps integer.  Maximum number of quasi-Newton steps during each iteration
#' @param WApprox logical. Should subsampling approximation be used in updating W?
#' @param UVApprox logical. Should subsampling approximation be used in updating U and V?
#' @param UVEstMethod character, either 'CG' for conjugate gradient (default) or 'QN' for quasi-Newton
#' @param BKNInit should initialization be done via Ball Karrer, and Newman (2011)?
#' @param randUVInit Should the initialization of U and V be random?
#' @param a_u, b_u, a_s, b_s, a_v, b_v, a_r, b_r, a_0 numeric.  Hyperparameters. 
eClustGEM = function(A,K=2,p=2,maxIter=1e3,eps=1e-4,QNSteps=25,NMSteps=50,
                     HJSteps=5,CGSteps=25,
                     WApprox=ecount(A)>1e3,nSubsamp = max(vcount(A),500),
                     UVApprox=FALSE,
                     UVEstMethod = c("NM","QN","HJ","CG"),
                     SCInit=FALSE,BKNInit=TRUE,randUVInit = FALSE,
                     a_u = 1,b_u = 1,a_s=1, b_s=1,
                     a_v = a_u,b_v = b_u,a_r=a_s,b_r=b_s,a_0 = 1){
  ### Create objects
  M = ecount(A)
  n = vcount(A)
  EE = cbind(ends(A,1:ecount(A),FALSE))
  IVApprox = matrix(0,nSubsamp,K)
  
  cat("\n Indexing edgelist \n")
  temp = indexEdges(EE,n)
  Mi1 = temp$Mi1[,1:max(temp$Mi1Index)]
  Mi2 = temp$Mi2[,1:max(temp$Mi2Index)]
  Mi1Index = temp$Mi1Index
  Mi2Index = temp$Mi2Index
  rm(temp)
  
  #--- Helper functions:
  {
    QGradSU = function(S,U){
      UGrad = cbind(-tauS*S,-tauU*U)
      
      eUiWk = exp(S%*%matrix(1,1,K) + tcrossprod(U,W))
      fuk = colSums(eUiWk)
      
      for(k in 1:K){
        UGrad = UGrad + 
          tcrossprod(Pmki1[,k] - drop(Pk[k]/fuk[k])*eUiWk[,k],c(1,W[k,]))
      }
      
      return(UGrad)
    }
    
    QSU = function(S,U){
      SpUWt = S%*%matrix(1,1,K) + tcrossprod(U,W)
      
      eUiWk = exp(SpUWt)
      fuk = colSums(eUiWk)
      
      ret = sum(Pmk*SpUWt[EE[,1],]) - Pk%*%log(fuk) - 0.5*tauU*sum(U^2) - 0.5*tauS*sum(S^2)
      
      return(ret)
    }
    
    .wrappers = list()
    .wrappers$U = function(x){
      drop(QSU(S = matrix(x[1:n],n,1),
               U=matrix(x[-c(1:n)],n,p)))
    }
    .wrappers$UGrad = function(x){
      temp = QGradSU(S = matrix(x[1:n],n,1),
                     U=matrix(x[-c(1:n)],n,p))
      return(c(temp))
    }
    .wrappers$V = function(x){
      QRV(RR=x[1:n],V=matrix(x[-c(1:n)],n,p),W,pmk=Pmk,Pk,tauR,tauV,EE)
    }
    .wrappers$UV = function(x){
      QUV(SS = x[1:n], RR = x[n+n*p + 1:n], 
          U=matrix(x[n + 1:(n*p)],n,p),
          V=matrix(x[2*n+n*p + 1:(n*p)],n,p),
          W,pmk=Pmk,Pk,Pmki1,Pmki2,tauS,tauR,tauU,tauV,EE,UVApprox,IVApprox)
    }
    .wrappers$UVGrad = function(x){
      dQUV(SS = x[1:n], RR = x[n+n*p + 1:n], 
           U=matrix(x[n + 1:(n*p)],n,p),
           V=matrix(x[2*n+n*p + 1:(n*p)],n,p),
           W,pmk=Pmk,Pk,Pmki1,Pmki2,tauS,tauR,tauU,tauV,EE,UVApprox,IVApprox)
    }
    .wrappers$W = function(x){
      QW(S,R,U,V,W=matrix(x,K,p),pmk=Pmk,Pk,EE,WApprox,IVApprox)
    }
    .wrappers$WGrad = function(x){
      c(dQW(S,R,U,V,W=matrix(x,K,p),Pmk,Pk,EE,WApprox,IVApprox))
    }
    
  }
  #--- End helper functions
  
  ### Initialize
  initSampleM = sample(M,min(5e3,M))
  S = scale(degree(A,mode="out") + 1)
  R = scale(degree(A,mode="in") + 1)
  attributes(R) = attributes(S) = NULL
  
  if(BKNInit){
    cat("\n Initializing via BKN \n")
    starter0 = BKN_eClust(A,K=K, maxIter = 1e2,eps = 1e-3,Plot=FALSE)
    starter = bkn2dir(A,starter0)
    Pmk = as.matrix(starter[,-c(1:2)])
    Pk = colSums(Pmk)
    Pmki1 = getPmki(Pmk,Mi1Index,Mi1)
    Pmki2 = getPmki(Pmk,Mi2Index,Mi2)
    if(p == 2){
      W = seq(0,2*pi,l=K+1)[-1]
      W = cbind(cos(W),sin(W))
    }
    if(p == 3){
      temp = list()
      temp$phi = ( (sqrt(5) + 1)/2 - 1 )*pi*1:K
      temp$theta = asin(-1+2*1:K/K)
      W = cbind(sin(temp$phi)*cos(temp$theta),
                sin(temp$phi)*sin(temp$theta),
                cos(temp$phi))
      rm(temp)
      # plot3D::scatter3D(W[,1],W[,2],W[,3])
      # plot3Drgl::plotrgl()
      # summary(apply(W,1,function(x)sqrt(crossprod(x))))
      # round(dist(W),2)
    }
    if(is.numeric(randUVInit)){
      U = matrix(rnorm(n*p,sd=randUVInit),n,p)
      V = matrix(rnorm(n*p,sd=randUVInit),n,p)
      tauV = 1/randUVInit^2
      tauU = 1/randUVInit^2
      tauS = tauR = 1
    }else{ 
      MM.Spec = getMMSpec(A,apply(starter0$Q,1,which.max),Jeffreys = TRUE)
      U = MM.Spec[,1:K]%*%W
      if(is_directed(A)){ V = U}
      
      temp = function(x){
        evalMargLogLik(S*x,R*x,U,V,W,Pk/sum(Pk),EE,M>5e3,initSampleM)
      }
      Opt0 = list(objective=NaN); count = 0; UB = 1e3; LB = 1e-6
      while(is.na(Opt0$obj) & count < 15){
        options(warn=-1)
        try({
          Opt0 = optimize(temp,c(LB,UB),maximum=TRUE)
        },silent=TRUE)
        options(warn=0)
        UB = UB/2; LB = LB*2; count = count + 1
      }
      if(!is.na(Opt0$obj)){S = Opt0$maximum*S;R = Opt0$maximum*R}
      tauS = 1/var(S); tauR = 1/var(R)
      
      temp = function(x){
        evalMargLogLik(S,R,U*x,V*x,W,Pk/sum(Pk),EE,M>5e3,initSampleM)
      }
      Opt0 = list(objective=NaN); count = 0; UB = 1e3; LB = 1e-6
      while(is.na(Opt0$obj) & count < 15){
        options(warn=-1)
        try({
          Opt0 = optimize(temp,c(LB,UB),maximum=TRUE)
        },silent=TRUE)
        options(warn=0)
        UB = UB/2; LB = LB*2; count = count + 1
      }
      if(!is.na(Opt0$obj))V = U = Opt0$maximum*U
      tauU = tauV = 1/var(c(U))
    }
    a = Pk/sum(Pk)
  }else{
    tauV = 1/randUVInit^2
    tauU = 1/randUVInit^2
    tauR = tauS = 1
    if(p == 2){
      W = seq(0,2*pi,l=K+1)[-1]
      W = cbind(cos(W),sin(W))
    }
    if(p == 3){
      temp = list()
      temp$phi = ( (sqrt(5) + 1)/2 - 1 )*pi*1:K
      temp$theta = asin(-1+2*1:K/K)
      W = cbind(sin(temp$phi)*cos(temp$theta),
                sin(temp$phi)*sin(temp$theta),
                cos(temp$phi))
      rm(temp)
    }
    if(p > 3){
      W = matrix(rnorm(K*p),K,p)
    }
    U = matrix(rnorm(n*p,sd=randUVInit),n,p)
    if(is_directed(A)) V = matrix(rnorm(n*p,sd=randUVInit),n,p)
    S = rnorm(n); R = rnorm(n)
    a = rep(1/K,K)
    Pmk = computePmk(S,R,U,V,W,a,EE)
    Pk = colSums(Pmk)
    Pmki1 = getPmki(Pmk,Mi1Index,Mi1)
    Pmki2 = getPmki(Pmk,Mi2Index,Mi2)
  }
  
  logPost = numeric(maxIter)
  logPost[1] = evalMargPost(S,R,U,V,W,tauS,tauR,tauU,tauV,alph=a,
                            EE,a_s,b_s,a_r,b_r,a_u,a_v,b_u,b_v,a_0,FALSE,initSampleM)
  
  cat("\n Finding reasonable starting U and V \n")
  if(WApprox | UVApprox) for(k in 1:K) IVApprox[,k] = sample(M,nSubsamp,replace=T,prob=Pmk[,k])
  
  if(UVEstMethod == "QN"){
    Opt1 = optim(par=c(S,c(U)),
                 fn = .wrappers$U,
                 gr = .wrappers$UGrad,
                 method="BFGS",
                 control = list(maxit=QNSteps*2,fnscale=-1))
    S = Opt1$par[1:n]
    U = matrix(Opt1$par[n + 1:(n*p)],n,p)
    ### Update V
    Opt2 = optim(par=c(R,c(V)),
                 fn = .wrappers$V,
                 method="Nelder-Mead",
                 control = list(maxit=NMSteps*2,fnscale=-1))
    R = Opt2$par[1:n]
    V = matrix(Opt2$par[n + 1:(n*p)],n,p)
  }
  if(UVEstMethod == "CG"){
    Opt4 = optim(par=c(S,c(U),R,c(V)),
                 fn = .wrappers$UV,
                 gr = .wrappers$UVGrad,
                 method="CG",
                 control = list(maxit=CGSteps*2,fnscale=-1))
    S = Opt4$par[1:n]
    U = matrix(Opt4$par[n + 1:(n*p)],n,p)
    R = Opt4$par[n + n*p + 1:n]
    V = matrix(Opt4$par[2*n +n*p + 1:(n*p)],n,p)
  }
  S = S - mean(S)
  R = R - mean(R)
  
  
  cat("\n Beginning EM algorithm \n")
  pb = txtProgressBar(0,maxIter-1,style=3)
  for(it in 2:maxIter){
    
    ### Update Pmk and Pk (E-step)
    Pmk = computePmk(S,R,U,V,W,alph=a,EE)
    Pk = colSums(Pmk)
    Pmki1 = getPmki(Pmk,Mi1Index,Mi1)
    Pmki2 = getPmki(Pmk,Mi2Index,Mi2)
    
    if(WApprox | UVApprox) for(k in 1:K) IVApprox[,k] = sample(M,nSubsamp,replace=T,prob=Pmk[,k])
    
    ### Update U and V
    if(UVEstMethod == "QN"){
      Opt1 = optim(par=c(S,c(U)),
                   fn = .wrappers$U,
                   gr = .wrappers$UGrad,
                   method="BFGS",
                   control = list(maxit=QNSteps,fnscale=-1))
      S = Opt1$par[1:n]
      U = matrix(Opt1$par[n + 1:(n*p)],n,p)
      ### Update V
      Opt2 = optim(par=c(R,c(V)),
                   fn = .wrappers$V,
                   method="Nelder-Mead",
                   control = list(maxit=NMSteps,fnscale=-1))
      R = Opt2$par[1:n]
      V = matrix(Opt2$par[n + 1:(n*p)],n,p)
    }
    if(UVEstMethod == "CG"){
      Opt4 = optim(par=c(S,c(U),R,c(V)),
                   fn = .wrappers$UV,
                   gr = .wrappers$UVGrad,
                   method="CG",
                   control = list(maxit=CGSteps,fnscale=-1))
      S = Opt4$par[1:n]
      U = matrix(Opt4$par[n + 1:(n*p)],n,p)
      R = Opt4$par[n + n*p + 1:n]
      V = matrix(Opt4$par[2*n +n*p + 1:(n*p)],n,p)
    }
    S = S - mean(S)
    R = R - mean(R)
    
    ### Update tau_S, tau_R, tau_U, and tau_V
    tauS = (a_s + n - 2)/(b_s + sum(S^2))
    tauR = (a_r + n - 2)/(b_r + sum(R^2))
    tauU = (a_u + n*p - 2)/(b_u + sum(U^2))
    tauV = (a_u + n*p - 2)/(b_u + sum(V^2))
    
    
    ### Update W
    Opt3 = optim(par = c(W),
                 fn = .wrappers$W,
                 gr = .wrappers$WGrad,
                 method = "BFGS",
                 control = list(maxit=QNSteps,fnscale=-1))
    W = matrix(Opt3$par,K,p)
    
    ### Update alpha
    a = (a_0 - 1 + Pk)/(K*a_0 + M - K)
    
    ### Compute marginal log posterior (summed over Z)
    logPost[it] = evalMargPost(S,R,U,V,W,tauS,tauR,tauU,tauV,alph=a,
                               EE,a_s,b_s,a_r,b_r,a_u,a_v,b_u,b_v,a_0,FALSE,initSampleM)
    
    ###Check for convergence
    if( ((logPost[it] - logPost[it-1])/abs(logPost[it-1]) < eps ) &
        ( logPost[it] > logPost[it-1] ) ) break
    
    setTxtProgressBar(pb,it)
  }
  
  if(logPost[maxIter] == 0){
    logPost = logPost[1:min(which(logPost == 0)) - 1]
  }
  plot(logPost,type='b',pch=16,
       xlab="Iteration",ylab="Log posterior")
  
  
  return(list(estimates = list(S=S,
                               R=R,
                               U = U,
                               V = V,
                               W = W,
                               Pmk = Pmk,
                               alpha = a,
                               tauU = tauU,
                               tauV = tauV,
                               tauS = tauS,
                               tauR = tauR),
              algoParms = list(logPost = logPost[1:it],
                               eps = eps,
                               numIter = it,
                               maxIter = maxIter,
                               Approx = c(W=WApprox,UV=UVApprox),
                               nSubsamp = nSubsamp),
              userInputs = list(K = K,
                                p = p,
                                a_u = a_u,
                                b_u = b_u,
                                a_v = a_v,
                                b_v = b_v,
                                a_0 = a_0)))
}


# BIC ---------------------------------------------------------------------

#' BIC for edge clustering
#' 
#' Compute the BIC for edge clustering object
#' 
#' @param A igraph object
#' @param eCl object returned from 'eClustGEM'
BIC.eClust = function(A,eCl){
  -2*evalMargLogLik(eCl$estimates$S,eCl$estimates$R,
                    eCl$estimates$U,eCl$estimates$V,
                    eCl$estimates$W,eCl$estimates$alpha,
                    as_edgelist(A,FALSE),FALSE,1:2) +
    log(ecount(A))*( 2^(is_directed(A))*(vcount(A)+length(eCl$estimates$U)) + length(eCl$estimates$W) + 
                      length(eCl$estimates$alpha)-1 )
}


#' AIC for edge clustering
#' 
#' Compute the AIC for edge clustering object
#' 
#' @param A igraph object
#' @param eCl object returned from 'eClustGEM'
AIC.eClust = function(A,eCl){
  -2*evalMargLogLik(eCl$estimates$S,eCl$estimates$R,
                    eCl$estimates$U,eCl$estimates$V,
                    eCl$estimates$W,eCl$estimates$alpha,
                    as_edgelist(A,FALSE),FALSE,1:2) +
    2*( 2^(is_directed(A))*(vcount(A)+length(eCl$estimates$U)) + length(eCl$estimates$W) + 
                       length(eCl$estimates$alpha)-1 )
}

#' AICc for edge clustering
#' 
#' Compute the AICc for edge clustering object
#' 
#' @param A igraph object
#' @param eCl object returned from 'eClustGEM'
AICc.eClust = function(A,eCl){
  nParam = ( 2^(is_directed(A))*(vcount(A)+length(eCl$estimates$U)) + length(eCl$estimates$W) + 
               length(eCl$estimates$alpha)-1 )
  AIC.eClust(A,eCl) + (2*nParam^2 + 2*nParam)/(log(ecount(A)) - nParam - 1)
}


# Ball, Karrer, and Newman 2011 -------------------------------------------

BKN_eClust = function(A,K=2, maxIter = 1e3,eps = 1e-4,Plot=TRUE){
  if(is_directed(A)){
    A1 <- as.undirected(A,"collapse")
  }else{
    A1 <- A
  }
  
  AMat = as_adjacency_matrix(A1,sparse=TRUE)
  ELshort = as_edgelist(A1,FALSE)
  EL = cbind(AMat@i+1,rep(seq_along(diff(AMat@p)),diff(AMat@p)))
  n = vcount(A)
  M = ecount(A)
  
  evalBKMLik = function(th){
    thTht = tcrossprod(th)
    sum(log(thTht[ELshort])) - sum(thTht[upper.tri(thTht)])
  }
  
  theta = matrix(runif(n*K),n,K)
  # QQ = vegetarian::normalize.rows(theta[EL[,1],]*theta[EL[,2],])
  QQ = as(theta[EL[,1],]*theta[EL[,2],],"dgeMatrix")
  QQ = Diagonal(x=1/rowSums(QQ))%*%QQ
  
  QList = list()
  for(k in 1:K) QList[[k]] = AMat
  
  logLik = numeric(maxIter)
  logLik[1] = evalBKMLik(theta)
  
  
  pb = txtProgressBar(0,maxIter,style=3)
  for(it in 2:maxIter){
    
    # QQ = vegetarian::normalize.rows(theta[EL[,1],]*theta[EL[,2],])
    QQ = theta[EL[,1],]*theta[EL[,2],]
    QQ = Diagonal(x=1/rowSums(QQ))%*%QQ
    for(k in 1:K){
      QList[[k]]@x = QQ[,k]
      theta[,k] = rowSums(QList[[k]])
    }
    theta = theta %*%Diagonal(x = 1/sqrt(colSums(theta)))
    
    
    logLik[it] = evalBKMLik(theta)
    
    
    setTxtProgressBar(pb,it)
    if((logLik[it] - logLik[it-1])/abs(logLik[it-1]) < eps) break
  }
  logLik = logLik[1:it]
  if(Plot){
    plot(logLik,type='b',pch=16,
         xlab="Iteration",ylab="Log likelihood")
  }
  
  
  return(list(theta=theta,
              Q=QQ[-which(EL[,1] > EL[,2]),],
              logLik=logLik,
              nIter = it))
}

bkn2dir = function(A,bknObj){
  A1 <- as.undirected(A,"collapse")
  K = ncol(bknObj$Q)
  A0.Mat = as_adjacency_matrix(A,"both",sparse=TRUE)
  ELOrig = as.data.frame(as_edgelist(A,F))
  EL = data.frame(Source = A0.Mat@i+1,
                  Target = rep(seq_along(diff(A0.Mat@p)),diff(A0.Mat@p)),
                  matrix(0,ecount(A),K,dimnames=list(NULL,paste("cluster",1:K,sep="_"))))
  
  for(k in 1:K){
    A1 <- set_edge_attr(A1,paste0("cluster",k),value=bknObj$Q[,k])
    A1.Mat = as_adjacency_matrix(A1,"both",attr=paste0("cluster",k),sparse=TRUE)
    A0 <- A0.Mat*A1.Mat
    EL[,2+k] = A0@x
  }
  
  colnames(ELOrig) = c("Source","Target")
  ELNew = left_join(ELOrig,EL,by = c("Source","Target"))
  
  return(ELNew)
}
