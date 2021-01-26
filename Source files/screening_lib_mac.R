# R library for screening methods
library(mvtnorm)
library(igraph)
library(Matrix)
library(glmnet)
library(parallel)
library(lattice)
library(geoR)

#my.cols = rgb(read.table(file="~/Dropbox/Code/C++/results/Rcode/my.heat.map"))

sapply_pb <- function(X, FUN, ...){
  env <- environment()
  pb_Total <- length(X)
  counter <- 0
  pb <- txtProgressBar(min = 0, max = pb_Total, style = 3)
  
  wrapper <- function(...){
    curVal <- get("counter", envir = env)
    assign("counter", curVal +1 ,envir=env)
    setTxtProgressBar(get("pb", envir=env), curVal +1)
    FUN(...)
  }
  res <- sapply(X, wrapper, ...)
  close(pb)
  res
}


mcsapply <- function(X, FUN, ...){
  res <- simplify2array(mclapply(X, FUN, ...))
  return(res)
}

simul_dat = function(dat,trueBeta,R2=0.9,sigmaY=NULL){
  eta = dat$X%*%trueBeta
  if(is.null(sigmaY)){
    sigmaY = sqrt(var(eta)*(1-R2)/R2)
  }
  else{
    R2 = var(eta)
    R2 = R2/(R2+sigmaY^2)
  }
  
  dat$Y = eta + rnorm(nrow(dat$X),sd=sigmaY)
  dat$sigmaY = sigmaY
  dat$trueBeta = trueBeta
  dat$R2 = 0.9
  return(dat)
}
  


simul_X_indep = function(n, p, mu = 0, sigma = 1){
  #Simulate design matrix from an independent model
  return(list(X=matrix(rnorm(n*p,mean=mu,sd=sigma),nrow=n,ncol=p),mu=mu,sigma=sigma))
}

simul_X_comp_sym = function(n, p, mu = 0, sigma = 1, rho = 0.5){
  #Simulate design matrix from an compound symmetry model
  if(p<=100){
    Sigma = sigma*sigma*(diag(p)*(1-rho) + matrix(rho,nrow=p,ncol=p))
    X = rmvnorm(n = n, mean = rep(mu,p), sigma=Sigma)
  }
  else{
    if(rho>0)
      X = matrix(rnorm(n*p,mean=mu,sd=sigma),nrow=n,ncol=p) + matrix(rnorm(n,mean=0,sd=sigma*sqrt(rho)),nrow=n,ncol=p)
    else{
    if(rho==0)
      X = matrix(rnorm(n*p,mean=mu,sd=sigma),nrow=n,ncol=p)
    else
      stop("Dimension is too large and rho is negative!")
    }
  }
  
  return(list(X = X, n=n,p=p,mu=mu,sigma=sigma,rho=rho))
} 

simul_X_block_diag_random_graph = function(n, p, blockSize = 10, 
                              connectProb = 0.2, connectSignal=1){
  #Simulate design matrix from block diagnal structures
  #Within each block, we use the random graph to construct the 
  #precision matrix
  
  m = floor(p/blockSize)
  X = NULL
  
  for(i in 1:m){
    graphStruct = erdos.renyi.game(blockSize,connectProb)
    blockMat = as.matrix(get.adjacency(graphStruct))*connectSignal
    blockCovMat = solve(convert_sym_pd(blockMat))
    blockX = rmvnorm(n,sigma=blockCovMat)
    if(i == 1)
      covMat = blockCovMat
    else
      covMat = as.matrix(bdiag(covMat,blockCovMat))
    
    X = cbind(X,blockX)
    cat("block ", i, "\n")
    flush.console()
  }
  
  perm_idx = sample(1:p,p,replace=FALSE)
  od_idx =order(perm_idx)
  
  return(list(X=X[,perm_idx], covMat=covMat[perm_idx,perm_idx], 
              perm_idx=perm_idx,od_idx=od_idx,blockSize=blockSize,n=n,
              p=p,connectProb=connectProb,connectSignal=connectSignal))
}


simul_X_block_diag_AR1 = function(n, p, blockSize = 10, rho=0.9,sigma=1){
  #Simulate design matrix from block diagnal structures
  #Within each block, we have the AR1 correlation structure
  m = floor(p/blockSize)
  X = NULL
  sigma2 = sigma*sigma;
  group = NULL
  for(i in 1:m){
    blockCovMat = corr_AR1(rho=rho,d=blockSize)*sigma2
    blockX = rmvnorm(n,sigma=blockCovMat)
    if(i == 1)
      covMat = blockCovMat
    else
      covMat = as.matrix(bdiag(covMat,blockCovMat))
    X = cbind(X,blockX)
    group = c(group, rep(i,length=ncol(blockX)))
    cat("block ", i, "\n")
    flush.console()
    
  }
  perm_idx = sample(1:p,p,replace=FALSE)
  od_idx =order(perm_idx)
  return(list(X=X[,perm_idx], covMat=covMat[perm_idx,perm_idx], group = group[perm_idx],
              perm_idx=perm_idx,od_idx=od_idx,n=n,p=p,blockSize=blockSize,
              rho=rho,sigma=sigma))
}


simul_X_block_diag_comp_sym = function(n, p, blockSize = 10, rho=0.9,sigma=1){
  #Simulate design matrix from block diagnal structures
  #Within each block, we have the AR1 correlation structure
  m = floor(p/blockSize)
  X = NULL
  sigma2 = sigma*sigma;
  group = NULL
  for(i in 1:m){
    
    blockCovMat = sigma2*(diag(blockSize)*(1-rho) + matrix(rho,nrow=blockSize,ncol=blockSize))
    blockX = rmvnorm(n,sigma=blockCovMat)
    if(i == 1)
      covMat = blockCovMat
    else
      covMat = as.matrix(bdiag(covMat,blockCovMat))
    X = cbind(X,blockX)
    group = c(group, rep(i,length=ncol(blockX)))
    cat("block ", i, "\n")
    flush.console()
  }
  perm_idx = sample(1:p,p,replace=FALSE)
  od_idx =order(perm_idx)
  return(list(X=X[,perm_idx], covMat=covMat[perm_idx,perm_idx], group = group[perm_idx],
              perm_idx=perm_idx,od_idx=od_idx,n=n,p=p,blockSize=blockSize,
              rho=rho,sigma=sigma))
}

fast_simul_X_block_diag_AR1 = function(n, p, blockSize = 10, rho=0.9,sigma=1,ifperm=FALSE){
  #Simulate design matrix from block diagnal structures
  #Within each block, we have the AR1 correlation structure
  m = floor(p/blockSize)
  X = NULL
  sigma2 = sigma*sigma;
  group = rep(1:m,each=blockSize)
  blockCovMat = corr_AR1(rho=rho,d=blockSize)*sigma2
  L = chol(blockCovMat)
  list_L = lapply(1:m,function(i) L)
  big_L=bdiag(list_L)
  Z = matrix(rnorm(n*p),nrow=n,ncol=p)
  X = as.matrix(Z%*%big_L)
  if(ifperm){
  perm_idx = sample(1:p,p,replace=FALSE)
  }
  else
  {
    perm_idx = 1:p
  }
  od_idx =order(perm_idx)
  return(list(X=X[,perm_idx], group = group[perm_idx],
              perm_idx=perm_idx,od_idx=od_idx,n=n,p=p,blockSize=blockSize,
              rho=rho,sigma=sigma))
}



fast_simul_X_block_diag_comp_sym = function(n, p, blockSize = 10, rho=0.9,sigma=1,ifperm=FALSE){
  #Simulate design matrix from block diagnal structures
  #Within each block, we have the AR1 correlation structure
  m = floor(p/blockSize)
  X = NULL
  sigma2 = sigma*sigma;
  group = rep(1:m,each=blockSize)
  blockCovMat = (diag(blockSize)*(1-rho)+matrix(rho,nrow=blockSize,ncol=blockSize))*sigma2
  L = chol(blockCovMat)
  list_L = lapply(1:m,function(i) L)
  big_L=bdiag(list_L)
  Z = matrix(rnorm(n*p),nrow=n,ncol=p)
  X = as.matrix(Z%*%big_L)
  if(ifperm){
    perm_idx = sample(1:p,p,replace=FALSE)
  }
  else{
    perm_idx = 1:p
  }
  od_idx =order(perm_idx)
  return(list(X=X[,perm_idx], group = group[perm_idx],
              perm_idx=perm_idx,od_idx=od_idx,n=n,p=p,blockSize=blockSize,
              rho=rho,sigma=sigma))
}


simul_spat_dat_comp_sym = function(n,m,rho,sigma=1,seed=NULL){
  if(!is.null(seed)){
    set.seed(seed)
  }
  p = m*m
  xgrid = as.matrix(expand.grid(seq(0,1,length=m),seq(0,1,length=m)))
  gridsize = sqrt(mean((xgrid[1,]-xgrid[2,])^2))
  X = matrix(rnorm(n*p,sd=sqrt(1-rho)),nrow=n,ncol=p)+matrix(rnorm(n,sd=rho),nrow=n,ncol=p)
  # beta = ifelse(abs(xgrid[,1]-0.5)+abs(xgrid[,2]-0.5)<0.1,signal,0)
  # delta = as.numeric(beta!=0)
  # trueidx = which(delta==1)
  # mu = X%*%beta
  # var_mu = var(mu[,1])
  # sigma2 = var_mu*(1/R2-1)
  # Y = X[,trueidx]%*%beta[trueidx]+rnorm(n,sd=sqrt(sigma2))
  return(list(X=X,p=p,m=m,n=n,xgrid=xgrid,
              gridsize=gridsize,sigma2=sigma*sigma))
}

simul_spat_dat_group_RF = function(n,m,rho,sigma2=1,seed=NULL,group){
  if(!is.null(seed)){
    set.seed(seed)
  }
  p = m*m
  
  xgrid = as.matrix(expand.grid(seq(0,1,length=m),seq(0,1,length=m)))
  gridsize = sqrt(mean((xgrid[1,]-xgrid[2,])^2))
  phi = gridsize/(-log(rho))
  uni_group = unique(group)
  X = matrix(NA, nrow=n,ncol=p)
  for(i in 1:length(uni_group)){
    group_idx=which(group==uni_group[i])
    temp=mclapply(1:n,function(i) return(grf(grid=xgrid[group_idx,],cov.pars=c(sigma2, phi),messages = FALSE)$data),mc.cores=7)
    X[,group_idx] = t(simplify2array(temp, higher = TRUE))
  }
  
  
  #X = matrix(rnorm(n*p,sd=sqrt(1-rho)),nrow=n,ncol=p)+matrix(rnorm(n,sd=rho),nrow=n,ncol=p)
  # beta = ifelse(abs(xgrid[,1]-0.5)+abs(xgrid[,2]-0.5)<0.1,signal,0)
  # delta = as.numeric(beta!=0)
  # trueidx = which(delta==1)
  # mu = X%*%beta
  # var_mu = var(mu[,1])
  # sigma2 = var_mu*(1/R2-1)
  # Y = X[,trueidx]%*%beta[trueidx]+rnorm(n,sd=sqrt(sigma2))
  return(list(X=X,p=p,m=m,n=n,xgrid=xgrid,
              gridsize=gridsize,sigma2=sigma2))
}



simul_spat_dat_RF = function(n,m,rho,sigma2=1,seed=NULL){
  if(!is.null(seed)){
    set.seed(seed)
  }
  p = m*m
  
  xgrid = as.matrix(expand.grid(seq(0,1,length=m),seq(0,1,length=m)))
  gridsize = sqrt(mean((xgrid[1,]-xgrid[2,])^2))
  phi = gridsize/(-log(rho))
  temp=mclapply(1:n,function(i) return(grf(grid=xgrid,cov.pars=c(sigma2, phi),messages = FALSE)$data),mc.cores=7)
  X = t(simplify2array(temp, higher = TRUE))
  
  #X = matrix(rnorm(n*p,sd=sqrt(1-rho)),nrow=n,ncol=p)+matrix(rnorm(n,sd=rho),nrow=n,ncol=p)
  # beta = ifelse(abs(xgrid[,1]-0.5)+abs(xgrid[,2]-0.5)<0.1,signal,0)
  # delta = as.numeric(beta!=0)
  # trueidx = which(delta==1)
  # mu = X%*%beta
  # var_mu = var(mu[,1])
  # sigma2 = var_mu*(1/R2-1)
  # Y = X[,trueidx]%*%beta[trueidx]+rnorm(n,sd=sqrt(sigma2))
  return(list(X=X,p=p,m=m,n=n,xgrid=xgrid,
              gridsize=gridsize,sigma2=sigma2))
}


convert_sym_pd = function(X,eps=1){
  #convert a symmetric matrix into positive denfinite matrix
  #without changing the off diagnal elements
  temp = eigen(X)
  d=abs(min(temp$values))+eps
  return(X+diag(nrow(X))*d)
}

corr_AR1 = function(rho=0.9,d=10){
  corrMat = matrix(NA,nrow=d,ncol=d)
  corrMat = rho^abs(row(corrMat)-col(corrMat))
  return(corrMat)
}



block_screening = function(X,thresh){
  adjmat = (abs(fast_cor_mat(X))>thresh)
  edgeList = which(adjmat==TRUE,arr.ind=TRUE)
  tempIdx =  which(edgeList[,1]>edgeList[,2])
  edges = c(t(edgeList[tempIdx,]))
  g = graph(edges,n=ncol(X),directed=FALSE)
  res=clusters(g)
  return(res$membership)
}



block_screening_positive = function(X,thresh){
  adjmat = (fast_cor_mat(X)>thresh)
  edgeList = which(adjmat==TRUE,arr.ind=TRUE)
  tempIdx =  which(edgeList[,1]>edgeList[,2])
  edges = c(t(edgeList[tempIdx,]))
  g = graph(edges,n=ncol(X),directed=FALSE)
  res=clusters(g)
  return(res$membership)
}



cov.group = function(X, up_thresh=0.9,step_thresh=0.05,low_thresh=0.3,verbose=FALSE){
  cor_mat = abs(fast_cor_mat(X))
  p = ncol(X)
  num_groups = p
  if(up_thresh==low_thresh)
    thresh_list = low_thresh
  else
    thresh_list = seq(low_thresh,up_thresh,step_thresh)
  bic_list = rep(NA, length=length(thresh_list))
  num_list = rep(NA,length=length(thresh_list))
  max_num = rep(NA, length=length(thresh_list))
  res_list = list()
  for(i in 1:length(thresh_list)){
    thresh = thresh_list[i]
    adjmat = (cor_mat>thresh)
    edgeList = which(adjmat==TRUE,arr.ind=TRUE)
    tempIdx =  which(edgeList[,1]>edgeList[,2])
    edges = c(t(edgeList[tempIdx,]))
    g = graph(edges,n=ncol(X),directed=FALSE)
    res_list[[i]]=clusters(g)
    max_num[i]=max(table(res_list[[i]]$membership))
    num_groups = length(unique(res_list[[i]]$membership))
    num_list[i] = num_groups
    bic_list[i] = compute.block.cov.EBIC(X,res_list[[i]]$membership)
    if(verbose==TRUE){
      cat(num_groups,thresh,bic_list[i],"\n")
    }
    if(num_groups==p)
      break
    
    flush.console()
  }
  
  #idx_set = which(num_list>1)
  idx_set = which(max_num<nrow(X))
  #idx_set = 1:length(thresh_list)
  min_bic = which.min(bic_list[idx_set])
  return(res_list[[idx_set[min_bic]]]$membership)
}

compute.block.cov.EBIC = function(X,group,gam = 1.0){
  group_idx = unique(group)
  X_bar = apply(X,2,mean)
  n = nrow(X)
  log_p = log(ncol(X))
  log_n = log(n)
  ones = matrix(1,nrow=n,ncol=1)
  bic = 0
  for(i in 1:length(group_idx)){
    group_list = which(group==group_idx[i])
    p_g = length(group_list)
    if(p_g>1){
      group_cov = fast.cov(X[,group_list])
      bic = bic + as.numeric(determinant(group_cov,logarithm = TRUE)$modulus)+p_g 
      + p_g*(p_g+1)/2*(log_n+4*gam*log_p)/n
      #bic = bic + log(det(group_cov))+p_g + p_g*(p_g+1)/2*log_n/n
    }
    else{
      group_cov = var(X[,group_list])
      bic = bic + (log(group_cov)+1)  + (log_n+4*gam*log_p)/n
    }
  }
  return(bic)
}

normalize = function(x){
  return(sqrt(length(x))/sqrt(length(x)-1)*(x-mean(x))/sd(x))
}

fast_cor_mat = function(X){
  Z = apply(X,2,normalize)
  return(t(Z)%*%Z/nrow(X))
}


set_beta_block =function(p,blockSize=100){
  trueBeta <- rep(0, p)
  trueBeta[1] =0.5
  trueBeta[2] = -0.5
  trueBeta[blockSize+1] = 0.5
  trueBeta[blockSize+2] = -0.5
  trueBeta[2*blockSize+1] = 0.5
  trueBeta[3*blockSize+1] = 0.5
  trueBeta[4*blockSize+1] = -0.5
  trueBeta[5*blockSize+1] =0.5
  trueBeta[6*blockSize+1] = -0.5
  trueBeta[7*blockSize+1] = 0.5
  return(trueBeta)
}




GD_screening = function(Y,X,block_ID,tol=1.0e-4,rate=0.1,maxit=10,maxit2=5,maxit3=100){

  #maxit:  selected groups
  #maxit2: screening iteration; if =1, then no iteration
  #maxit3: boosting iteration
  
  z = apply(X,2,normalize)
  Y_mean=Y-mean(Y)
  
  ####GDscreen
  number.groups=length(unique(block_ID))
  
  select_index=NULL
  track1=NULL  ##track which variable is screened (before cleaning) at each step
  track2=NULL  ##track which variable is selected after clearning at each step
  beta=rep(0,p)
  key=0
  number_select_group=0
  repeat{
    key=key+1
    beta_old=beta
    Y_hat=z%*%beta
    
    residual=Y-Y_hat
    
    group=0
    beta2=rep(0,p)
    GD=rep(0,number.groups)
    select=rep(0,number.groups)
    while (group<number.groups){
      group= group+1
      z_temp=as.matrix(z[,block_ID==group],nrow=length(Y))
      p_temp=sum(block_ID==group)
      if(ncol(z_temp)>2){
        aa<-cv.glmnet(z_temp, residual,standardize=FALSE,family="gaussian",alpha=0.5)  #default alpha=1 Lasso
        aa3<-glmnet(z_temp, residual,standardize=FALSE,family="gaussian",alpha=0.5)  #default alpha=1 Lasso
        final=which(aa$lambda==aa$lambda.min)
        beta_lasso <- coef(aa3, s=aa$lambda[final]) 
        beta_lasso2=matrix(beta_lasso, ncol=1)
        beta_lasso3=beta_lasso2[-1,]
        beta_temp=beta_lasso3
      }
      else{
        beta_temp=t(z_temp)%*%residual/length(Y)
      }
        
        beta2[block_ID==group]=beta_temp
      
        GD[group]=(t(t(z_temp)%*%matrix(residual,ncol=1))%*%beta_temp)/sqrt(p_temp)
        select[group]=sum(abs(beta_temp)>0)
    }
    
    select[select==0]=p  ##account for zero selection
    # group_select=rev(order(GD/sqrt(select)))[1:max(floor((maxit-number_select_group)*1),1)] #arbitrary threshold
    group_select=rev(order(GD))[1:max(floor((maxit-number_select_group)*1),1)] #arbitrary threshold
    
    beta2[block_ID%!in%group_select]=0
    beta_key=which((abs(beta2)>0)==1)
    beta_key3=unique(c(beta_key,select_index))
    
    track1[[key]]=beta_key3  ##track which variable is selected at each step
    z_temp=z[,beta_key3]  ##build covariate for clearning step
    group_temp=block_ID[beta_key3]   ##build covariate  group index for clearning step
    
    
    ####GDboost for cleaning step (Sparse Group Lasso maybe better; however, exisiting R packages for SGL seem not good)
    p_select_number=length(beta_key3)
    group_select_key=unique(group_temp)
    group_select_number=length(group_select_key)
    
    beta_boost=rep(0,p_select_number)
    select_index2=NULL
    key2=0
    repeat{
      key2=key2+1
      Y_hat2=z_temp%*%beta_boost
      residual2=Y_mean-Y_hat2
      
      group2=0
      beta22=rep(0,p_select_number) ##group-wise update for beta
      GD2=rep(0,group_select_number)
      select2=rep(0,group_select_number)
      
      while (group2<group_select_number){
        
        group2= group2+1
        group3=group_select_key[group2]
        z_temp2=z_temp[,group_temp==group3]
        p_temp2=sum(group_temp==group3)
        
        if (p_temp2>1) {
          
          aa<-cv.glmnet(z_temp2, residual2,standardize=FALSE,family="gaussian",alpha=0.5)  #default alpha=1 Lasso
          aa3<-glmnet(z_temp2, residual2,standardize=FALSE,family="gaussian",alpha=0.5)  #default alpha=1 Lasso
          
          final=which(aa$lambda==aa$lambda.min)
          beta_lasso <- coef(aa3, s=aa$lambda[final]) 
          beta_lasso2=matrix(beta_lasso, ncol=1)
          beta_lasso3=beta_lasso2[-1,]
          
          
          beta_temp=beta_lasso3
        } else {
          z_temp2=matrix(z_temp2,ncol=1)
          beta_temp= solve(t(z_temp2)%*%z_temp2)%*%(t(z_temp2)%*%matrix(residual2,ncol=1))
        }
        
        
        beta22[group_temp==group3]=beta_temp
        #GD2[group2]=t(t(z_temp2)%*%matrix(residual2,ncol=1))%*%beta_temp
        GD2[group2]=(t(t(z_temp2)%*%matrix(residual2,ncol=1))%*%beta_temp)/sqrt(p_temp2)
        select2[group2]=sum(abs(beta_temp)>0)
      }
      
      select2[select2==0]=p
      #group_select2=which.max(GD2/sqrt(select2)) #better thangd2
      group_select2=which.max(GD2)
      # track=c(track,group_select2)
      beta22[group_temp%!in%group_select_key[group_select2]]=0
      
      beta_boost=beta_boost+rate*beta22
      
      beta=rep(0,p)
      beta[beta_key3]=beta_boost  ##new beta for current iteration
      
      
      ##summary of current iteration
      
      select_group=unique(block_ID[which((abs(beta)>0)==1)]) ##selected group
      select_index=which((abs(beta)>0)==1) ##selected variable
      number_select_p=sum(select_index)  ##number of selected variable
      number_select_group=length(select_group) ####number of selected group
      
      #if (max(abs(beta-beta_old))<tol) break
      if(key2>=maxit3) break
    }
    track2[[key]]=select_index
    if (max(abs(beta-beta_old))<tol) break
    if (number_select_group>=maxit)  break
    if(key>=maxit2) break
  }
  
  
  return(beta)
}


'%!in%' = function(x,y){
  !('%in%'(x,y))
}

simple.lars = function(y,X,method="BIC",K=5,gamm=0.5){
  elapsed = proc.time()[3]
  res = list()
  p = ncol(X)
  n = nrow(X)
  
  larsres = lars(x=X,y=y,intercept=FALSE,use.Gram=FALSE)
  if(method=="cv"){
    cvres = cv.lars(x=X,y=y,K=K,plot.it=FALSE,index=seq(0,1,length=nrow(larsres$beta)))  
    temp = matrix(as.numeric(larsres$beta[which.min(cvres$cv),]))
  }
  if(method=="Cp"){
    temp = matrix(as.numeric(larsres$beta[which.min(larsres$Cp),]))
  }
  if(method=="BIC"){
    parlen = apply(larsres$beta!=0,1,sum)
    bic = n*log(larsres$RSS/n)+log(n)*parlen+2*gamm*(lgamma(p+1)-lgamma(parlen+1)-lgamma(p-parlen+1))
    temp = matrix(as.numeric(larsres$beta[which.min(bic),]))
  }
  res = list()
  res$selected_idx = which(temp!=0)
  res$betacoef = temp
  res$elapsed = proc.time()[3]-elapsed
  
  #  }
  #   else{
  #     res$ix = 1
  #     res$betacoef = sum(X*y)/sum(X^2)
  #   }
  return(res)
}


simple_ISIS = function(Y,X){
  elapsed = proc.time()[3]
  z = apply(X,2,normalize)
  res=try(SIS(z, Y, family="gaussian", penalty="lasso", tune="bic",varISIS="aggr", seed=41),silent=TRUE)
  beta_ISIS=rep(0,p)
  if(class(res)[1]!="try-error"){
      beta_ISIS[res$ix]=res$coef.est[-1]
  }
  return(list(betaEst = beta_ISIS,elapsed=as.numeric(proc.time()[3]-elapsed)))
}

simple_GD = function(Y,X,thresh=0.3){
  elapsed = proc.time()[3]
  blockID = block_screening(X,thresh=thresh)
  betaEst = GD_screening(Y,X, blockID)
  return(list(betaEst = betaEst, elapsed=as.numeric(proc.time()[3]-elapsed)))
}


pred_err = function(Y,X,beta_coef){
  epsilon = Y - X%*%beta_coef
  return(mean(epsilon^2))
}

summary_fitting = function(trueBeta,estBeta,method="",datset="",compTime=NULL){
  FP = length(which(trueBeta==0 & estBeta!=0))
  FN = length(which(trueBeta!=0 & estBeta==0))
  Se = length(which(trueBeta!=0 & estBeta!=0))/length(which(trueBeta!=0))
  Sp = length(which(trueBeta==0 & estBeta==0))/length(which(trueBeta==0))
  Mse = mean((trueBeta-estBeta)^2)
  if(!is.null(compTime))
    output = data.frame(list(Dataset=datset,Method=method,FP=FP,FN=FN,Se=Se,Sp=Sp,Mse=Mse,CompTime=compTime))
  else
    output = data.frame(list(Dataset=datset,Method=method,FP=FP,FN=FN,Se=Se,Sp=Sp,Mse=Mse))
  return(output)
}


simple_Bayes = function(Y,X,alpha=0.05,starting.idx = NULL,
                                                update.lambda=FALSE,max.step=Inf,
                                                k = 1,tau2=2,sigma2=1){
  elapsed = proc.time()[3]
  p = ncol(X)
  n = nrow(X)
  if(is.null(starting.idx))
    starting.idx = 1:p
  prob.thresh = qnorm((1-alpha*0.5))
  idx0 = starting.idx
  if(length(idx0)>0){
    Wfit = Bayes_marginal_screening(Y,X,includingIdx = idx0,k=k,tau2=tau2,sigma2=sigma2)$selection
    if(!all(abs(Wfit[,5])>prob.thresh)){
      if(nrow(Wfit)>n)
        res = Mclust(Wfit[,4],G=2)
      else
        res = Mclust(Wfit[,4],G=2)
      cluster1 = which(res$classification==1)
      cluster2 = which(res$classification==2)
      if(mean(Wfit[cluster1,5])>mean(Wfit[cluster2,5]))
        idx1 = idx0[cluster1]
      else
        idx1 = idx0[cluster2]
    }
  }
  step = 1
  while((length(idx0)>length(idx1))&(step<max.step)){
    idx0 = idx1
    if(length(idx0)>0){
      Wfit = Bayes_marginal_screening(Y,X,,includingIdx = idx0,k=k,tau2=tau2,sigma2=sigma2)$selection
      if(!all(abs(Wfit[,5])>prob.thresh)){
        if(nrow(Wfit)>n)
          res = Mclust(Wfit[,4],G=2)
        else
          res = Mclust(Wfit[,4],G=2)
        cluster1 = which(res$classification==1)
        cluster2 = which(res$classification==2)
        if(mean(Wfit[cluster1,5])>mean(Wfit[cluster2,5]))
          idx1 = idx0[cluster1]
        else
          idx1 = idx0[cluster2]
      }
      
    }
    step=step+1
  }
  
  betaEst = rep(0,length=p)
  betaEst[idx1] = solve(t(X[,idx1])%*%X[,idx1])%*%t(X[,idx1])%*%Y
  return(list(betaEst = betaEst,elapsed = as.numeric(proc.time()[3]-elapsed)))
  
}

Bayes_marginal_screening = function(Y,X,computingIdx=NULL,includingIdx=1:ncol(X),
                                    k=1,tau2=2,sigma2=1){
    if(is.null(computingIdx))
      computingIdx = 1:length(includingIdx)
    n = nrow(X)
    start = proc.time()[3]
    Xmat = X[,includingIdx[computingIdx]]%*%t(X[,includingIdx[computingIdx]])
    
    OmegaAll = solve(Xmat+sigma2/tau2*diag(n))
    
    denbeta = sapply_pb(computingIdx,function(j){
      temp = X[,includingIdx[j]]%*%OmegaAll
      kappa2 = temp%*%X[,includingIdx[j]]
      nu = temp%*%Y 
      A = sum(temp*temp)
      return(c(nu,kappa2,A))
    })
    nu = denbeta[1,]
    kappa2=tau2*(1-denbeta[2,])
    A = sigma2*denbeta[3,]
    
    Zstat = abs(nu/sqrt(A))
    
    elapsed = proc.time()[3]-start
    
    return(list(selection=data.frame(index=includingIdx[computingIdx],
                      nu=nu,
                      kappa2=kappa2,
                      absnu = abs(nu),
                      Zstat = Zstat),elapsed=elapsed))
}

community_detection_score = function(A,K=2){
  #Community detection by score 
  # A is an adjacent matrix with diagnal term being equal to zero
  eig =eigen(A)
  od = order(abs(eig$values),decreasing=TRUE)
  Tn = log(nrow(A))
  R_hat = eig$vectors[,od[2:K]]/(eig$vectors[,od[1]]+1e-8)
  R_hat = ifelse(R_hat>Tn,Tn,R_hat)
  R_hat = ifelse(R_hat<(-Tn),Tn,-R_hat)
  res = kmeans(R_hat,centers=K)
  return(res)
}



Screen3S=function(y,x,family="gaussian",nsize,lambda){
  
  elapsed = proc.time()[3]
  n=nrow(x)
  p=ncol(x)
  fullset=1:p
  x<- apply(x,2, function(x) (x-mean(x))/sd(x))
  
  ##Step1) Pre-screening screening
  ## find the q variables with the largest magnitude
  q=n 
  beta0=rep(0,p)
  for(j in 1:p){
    Xtemp=cbind(x[,j]);
    glm=glm(y~Xtemp,family=family)
    beta0[j]=abs(glm$coeff)[2];
  }
  
  subset=rep(0,q) 
  sort(beta0,TRUE) ->aa
  for ( i in 1:q){
    subset[i]=which(beta0==aa[i])
  }
  
  #Step 2) Fantope projection
  subset=sort(subset)
  out <- fps(cor(x[,subset]), ndim = nsize)
  #cvout <- cv(out, x[,subset], FUN = cor, verbose = 1)
  #v1 <- coef(out, lambda=cvout$lambda.cv) 
  v <- coef(out, lambda=lambda) 
  
  # select the conditioning set
  maxi=function(x){max(x,na.rm=TRUE)}
  vv=apply(abs(v),1,maxi)
  
  id=NULL
  for(jj in 1:nsize){
    id1=which(abs(vv)==sort(abs(vv), TRUE)[jj])
    id=c(id, id1)
  }
  Cset=subset[unique(id)[1:nsize]]
  
  ## Step3) Perform CSIS
  scanset=fullset[-Cset];
  cSize=length(Cset);
  beta=rep(0,p)
  
  for(j in scanset){
    Xtemp=cbind(x[,Cset],x[,j]);
    glm=glm(y~Xtemp,family=family)
    beta[j]=glm$coeff[cSize+2];
  }
  
  ##replace beta as NA if variables are selected as Cset in Step 2 
  #beta[Cset]= NA
  temp = beta
  temp[Cset] = 0
  selection_order = c(Cset,selected.set(temp,cut_num = p-cSize))
  
  list(screen_stat=abs(beta), selected_order=selection_order, coef=beta, rank=p+1-length(Cset)-rank(beta,na.last="keep"), Cset=Cset,
       elapsed=proc.time()[3]-elapsed)
}




random.group = function(p,G){
  group=sample(1:G,p,replace=TRUE)
  return(group)
}

multi.screening = function(Y,X,group=NULL,family="gaussian"){
  elapsed = proc.time()[3]
  n = nrow(X)
  p = ncol(X)
  if(is.null(group)){
    group = 1:p
  }
  group_idx = unique(group)
  betacoef = rep(NA, length=p)
  for(i in 1:length(group_idx)){
    group_list = which(group==group_idx[i])
    fit = glm(Y~X[,group_list],family=family)
    betacoef[group_list] = coef(fit)[-1]
  }
  selected_order = selected.set(betacoef,cut_num=p)
  return(list(screen_stat=betacoef,selected_order=selected_order,
              elapsed=proc.time()[3]-elapsed))
}

fast.multi.screening = function(Y,X,Z=NULL,group=NULL,family="gaussian",
                                mc.cores=7){
  elapsed = proc.time()[3]
  n = nrow(X)
  p = ncol(X)
  if(is.null(group)){
    group = 1:p
  }
  group_idx = unique(group)
  betacoef = rep(NA, length=p)
  
  
  group_list = mclapply(1:length(group_idx),function(i){ 
    return(which(group==group_idx[i]))
  },mc.cores = mc.cores)
  
  betacoef_list= mclapply(1:length(group_list),function(i){
    if(is.null(Z)){
      fit = glm(Y~X[,group_list[[i]]],family=family)  
      return(coef(fit)[-1])
    }
    else{
      fit = glm(Y~Z+X[,group_list[[i]]],family=family) 
      return(coef(fit)[-c(1:(ncol(Z)+1))])
    }
  },mc.cores = mc.cores)
  
  betacoef[unlist(group_list)] = unlist(betacoef_list)
  
  selected_order = selected.set(betacoef,cut_num=p)
  
  
  return(list(screen_stat=betacoef, selected_order=selected_order, group=group,
              elapsed=proc.time()[3]-elapsed))
}




fast.multi.screening.glmnet = function(Y,X,Z=NULL,group=NULL,group_list=NULL,family="gaussian",
                                mc.cores=7,lambda=0.0){
  elapsed = proc.time()[3]
  n = nrow(X)
  p = ncol(X)
  if(is.null(group) && is.null(group_list)){
    group = 1:p
  }
  group_idx = unique(group)
  #betacoef = rep(NA, length=p)
  betacoef = rep(0, length=p)
  
  if(is.null(group_list)){
    group_list = mclapply(1:length(group_idx),function(i){ 
      return(which(group==group_idx[i]))
    },mc.cores = mc.cores)
  }
  
  
  
  if(is.null(Z)){
  betacoef_list= mclapply(1:length(group_list),function(i){
    fit=glmnet(x = X[,group_list[[i]]], y = Y, family = family, 
           alpha = 0, lambda = lambda, intercept = TRUE)
    #fit = glm(Y~X[,group_list[[i]]],family=family)  
    return(coef(fit)[-1])
  },mc.cores = mc.cores)
  }
  else
  {
    betacoef_list= mclapply(1:length(group_list),function(i){
      fit=glmnet(x = cbind(Z,X[,group_list[[i]]]), y = Y, family = family, 
                 alpha = 0, lambda = lambda, intercept = TRUE)
      #fit = glm(Y~X[,group_list[[i]]],family=family)  
      return(coef(fit)[-c(1:(ncol(Z)+1))])
    },mc.cores = mc.cores)
  }
  
  #print(betacoef_list)
  #print(group_list)
  #betacoef[unlist(group_list)] = unlist(betacoef_list)
  #count = rep(0,length=p)
  for(i in 1:length(group_list)){
     betacoef[group_list[[i]]] = ifelse(abs(betacoef_list[[i]])>abs(betacoef[group_list[[i]]]),betacoef_list[[i]],betacoef[group_list[[i]]])
    #count[group_list[[i]]] = count[group_list[[i]]] + 1
  }
  #pos_idx = which(count>1)
  #betacoef[pos_idx] = betacoef[pos_idx]/count[pos_idx]
  
  
  selected_order = selected.set(betacoef,cut_num=p)
  
  
  return(list(screen_stat=betacoef, selected_order=selected_order,
              elapsed=proc.time()[3]-elapsed))
}

fast.multi.screening.glmnet.RSS = function(Y,X,Z=NULL,group=NULL,group_list=NULL,family="gaussian",
                                       mc.cores=7,lambda=0.0){
  elapsed = proc.time()[3]
  n = nrow(X)
  p = ncol(X)
  if(is.null(group) && is.null(group_list)){
    group = 1:p
  }
  group_idx = unique(group)
  #betacoef = rep(NA, length=p)
  betacoef = rep(0, length=p)
  
  if(is.null(group_list)){
    group_list = mclapply(1:length(group_idx),function(i){ 
      return(which(group==group_idx[i]))
    },mc.cores = mc.cores)
  }
  
  
  
  if(is.null(Z)){
    betacoef_list= mclapply(1:length(group_list),function(i){
      fit=glmnet(x = cbind(1,X[,group_list[[i]]]), y = Y, family = family, 
                 alpha = 0, lambda = lambda, intercept = FALSE)
      #fit = glm(Y~X[,group_list[[i]]],family=family)  
      return(coef(fit)[-c(1,2)]*fit$dev.ratio)
    },mc.cores = mc.cores)
  }
  else
  {
    betacoef_list= mclapply(1:length(group_list),function(i){
      fit=glmnet(x = cbind(Z,X[,group_list[[i]]]), y = Y, family = family, 
                 alpha = 0, lambda = lambda, intercept = TRUE)
      #fit = glm(Y~X[,group_list[[i]]],family=family)  
      return(coef(fit)[-c(1:(ncol(Z)+1))]*fit$dev.ratio)
    },mc.cores = mc.cores)
  }
  
  #print(betacoef_list)
  #print(group_list)
  #betacoef[unlist(group_list)] = unlist(betacoef_list)
  #count = rep(0,length=p)
  for(i in 1:length(group_list)){
    betacoef[group_list[[i]]] = ifelse(abs(betacoef_list[[i]])>abs(betacoef[group_list[[i]]]),betacoef_list[[i]],betacoef[group_list[[i]]])
    #count[group_list[[i]]] = count[group_list[[i]]] + 1
  }
  #pos_idx = which(count>1)
  #betacoef[pos_idx] = betacoef[pos_idx]/count[pos_idx]
  
  
  selected_order = selected.set(betacoef,cut_num=p)
  
  
  return(list(screen_stat=betacoef, selected_order=selected_order,
              elapsed=proc.time()[3]-elapsed))
}



fast.cond.screening = function(Y,X,cond_set=NULL,family="gaussian",
                                mc.cores=7){
  elapsed = proc.time()[3]
  n = nrow(X)
  p = ncol(X)
  betacoef = rep(NA, length=p)
  
  all_set = 1:p
  screen_set = setdiff(all_set,cond_set)
  betacoef_list= mclapply(1:length(screen_set),function(i){
    fit = glm(Y~X[,c(cond_set,screen_set[i])],family=family)  
    return(coef(fit)[-c(1:(1+length(cond_set)))])
  },mc.cores = mc.cores)
  
  betacoef[cond_set] = Inf
  betacoef[screen_set] = unlist(betacoef_list)
  
  selected_order = selected.set(betacoef,cut_num=p)
  
  
  return(list(screen_stat=betacoef, selected_order=selected_order,
              elapsed=proc.time()[3]-elapsed))
}


simple.holp = function(Y,X,family="gaussian"){
  elapsed=proc.time()[3]
  #temp=screening(X, Y, method = 'holp', num.select = ncol(X),family=family)
  #res=list(selected_order=temp$screen,elapsed=proc.time()[3]-elapsed)
   temp = glmnet(X,Y,family=family,alpha=0,lambda=1,intercept = FALSE)
   res=list(screen_stat=as.numeric(temp$beta),selected_order=order(abs(temp$beta),decreasing = TRUE),elapsed=proc.time()[3]-elapsed)
  return(res)
}

self.holp = function(Y,X,Z=NULL,family="gaussian",lambda=0.001){
  elapsed=proc.time()[3]
  model = glmnet(x = cbind(Z,X), y = Y, family = family, 
                 alpha = 0, lambda = lambda, intercept = FALSE)
  coefs = coef(model)[-c(1:(ncol(Z)+1))]
  ranking = sort(abs(coefs), index.return = TRUE, 
                 decreasing = TRUE)
  
  
  res=list(selected_stat=coefs,selected_order=ranking$ix,
           elapsed=proc.time()[3]-elapsed)
  return(res)
}


iterative.holp = function(Y,X,family="gaussian",cut_ratio=0.5){
  elapsed=proc.time()[3]
  p = ncol(X)
  res = screening(X, Y, method = 'holp', num.select = p,family=family)
  holp_order = res$screen
  p_start = floor(cut_ratio*p)
  
    while(p_start>1){
      res = screening(X[,holp_order[1:p_start]], Y, method = 'holp', num.select = p_start,family=family)
      holp_order[1:p_start] = holp_order[1:p_start][res$screen]
      p_start = floor(cut_ratio*p_start)
      cat("iter",iter,p_start,"\n")
      iter = iter + 1
      flush.console()
    }
  return(list(selected_order=holp_order,elapsed=proc.time()[3]-elapsed))
}

iholp = function(Y,X,family="gaussian",ebic.gamma=1){
  elapsed=proc.time()[3]
  p = ncol(X)
  res = screening(X, Y, method = 'holp', family=family,ebic=TRUE,
                  ebic.gamma=ebic.gamma)
  return(list(selected_order=res$screen,elapsed=proc.time()[3]-elapsed))
}

fast.multi.screening.wald = function(Y,X,Z=NULL,group=NULL,family="gaussian",
                                mc.cores=7){
  elapsed = proc.time()[3]
  n = nrow(X)
  p = ncol(X)
  if(is.null(group)){
    group = 1:p
  }
  group_idx = unique(group)
  betacoef = rep(NA, length=p)
  
  
  group_list = mclapply(1:length(group_idx),function(i){ 
    return(which(group==group_idx[i]))
  },mc.cores = mc.cores)
  
  betacoef_list= mclapply(1:length(group_list),function(i){
    if(is.null(Z)){
      fit = summary(glm(Y~X[,group_list[[i]]],family=family))  
      return(fit$coefficients[-1,3])
    }
    else{
      fit = summary(glm(Y~Z+X[,group_list[[i]]],family=family))  
      return(fit$coefficients[-c(1:(ncol(Z)+1)),3])
    }
  },mc.cores = mc.cores)
  
  betacoef[unlist(group_list)] = unlist(betacoef_list)
  
  selected_order = selected.set(betacoef,cut_num=p)
  
  
  return(list(screen_stat=betacoef, selected_order=selected_order,
              elapsed=proc.time()[3]-elapsed))
}

vec.multi.screening = function(Y,X,group=NULL,family="gaussian"){
  elapsed = proc.time()[3]
  n = nrow(X)
  p = ncol(X)
  if(is.null(group)){
    group = 1:p
  }
  group_idx = unique(group)
  betacoef = rep(NA, length=p)
  
  
  group_list = lapply(1:length(group_idx),function(i){ 
    return(which(group==group_idx[i]))
  })
  
  betacoef_list= lapply(1:length(group_list),function(i){
    fit = glm(Y~X[,group_list[[i]]],family=family)  
    return(coef(fit)[-1])
  })
  
  betacoef[unlist(group_list)] = unlist(betacoef_list)
  
  #   for(i in 1:length(group_idx)){
  #     group_list = which(group==group_idx[i])
  #     fit = glm(Y~X[,group_list],family=family)
  #     betacoef[group_list] = coef(fit)[-1]
  #   }
  selected_order = selected.set(betacoef,cut_num=p)
  return(list(screen_stat=betacoef,selected_order=selected_order,
              elapsed=proc.time()[3]-elapsed))
}

prior.cond.multi.screening = function(Y,X,group=NULL,Cset=NULL,family="gaussian",prior_weight = 0.5){
  elapsed = proc.time()[3]
  n = nrow(X)
  p = ncol(X)
  if(is.null(group)){
    group = 1:p
  }
  
  group_idx = unique(group)
  betacoef = rep(NA, length=p)
  for(i in 1:length(group_idx)){
    group_list = which(group==group_idx[i])
    fit = glm(Y~X[,union(group_list,Cset)],family=family)
    pred = predict(fit)
    Y_prior = (1.0-prior_weight)*Y+ prior_weight*pred
    fit_prior = glm(Y_prior~X[,group_list],family=family)
    betacoef[group_list] = coef(fit_prior)[-1]
  }
  selected_order = selected.set(betacoef,cut_num=p)
  return(list(screen_stat=betacoef,selected_order=selected_order,
              elapsed=proc.time()[3]-elapsed))
}



cond.multi.screening = function(Y,X,group=NULL,Cset=NULL,family="gaussian"){
  elapsed = proc.time()[3]
  n = nrow(X)
  p = ncol(X)
  if(is.null(group)){
    group = 1:p
  }
  
  group_idx = unique(group)
  betacoef = rep(NA, length=p)
  for(i in 1:length(group_idx)){
    group_list = which(group==group_idx[i])
    fit = glm(Y~X[,union(group_list,Cset)],family=family)
    betacoef[union(group_list,Cset)] = coef(fit)[-1]
  }
  if(is.null(Cset))
    temp = betacoef
  else
    temp = betacoef[-Cset]
  selected_order = c(Cset,selected.set(temp,cut_num=p-length(Cset)))
  return(list(screen_stat=betacoef,selected_order=selected_order,
              elapsed=proc.time()[3]-elapsed))
}


cov.multi.screening = function(Y,X,family="gaussian"){
  elapsed = proc.time()[3]
  group=cov.group(X,up_thresh=0.6,step_thresh=0.15,low_thresh=0.3)
  multi_res = multi.screening(Y,X,group=group,family=family)
  
  return(list(screen_stat=multi_res$screen_stat,selected_order=multi_res$selected_order, elapsed=proc.time()[3]-elapsed))
}

random.fast.multi.screening = function(Y,X,Z = NULL,family="gaussian",
                                       rep=20, group_size=NULL,verbose=FALSE){
  
  elapsed = proc.time()[3]
  p = ncol(X)
  n = nrow(X)
  
  if(is.null(group_size)){
    group_size = n/5
  }
  G = floor(p/group_size)
  #screen_stat = 0
  screen_stat_rep = matrix(NA, nrow=ncol(X),ncol=rep)
  for(i in 1:rep){
    group = sample(1:G,p,replace=TRUE)
    multi_res = fast.multi.screening(Y,X,Z = Z, group=group,family=family)
    #selected_set = order(abs(multi_res$screen_stat),decreasing = TRUE)[1:n]
    #include_prob = include_prob + as.numeric(is.element(1:p,selected_set))
    #screen_stat = screen_stat + multi_res$screen_stat
    screen_stat_rep[,i] = multi_res$screen_stat
    #screen_stat = ifelse(abs(screen_stat)>abs(multi_res$screen_stat),screen_stat,multi_res$screen_stat)
    if(verbose==TRUE){
      cat(i,"\n")
      flush.console()
    }
  }
  #screen_stat = screen_stat/rep
  screen_stat = apply(screen_stat_rep,1,max)
  selected_order = selected.set(screen_stat,cut_num = p)
  return(list(screen_stat=screen_stat,selected_order=selected_order,
              elapsed = proc.time()[3]-elapsed))
}

random.fast.multi.screening.glmnet.RSS = function(Y,X,family="gaussian",
                                       rep=20, group_size=NULL,verbose=FALSE){
  
  elapsed = proc.time()[3]
  p = ncol(X)
  n = nrow(X)
  
  if(is.null(group_size)){
    group_size = n/5
  }
  G = floor(p/group_size)
  #screen_stat = 0
  screen_stat_rep = matrix(NA, nrow=ncol(X),ncol=rep)
  for(i in 1:rep){
    group = sample(1:G,p,replace=TRUE)
    multi_res = fast.multi.screening.glmnet.RSS(Y,X,group=group,family=family)
    #selected_set = order(abs(multi_res$screen_stat),decreasing = TRUE)[1:n]
    #include_prob = include_prob + as.numeric(is.element(1:p,selected_set))
    #screen_stat = screen_stat + multi_res$screen_stat
    screen_stat_rep[,i] = multi_res$screen_stat
    #screen_stat = ifelse(abs(screen_stat)>abs(multi_res$screen_stat),screen_stat,multi_res$screen_stat)
    if(verbose==TRUE){
      cat(i,"\n")
      flush.console()
    }
  }
  #screen_stat = screen_stat/rep
  screen_stat = apply(screen_stat_rep,1,max)
  selected_order = selected.set(screen_stat,cut_num = p)
  return(list(screen_stat=screen_stat,selected_order=selected_order,
              elapsed = proc.time()[3]-elapsed))
}


spat.random.fast.multi.screening = function(Y,X,xgrid,family="gaussian",
                                       rep=20, group_size=NULL,verbose=FALSE){
  
  elapsed = proc.time()[3]
  p = ncol(X)
  n = nrow(X)
  
  if(is.null(group_size)){
    group_size = n/5
  }
  G = floor(p/group_size)
  #screen_stat = 0
  screen_stat_rep = matrix(NA, nrow=ncol(X),ncol=rep)
  for(i in 1:rep){
    #group = sample(1:G,p,replace=TRUE)
    group =kmeans(xgrid,centers = G)$cluster
    multi_res = fast.multi.screening.wald(Y,X,group=group,family=family)
    #selected_set = order(abs(multi_res$screen_stat),decreasing = TRUE)[1:n]
    #include_prob = include_prob + as.numeric(is.element(1:p,selected_set))
    #screen_stat = screen_stat + multi_res$screen_stat
    screen_stat_rep[,i] = multi_res$screen_stat
    #screen_stat = ifelse(abs(screen_stat)>abs(multi_res$screen_stat),screen_stat,multi_res$screen_stat)
    if(verbose==TRUE){
      cat(i,"\n")
      flush.console()
    }
  }
  #screen_stat = screen_stat/rep
  screen_stat = apply(screen_stat_rep,1,median)
  selected_order = selected.set(screen_stat,cut_num = p)
  return(list(screen_stat=screen_stat,selected_order=selected_order,
              elapsed = proc.time()[3]-elapsed))
}

fast.random.multi.screening = function(Y,X,family="gaussian",
                                       rep=20, group_size=NULL,
                                       verbose=FALSE,
                                       mc.cores=7){
  
  elapsed = proc.time()[3]
  p = ncol(X)
  n = nrow(X)
  
  if(is.null(group_size)){
    group_size = n/5
  }
  G = floor(p/group_size)
  screen_stat = 0
  for(i in 1:rep){
    group = sample(1:G,p,replace=TRUE)
    multi_res = fast.multi.screening(Y,X,group=group,
                                     family=family,
                                     mc.cores=mc.cores)
    screen_stat = screen_stat + multi_res$screen_stat
    if(verbose==TRUE){
      cat(i,"\n")
      flush.console()
    }
  }
  screen_stat = screen_stat/rep
  selected_order = selected.set(screen_stat,cut_num = p)
  return(list(screen_stat=screen_stat,selected_order=selected_order,
              elapsed = proc.time()[3]-elapsed))
}



cond.random.multi.screening = function(Y,X,family="gaussian",rep=20, group_size=NULL,Cset=NULL,verbose=FALSE){
  
  elapsed = proc.time()[3]
  p = ncol(X)
  n = nrow(X)
  
  if(is.null(group_size)){
    group_size = n/5
  }
  G = floor(p/group_size)
  #screen_stat = 0
  screen_stat_rep = matrix(NA, nrow=ncol(X),ncol=rep)
  for(i in 1:rep){
    group = sample(1:G,p,replace=TRUE)
    multi_res = cond.multi.screening(Y,X,Cset=Cset,group=group,family=family)
    #selected_set = order(abs(multi_res$screen_stat),decreasing = TRUE)[1:n]
    #include_prob = include_prob + as.numeric(is.element(1:p,selected_set))
    #screen_stat = screen_stat + multi_res$screen_stat
    screen_stat_rep[,i] = multi_res$screen_stat
    if(verbose==TRUE){
    cat(i,"\n")
    flush.console()
    }
  }
  #screen_stat = screen_stat/rep
  screen_stat = apply(screen_stat_rep,1,median)
  if(is.null(Cset))
    selected_order = selected.set(screen_stat,cut_num = p)
  else{
    temp = screen_stat
    temp[Cset] = 0
    selected_order = c(Cset,selected.set(temp,cut_num = p)[1:(p-length(Cset))])
  }
  
  return(list(screen_stat=screen_stat,selected_order=selected_order,
              elapsed = proc.time()[3]-elapsed))
}


cond.screening = function(y,x,Cset,family="gaussian"){
  elapsed = proc.time()[3]
  
  n=dim(x)[1];p=dim(x)[2]
  fullset=1:p;
  cSize=length(Cset);
  scanset=fullset[-Cset];
  
  ## standardize covariates
  x<- apply(x,2, function(x) (x-mean(x))/sd(x))
  ##store the abs(beta)
  beta=vector(mode = "numeric", p) 
  
  ## regress y on x_j given Cset. 
  for(j in scanset){
    Xj<- x[,j]
    Xtemp=cbind(x[,Cset],Xj);
    glm=glm(y~Xtemp,family=family)
    beta[j]=glm$coeff[cSize+2];
  }
  
  ## if CSIS for the cSize are 0, CSIS is the same as the SIS
  if(cSize==0){
    for(j in scanset){
      Xj<- x[,j]
      Xtemp=cbind(Xj);
      glm=glm(y~Xtemp,family=family)
      beta[j]=glm$coeff[2];
    }
  }
  
  screen_stat = beta
  if(is.null(Cset)){
    selected_order = selected.set(screen_stat,cut_num = p)
  }
  else{
    temp = screen_stat
    temp[Cset] = 0
    selected_order = c(Cset,selected.set(temp,cut_num = p)[1:(p-cSize)])
  }
  
  
  #beta[Cset]=NA
  return(list(screen_stat=screen_stat, selected_order=selected_order, 
              elapsed=proc.time()[3]-elapsed))
}

simul.example.1 = function(n,p,rho=0.5,R2 = 0.5){
  betacoef = c(rep(3,5),-7.5,rep(0,length=p-6))
  X = matrix(rnorm(n*p,mean=0,sd=sqrt(1-rho)),nrow=n,ncol=p) + matrix(rnorm(n,mean=0,sd=sqrt(rho)),nrow=n,ncol=p)
  mu = X%*%betacoef
  var_mu = var(mu)
  sigma2 = var_mu*(1/R2-1)
  Y = rnorm(n,mean=mu,sd=sqrt(sigma2))
  return(list(Y=Y,X=X,betacoef=betacoef,sigma2=sigma2,
              trueset=which(betacoef!=0)))
}

simul.example.1.bin = function(n,p,rho=0.5){
  betacoef = c(rep(3,5),-7.5,rep(0,length=p-6))
  X = matrix(rnorm(n*p,mean=0,sd=1),nrow=n,ncol=p) + matrix(rnorm(n,mean=0,sd=sqrt(rho)),nrow=n,ncol=p)
  mu = X%*%betacoef
  Y=as.numeric(1/(1+exp(-mu))>runif(n))
  return(list(Y=Y,X=X,betacoef=betacoef,
              trueset=which(betacoef!=0)))
}

simul.example.2 = function(n,p,q=20, rho=0.9,R2 = 0.9,signal=0.5){
  betacoef = c(signal*(-1)^(1:q),rep(0,length=p-q))
  #Auto correaltion
  dat=fast_simul_X_block_diag_AR1(n,p,blockSize = 100,rho = rho,sigma = 1)
  dat$betacoef = betacoef
  mu = dat$X%*%betacoef
  var_mu = var(mu[,1])
  sigma2 = var_mu*(1/R2-1)
  dat$Y = dat$X%*%betacoef+rnorm(n,sd=sqrt(sigma2))
  dat$trueset=which(betacoef!=0)
  return(dat)
}

simul.example.2.bin = function(n,p,q=10, rho=0.9,signal=0.5){
 #Auto correaltion binary 
  betacoef = c(signal*(-1)^(1:q),rep(0,length=p-q))
  #Auto correaltion
  dat=fast_simul_X_block_diag_AR1(n,p,blockSize = 100,rho = rho,sigma = 1)
  dat$betacoef = betacoef
  mu = dat$X%*%betacoef
  dat$Y=as.numeric(1/(1+exp(-mu[,1]))>runif(n))
  dat$trueset=which(betacoef!=0)
  return(dat)
}

simul.example.3 = function(n,p,q=20, rho=0.6,R2 = 0.9,signal=0.5){
  betacoef = c(signal*(-1)^(1:q),rep(0,length=p-q))
  #Auto correaltion
  dat=fast_simul_X_block_diag_comp_sym(n,p,blockSize = 100,rho = rho,sigma = 1)
  dat$betacoef = betacoef
  mu = dat$X%*%betacoef
  var_mu = var(mu[,1])
  sigma2 = var_mu*(1/R2-1)
  dat$Y = mu+rnorm(n,sd=sqrt(sigma2))
  dat$trueset=which(betacoef!=0)
  return(dat)
}

simul.example.3.bin = function(n,p,q=10, rho=0.9,signal=0.5){
  betacoef = c(signal*(-1)^(1:q),rep(0,length=p-q))
  #Auto correaltion
  dat=fast_simul_X_block_diag_comp_sym(n,p,blockSize = 100,rho = rho,sigma = 1)
  dat$betacoef = betacoef
  mu = dat$X%*%betacoef
  dat$Y=as.numeric(1/(1+exp(-mu[,1]))>runif(n))
  dat$trueset=which(betacoef!=0)
  return(dat)
}


simul.example.4 = function(n,p,rho=0.95,R2 = 0.9,signal=0.5,group=NULL){
  #spatial correaltion compound symetric
  m = floor(sqrt(p))
  #dat=simul_spat_dat_comp_sym(n,m,rho=rho,sigma=1)
  
  if(is.null(group)){
    group=kmeans(dat$xgrid,centers=100)$cluster
  }
  
  
  #dat = simul_spat_dat_group_RF(n,m,rho,sigma2 = 1,group=group)
  dat = simul_spat_dat_RF(n,m,rho,sigma2 = 1)
  #betacoef = ifelse(abs(dat$xgrid[,1]-0.1)+abs(dat$xgrid[,2]-0.1)<0.04,signal,0)
  betacoef = rep(0,length=p)
  betacoef[group==1] = signal*(-1)^(1:sum(group==1))
  dat$betacoef = betacoef
  mu = dat$X%*%betacoef
  var_mu = var(mu[,1])
  sigma2 = var_mu*(1/R2-1)
  dat$Y = dat$X%*%betacoef+rnorm(n,sd=sqrt(sigma2))
  dat$trueset=which(betacoef!=0)
  dat$delta = (betacoef!=0)
  dat$group = group
  dat$sigma2 = sigma2
  return(dat)
}


simul.example.4.bin = function(n,p,rho=0.95,signal=0.5,group=NULL){
  #spatial correaltion compound symetric
  m = floor(sqrt(p))
  #dat=simul_spat_dat_comp_sym(n,m,rho=rho,sigma=1)
  
  if(is.null(group)){
    group=kmeans(dat$xgrid,centers=100)$cluster
  }
  
  
  #dat = simul_spat_dat_group_RF(n,m,rho,sigma2 = 1,group=group)
  dat = simul_spat_dat_RF(n,m,rho,sigma2 = 1)
  #betacoef = ifelse(abs(dat$xgrid[,1]-0.1)+abs(dat$xgrid[,2]-0.1)<0.04,signal,0)
  betacoef = rep(0,length=p)
  betacoef[group==1] = signal*(-1)^(1:sum(group==1))
  dat$betacoef = betacoef
  mu = dat$X%*%betacoef
  dat$Y=as.numeric(1/(1+exp(-mu[,1]))>runif(n))
  dat$trueset=which(betacoef!=0)
  
  return(dat)
}




simul.example.1.a = function(n,p,rho=0.5,R2 = 0.9){
  betacoef = c(10,rep(0,length=p-2),1)
  X = matrix(rnorm(n*(p-1),mean=0,sd=1),nrow=n,ncol=p-1) + matrix(rnorm(n,mean=0,sd=sqrt(rho)),nrow=n,ncol=p-1)
  X = cbind(X,rnorm(n))
  mu = X%*%betacoef
  var_mu = var(mu)
  sigma2 = var_mu*(1/R2-1)
  Y = rnorm(n,mean=mu,sd=sqrt(sigma2))
  return(list(Y=Y,X=X,betacoef=betacoef,sigma2=sigma2,
              trueset=which(betacoef!=0)))
}

simul.example.1.Barut_Ex3 = function(n,p,q=12,s0=100,rho=0.5,R2 = 0.9){
  betacoef = c(rep(c(1,1.3),floor(q/2)),rep(0,p-2*floor(q/2)))
  #print(length(betacoef))
  E1 = matrix(rnorm(n*floor(p/3),mean=0,sd=1),nrow=n,ncol=floor(p/3))
  #E0 = 
  E2 = matrix(sample(c(-1,1),n*floor(p/3),replace=TRUE)*rexp(n*floor(p/3),1),nrow=n,ncol=floor(p/3))
  p3 = p-2*floor(p/3)
  p31 = floor(p3/2)
  p32 = p3 - p31
  temp = c(rnorm(n*p31,mean=-1,sd=1),rnorm(n*p32,mean=1,sd=0.5))
  E3 = matrix(temp[sample(1:(n*p3),p3*n)],nrow=n,ncol=p3)
  a2 = rho/(1-rho)
  a = sqrt(a2)
  E0 = cbind(a*matrix(rnorm(n,mean=0,sd=1),nrow=n,ncol=s0),matrix(0,nrow=n,ncol=p-s0))
  #X = matrix(rnorm(n*(p-1),mean=0,sd=1),nrow=n,ncol=p-1) + matrix(rnorm(n,mean=0,sd=sqrt(rho)),nrow=n,ncol=p-1)
  
  X = cbind(E1,E2,E3)+E0
  X[,1:s0] = X[,1:s0]/sqrt(1+a2)
  
  X=sapply(1:ncol(X),function(i) normalize(X[,i]))
  #print(dim(X))
  mu = X%*%betacoef
  var_mu = var(mu)
  sigma2 = var_mu*(1/R2-1)
  Y = rnorm(n,mean=mu,sd=sqrt(sigma2))
  return(list(Y=Y,X=X,betacoef=betacoef,sigma2=sigma2,
              trueset=which(betacoef!=0)))
}

simul.new.example2 = function(n,p,q=2,signal=0.1,rho=0.5,R2 = 0.9,sigma2=1.0){
  betacoef = c(signal*rep(c(1,-1),floor(q/2)),rep(0,p-q))
  #print(length(betacoef))
  
  X =sqrt(rho)*matrix(rnorm(n*p),nrow=n,ncol=p) + sqrt(1-rho)*matrix(rnorm(n),nrow=n,ncol=p)
  
  
  #print(dim(X))
  mu = X%*%betacoef
  var_mu = var(mu)
  if(is.null(sigma2)){
    sigma2 = var_mu*(1/R2-1)
  }
  Y = rnorm(n,mean=mu,sd=sqrt(sigma2))
  return(list(Y=Y,X=X,betacoef=betacoef,sigma2=sigma2,
              trueset=which(betacoef!=0)))
}


simul.new.example3 = function(n,p,q=2,rho=0.5,R2 = 0.9,sigma2=1.0){
  betacoef = c(rep(c(0.5,-1),floor(q/2)),rep(0,p-q))
  #print(length(betacoef))
  
  X =sqrt(rho)*matrix(rnorm(n*p),nrow=n,ncol=p) + sqrt(1-rho)*matrix(rnorm(n),nrow=n,ncol=p)
  
  
  #print(dim(X))
  mu = X%*%betacoef
  var_mu = var(mu)
  if(is.null(sigma2)){
    sigma2 = var_mu*(1/R2-1)
  }
  Y = rnorm(n,mean=mu,sd=sqrt(sigma2))
  return(list(Y=Y,X=X,betacoef=betacoef,sigma2=sigma2,
              trueset=which(betacoef!=0)))
}

simul.new.example4 = function(n,p,q=2,rho=0.5,R2 = 0.9,sigma2=1.0){
  betacoef = c(rep(c(0.5,0.5),floor(q/2)),rep(0,p-q))
  #print(length(betacoef))
  
  X =sqrt(rho)*matrix(rnorm(n*p),nrow=n,ncol=p) + sqrt(1-rho)*matrix(rnorm(n),nrow=n,ncol=p)
  
  
  #print(dim(X))
  mu = X%*%betacoef
  var_mu = var(mu)
  if(is.null(sigma2)){
    sigma2 = var_mu*(1/R2-1)
  }
  Y = rnorm(n,mean=mu,sd=sqrt(sigma2))
  return(list(Y=Y,X=X,betacoef=betacoef,sigma2=sigma2,
              trueset=which(betacoef!=0)))
}

# simul.example.2 = function(n,p,blockSize=10,connectProb=0.1,connectSignal=1,R2=0.5){
#   betacoef = c(rep(3,5),-7.5,rep(0,length=p-6))
#   temp = simul_X_block_diag_random_graph(n,p,blockSize=blockSize,connectProb=connectProb,
#                                          connectSignal=connectSignal)
#   X = temp$X
#   mu = X%*%betacoef
#   var_mu = var(mu)
#   sigma2 = var_mu*(1/R2-1)
#   Y = rnorm(n,mean=mu,sd=sqrt(sigma2))
#   return(list(Y=Y,X=X,betacoef=betacoef,sigma2=sigma2,
#               trueset=which(betacoef!=0),covMat = temp$covMat,od_idx = temp$od_idx))
# }

simul.new.example = function(n,p,q=2,signal=0.1,rho=0.5,R2 = 0.9,sigma2=1.0){
  betacoef = c(signal*rep(c(1,-1),floor(q/2)),rep(0,p-q))
  #print(length(betacoef))
  
  X12 =sqrt(rho)*matrix(rnorm(n*q),nrow=n,ncol=q) + sqrt(1-rho)*matrix(rnorm(n),nrow=n,ncol=q)
  X= cbind(X12,matrix(rnorm(n*(p-q)),nrow=n,ncol=p-q))
  
  #print(dim(X))
  mu = X%*%betacoef
  var_mu = var(mu)
  if(is.null(sigma2)){
    sigma2 = var_mu*(1/R2-1)
  }
  Y = rnorm(n,mean=mu,sd=sqrt(sigma2))
  return(list(Y=Y,X=X,betacoef=betacoef,sigma2=sigma2,
              trueset=which(betacoef!=0)))
}



selected.set = function(screen_stat,cut_num){
  return(order(abs(screen_stat),decreasing = TRUE)[1:cut_num])
}



summary.screening = function(res,trueset, cut_num,method=""){
  # res should be a list variable with 
  # at least two components
  # 1. selected_order
  # 2. elapsed
  
  true_positives = sum(is.element(trueset,res$selected_order[1:cut_num]))
  model_size = max(which(is.element(res$selected_order,trueset)))
  sum_res = list(method=method,true_positives=true_positives,model_size=model_size,computing_time=round(res$elapsed,digits=3))
  return(sum_res)
}

summary.screening.pit = function(res,trueset, cut_num,method=""){
  # res should be a list variable with 
  # at least two components
  # 1. selected_order
  # 2. elapsed
  
  true_positives = sum(is.element(trueset,res$selected_order[1:cut_num]))
  model_size = max(which(is.element(res$selected_order,trueset)))
  sum_res = list(method=method,true_positives=true_positives,model_size=model_size,computing_time=round(res$elapsed,digits=3))
  return(sum_res)
}

prior.random.screening = function(Y,X,family="gaussian",rep=20, group_size=NULL,
                                  Cset=NULL,prior_weight=0.5,verbose=FALSE){
  elapsed = proc.time()[3]
  p = ncol(X)
  n = nrow(X)
  
  if(is.null(group_size)){
    group_size = n/5
  }
  G = floor(p/group_size)
  screen_stat = 0
  for(i in 1:rep){
    group = sample(1:G,p,replace=TRUE)
    multi_res = prior.cond.multi.screening(Y,X,Cset=Cset,group=group,family=family,
                                           prior_weight=prior_weight)
    
    #include_prob = include_prob + as.numeric(is.element(1:p,selected_set))
    screen_stat = screen_stat + multi_res$screen_stat
    if(verbose==TRUE){
    cat(i,"\n")
    flush.console()
    }
  }
  screen_stat = screen_stat/rep
  if(is.null(Cset))
    selected_order = selected.set(screen_stat,cut_num = p)
  else{
    temp = screen_stat
    temp[Cset] = 0
    selected_order = c(Cset,selected.set(temp,cut_num = p)[1:(p-length(Cset))])
  }
  
  return(list(screen_stat=screen_stat,selected_order=selected_order,
              elapsed = proc.time()[3]-elapsed))
}

pretty.print = function(list_var){
  return(noquote(format(list_var)))
}


fast.cov = function(X){
  X_bar = apply(X,2,mean)
  ones = matrix(1,nrow=nrow(X),ncol=1)
  X_center = X-ones%*%X_bar
  cov_mat = t(X_center)%*%X_center/(nrow(X)-1)
  return(cov_mat)
}

summary.sim.tab = function(alltab,methods,measures,true_num){
  med_tab = apply(alltab,c(1,2),median)
  rownames(med_tab) = methods
  colnames(med_tab) = measures
  avg_tab = apply(alltab,c(1,2),mean)
  rownames(avg_tab) = methods
  colnames(avg_tab) = measures
  std_tab = apply(alltab,c(1,2),sd)
  rownames(std_tab) = methods
  colnames(std_tab) = measures
  ucl_tab = apply(alltab,c(1,2),quantile,prob=0.75)
  lcl_tab = apply(alltab,c(1,2),quantile,prob=0.25)
  iqr_tab = ucl_tab-lcl_tab
  rownames(iqr_tab) = methods
  colnames(iqr_tab) = measures
  
  med_print_tab = matrix(paste(med_tab,"(",round(iqr_tab,digits=1),")",sep=""),
                     nrow=nrow(med_tab),ncol=ncol(med_tab))
  
  rownames(med_print_tab) = methods
  colnames(med_print_tab) = measures
  
  
  avg_print_tab = matrix(paste(round(avg_tab,digits=1),"(",round(std_tab,digits=1),")",sep=""),
                         nrow=nrow(avg_tab),ncol=ncol(avg_tab))
  
  rownames(avg_print_tab) = methods
  colnames(avg_print_tab) = measures
  
  include_prob =apply(alltab[,1,]==true_num,1,sum)/dim(alltab)[3]
  
  all_print_tab = cbind(round(include_prob,digits=2),
                        avg_print_tab[,1],med_print_tab[,1],
                        avg_print_tab[,2],med_print_tab[,2],
                        avg_print_tab[,3],med_print_tab[,3])
  colnames(all_print_tab) = c("Inclusion Prob","TP Mean (Sd)","TP Median (Iqr)",
                              "MMS Mean (Sd)","MMS Median (Iqr)",
                              "Comp.Time Mean (Sd)","Comp.Time Median (Iqr)")
  
  return(list(med_tab=med_tab,iqr_tab=iqr_tab,
              avg_tab=avg_tab,std_tab=std_tab,
              med_print_tab=med_print_tab,
              avg_print_tab=avg_print_tab,
              all_print_tab=all_print_tab))
}

summary.sim.tab.simple = function(alltab,methods,measures,true_num){
  med_tab = apply(alltab,c(1,2),median)
  rownames(med_tab) = methods
  colnames(med_tab) = measures
  avg_tab = apply(alltab,c(1,2),mean)
  rownames(avg_tab) = methods
  colnames(avg_tab) = measures
  std_tab = apply(alltab,c(1,2),sd)
  rownames(std_tab) = methods
  colnames(std_tab) = measures
  ucl_tab = apply(alltab,c(1,2),quantile,prob=0.75)
  lcl_tab = apply(alltab,c(1,2),quantile,prob=0.25)
  iqr_tab = ucl_tab-lcl_tab
  rownames(iqr_tab) = methods
  colnames(iqr_tab) = measures
  
  med_print_tab = matrix(paste(med_tab,"(",round(iqr_tab,digits=1),")",sep=""),
                         nrow=nrow(med_tab),ncol=ncol(med_tab))
  
  rownames(med_print_tab) = methods
  colnames(med_print_tab) = measures
  
  
  avg_print_tab = matrix(paste(round(avg_tab,digits=1),"(",round(std_tab,digits=1),")",sep=""),
                         nrow=nrow(avg_tab),ncol=ncol(avg_tab))
  
  rownames(avg_print_tab) = methods
  colnames(avg_print_tab) = measures
  
  include_prob =apply(alltab[,1,]==true_num,1,sum)/dim(alltab)[3]
  
  all_print_tab = cbind(round(include_prob,digits=2),
                        avg_print_tab[,1],med_print_tab[,1],
                        avg_print_tab[,2],med_print_tab[,2],
                        avg_print_tab[,3],med_print_tab[,3])
  colnames(all_print_tab) = c("Inclusion Prob","TP Mean (Sd)","TP Median (Iqr)",
                              "MMS Mean (Sd)","MMS Median (Iqr)",
                              "Comp.Time Mean (Sd)","Comp.Time Median (Iqr)")
  
  return(list(med_tab=med_tab,iqr_tab=iqr_tab,
              avg_tab=avg_tab,std_tab=std_tab,
              med_print_tab=med_print_tab,
              avg_print_tab=avg_print_tab,
              all_print_tab=all_print_tab))
}


find.partition.hclust = function(X,G=NULL,dist_method="manhattan"){
  n = nrow(X)
  p = ncol(X)
  if(is.null(G)){
    group_size = max(floor(sqrt(n)),1)
    G = floor(p/group_size)
  }
  cor_mat = NULL
  if(dist_method=="correlation"){
    cor_mat = fast_cor_mat(X)
    dist_mat = 1-abs(cor_mat)
    dist_vec = as.dist(dist_mat)
  }
  else{
    dist_vec = dist(t(X),method=dist_method)
  }
  hc = hclust(dist_vec)
  group = cutree(hc,k=G)
  return(list(group=group,hc=hc,G=G,dist_vec=dist_vec,cor_mat=cor_mat))
}


image.cor.mat = function(cor_mat,sub_sample=0.1,cols=my.cols){
  p = nrow(cor_mat)
  q = floor(sub_sample*p)+1
  shift = ceiling(p/q)
  plotdat = list()
  sub_idx = seq(1,p-shift+1,by=shift)
  plotdat$z = c(cor_mat[sub_idx,sub_idx])
  if(shift>1){
    for(i in 1:(shift-1)){
      plotdat$z = plotdat$z + c(cor_mat[sub_idx+i,sub_idx+i])
    }
  }
  plotdat$z = plotdat$z/shift
  
  temp = expand.grid(sub_idx,sub_idx)
  plotdat$x = temp[,1]
  plotdat$y = temp[,2]
  levelplot(z~x+y,data=plotdat,aspect="iso",col.regions=cols,at=seq(-1.01,1.01,length=length(cols)-1))
}


spatial_kernel_smooth = function(val,grids,neighbors=NULL,radius=3,n_cor=0.9){
  if(is.null(neighbors)){
    neightbors=find_neighbors(grids,radius)
  }
  rho = -log(n_cor)/sqrt(sum((grids[1,]-grids[2,])^2))
  smooth_val = mclapply(1:length(val),function(i){
    weighted.mean(val[neighbors[[i]][,1]],w=exp(-rho*neighbors[[i]][,2]))
  },mc.cores=7)
  return(unlist(smooth_val))
}


find_neighbors_2D = function(grids,radius = 0.05){
  names(grids) = c("x","y")
  x=unique(grids$x)
  y=unique(grids$y)
  x_neighbors = mclapply(x,function(a) which(abs(grids$x-a)<radius),mc.cores = 7)
  y_neighbors = mclapply(y,function(a) which(abs(grids$y-a)<radius),mc.cores = 7)

  
  names(x_neighbors) = paste(x)
  names(y_neighbors) = paste(y)
  
  neighbors=mclapply(1:nrow(grids),function(i){
    idx = intersect(x_neighbors[[paste(grids$x[i])]],y_neighbors[[paste(grids$y[i])]])
    distance = sqrt((grids$x[idx]-grids$x[i])^2+(grids$y[idx]-grids$y[i])^2)
    return(cbind(idx,distance))
    #temp[which(distance<radius)]
  },mc.cores = 7)
  return(neighbors)
}

show_image = function(img,all_coords,zlist,my.cols,col_lim = c(-0.5,0.8),
                      layout=NULL,bgcol="white",main=""){
  
  idx =which(is.element(all_coords$z,zlist))
  z = sprintf("z = %2dmm",3*(all_coords$z[idx]-25))
  fig=levelplot(img[idx]~all_coords$x[idx]+all_coords$y[idx] | z, at=seq(col_lim[1],col_lim[2],length=length(my.cols)-1),
                col.regions = my.cols,cuts = length(my.cols),
                par.settings=list(panel.background=list(col =bgcol),grid.pars = list(fontfamily = 'Times')), 
                xlab="x",ylab="y",
                layout=layout,aspect=1.0,main=main)
  return(fig)
}


merge.group = function(group){
  unique_group = unique(group)
  unique_group= setdiff(unique_group,0)
  region_coords = lapply(unique_group, function(i) which(group==i,arr.ind=TRUE))
  names(region_coords) = paste(unique_group)
  region_center = t(sapply(unique_group,function(i) apply(region_coords[[paste(i)]],2,mean)))
  dist_vec=dist(region_center)
  # hc=hclust(dist_vec)
  # memb=cutree(hc,k=512)
  dist_mat = as.matrix(dist_vec)
  idx_set = 1:nrow(dist_mat)
  merge_set = NULL
  i = idx_set[1]
  while(length(idx_set)>=1){
    target_idx = setdiff(idx_set,i)
    j = target_idx[which.min(dist_mat[i,target_idx])]
    idx_set = setdiff(idx_set, c(i,j))
    merge_set = rbind(merge_set,c(i,j))
    i = idx_set[1]
  }
  
  group_new = group
  for(i in 1:nrow(merge_set)){
    group_new[group==merge_set[i,2]] = merge_set[i,1]
  }
  return(group_new)
  
}


cor.merge.group = function(group,X,pix_index){
  unique_group = unique(group)
  unique_group = setdiff(unique_group,0)
  unique_group = sort(unique_group)
  temp=aggregate(t(X),by=list(group[pix_index]),FUN=mean)
  cor_mat = fast_cor_mat(t(temp[,1:nrow(X)+1]))
  dist_mat = 1 - abs(cor_mat)
  #dist_vec=dist(region_center)
  dist_vec = as.dist(dist_mat)
  # hc=hclust(dist_vec)
  # memb=cutree(hc,k=512)
  #dist_mat = as.matrix(dist_vec)
  idx_set = 1:nrow(dist_mat)
  merge_set = NULL
  i = idx_set[1]
  while(length(idx_set)>=1){
    target_idx = setdiff(idx_set,i)
    j = target_idx[which.min(dist_mat[i,target_idx])]
    idx_set = setdiff(idx_set, c(i,j))
    merge_set = rbind(merge_set,c(i,j))
    i = idx_set[1]
  }
  
  group_new = group
  for(i in 1:nrow(merge_set)){
    group_new[group==unique_group[merge_set[i,2]]] = unique_group[merge_set[i,1]]
  }
  return(group_new)
  
}

print.mean.sd.tab = function(mean_tab,sd_tab,format="%.2f (%.2f)"){
  p = ncol(mean_tab)
  print_tab = NULL
  for(i in 1:p)
    print_tab = cbind(print_tab,sprintf(format,mean_tab[,i],sd_tab[,i]))
  
  return(print_tab)
}

print.mean.tab = function(mean_tab,format="%.2f"){
  p = ncol(mean_tab)
  print_tab = NULL
  for(i in 1:p)
    print_tab = cbind(print_tab,sprintf(format,mean_tab[,i]))
  
  return(print_tab)
}

compute.summary.stat = function(tab,re_order,truelen=6){
  rtab = tab[re_order,]
  TPR = rtab[,1]/truelen
  PIT = rtab[,1]==truelen
  MMS = rtab[,2]
  Time = rtab[,3]
  outtab = cbind(PIT,MMS,TPR,Time)
  return(outtab)
}

compute.summary.stat.simple = function(tab,re_order,truelen=6){
  rtab = tab[re_order,]
  TPR = rtab[,1]/truelen
  PIT = rtab[,1]==truelen
  MMS = rtab[,2]
  outtab = cbind(PIT,MMS,TPR)
  return(outtab)
}



sum.med.IQR = function(outtab,method_names){
  medtab = apply(outtab,1:2,median)
  IQRtab = apply(outtab,1:2,IQR)
  meantab = apply(outtab,1:2,mean)
  sdtab = apply(outtab,1:2,sd)
  med_print_tab = print.mean.sd.tab(medtab,IQRtab)
  mean_print_tab = print.mean.sd.tab(meantab,sdtab)
  print_tab = cbind(mean_print_tab[,1],med_print_tab[,2:3],mean_print_tab[,4])
  colnames(print_tab) = colnames(medtab)
  rownames(print_tab) = method_names
  return(list(print_tab = print_tab,
              med_print_tab=med_print_tab,
              mean_print_tab=mean_print_tab,
              medtab=medtab,IQRtab=IQRtab))
}


sum.med = function(outtab,method_names){
  medtab = apply(outtab,1:2,median)
  meantab = apply(outtab,1:2,mean)
  med_print_tab = print.mean.tab(medtab)
  mean_print_tab = print.mean.tab(meantab)
  print_tab = cbind(mean_print_tab[,1],med_print_tab[,2:3])
  colnames(print_tab) = colnames(medtab)
  rownames(print_tab) = method_names
  return(list(print_tab = print_tab,
              med_print_tab=med_print_tab,
              mean_print_tab=mean_print_tab,
              medtab=medtab,IQRtab=IQRtab))
}


draw.spat.part = function(group,xgrid,graycols = gray(c(1,0.5)),
                          highlight=c(1,24),hcols="yellow",main="",hgrp=NULL){
  m = sqrt(length(group))
  group_mat = matrix(group,nrow=m,ncol=m)
  diff_mat_1 = sapply(1:m,function(i) abs(diff(group_mat[,i]))>0)
  diff_mat_2 = t(sapply(1:m,function(i) abs(diff(group_mat[i,]))>0))
  diff_mat = rbind(diff_mat_1,rep(FALSE,m)) | cbind(diff_mat_2,rep(rep(FALSE,m)))
  diff_mat[diff_mat==0] = NA
  group_center = aggregate(xgrid,by=list(group),FUN=mean)
  
  my.panel.fun = function(x,y,z,...){
    panel.levelplot(x,y,z,...)
    panel.text(group_center[,2],group_center[,3],
               labels=group_center[,1],cex=0.8)
  }
  
  fig_part_00 = levelplot(diff_mat~dat$xgrid[,1]+dat$xgrid[,2],
                          col.regions = graycols,cuts = length(graycols),
                          aspect = 1.0,xlab="",ylab="",
                          colorkey=FALSE,
                          scales=list(draw=FALSE),panel = my.panel.fun)
  
  if(is.null(hgrp)){
    hgrp = rep(NA, length=length(group))
    hgrp[is.element(group,highlight)] = 1
  }
  
  fig_highlight = levelplot(hgrp~dat$xgrid[,1]+dat$xgrid[,2],
                          col.regions = hcols,cuts = length(hcols),
                          aspect = 1.0,xlab="",ylab="",
                          colorkey=FALSE,
                          scales=list(draw=FALSE),main=main,
                          par.settings = list(panel.background=list(col="lightblue",alpha=0.5),par.main.text=list(font=1,just="left")))
  
  return(fig_highlight+as.layer(fig_part_00))
}


holp_group_partS = function(Y,X,Z,group,family="gaussian",cutoff=2700){
  holp_res = self.holp(Y=Y,X=X,Z=Z,family=family,lambda=lambda)
  holp_selected_idx = holp_res$selected_order[1:cutoff]
  group_1  = group
  group_1[-holp_selected_idx] = group[-holp_selected_idx]+max(group)
  #fast.multi.screening.glmnet(Y)
  return(group_1)
}

show.selected.region = function(selected_idx,AALregion,region_tab,count_cutoff = 40){
  selected_region_count = sort(table(AALregion[pix_index[selected_idx]]),decreasing = TRUE)
  selected_region_index = as.numeric(names(selected_region_count))
  subidx = which(selected_region_count>=count_cutoff)
  selected_region =  region_tab[selected_region_index[subidx],3]
  
  top_rank = sapply(1:length(subidx),function(i) min(which(AALregion[pix_index[selected_idx]]==selected_region_index[subidx[i]])))
  output = data.frame(selected_region = as.character(selected_region),
                      number_voxels = selected_region_count[subidx],
                      top_rank = top_rank)
  return(output)
}

hist.two = function(x,y,main="", xlab="Values",legends=c("Current Method","New Method"),
                    cex=3,mar=c(3,3,3,3),leg.loc="top",breaks=10,freq=FALSE){
  xlim = c(min(x,y),max(x,y))
  a=hist(x,plot=FALSE,breaks=breaks)
  b=hist(y,plot=FALSE,breaks=1.5*breaks)
  if(freq)
    ylim = c(min(a$counts,b$counts),max(a$counts,b$counts)*1.2)
  else
    ylim = c(min(a$density,b$density),max(a$density,b$density)*1.2)
  par(mar=mar)
  plot(a,col=rgb(0,0,1,1/4),xlim=xlim,ylim=ylim,main=main,xlab=xlab,cex.main=cex,cex.axis=cex,cex.lab=cex,freq=freq)
  plot(b,col=rgb(1,0,0,1/4),add=TRUE,cex.main=cex,cex.axis=cex,cex.lab=cex,freq=freq)
  legend(leg.loc,legend=legends,
         pch=15,col=c(rgb(0,0,1,1/4),rgb(1,0,0,1/4)),cex=cex)
  box()
  
}


hist.three = function(x,y,z,main="", xlab="Values",legends=c("hist1","hist2","hist3"),
                    cex=1,mar=c(3,3,3,3),leg.loc="top",breaks=10*rep(1,3),freq=FALSE,xlim=NULL){
  if(is.null(xlim)){
    xlim = c(min(x,y,z),max(x,y,z))
  }
  a=hist(x,plot=FALSE,breaks=breaks[1])
  b=hist(y,plot=FALSE,breaks=breaks[2])
  d=hist(z,plot=FALSE,breaks=breaks[3])
  if(freq)
    ylim = c(min(a$counts,b$counts,d$counts),max(a$counts,b$counts,d$counts)*1.2)
  else
    ylim = c(min(a$density,b$density,d$density),max(a$density,b$density,d$density)*1.2)
  par(mar=mar)
  plot(a,col=rgb(0,0,1,1/4),xlim=xlim,ylim=ylim,main=main,xlab=xlab,cex.main=cex,cex.axis=cex,cex.lab=cex,freq=freq)
  plot(b,col=rgb(1,0,0,1/4),add=TRUE,cex.main=cex,cex.axis=cex,cex.lab=cex,freq=freq)
  plot(d,col=rgb(0,1,0,1/4),add=TRUE,cex.main=cex,cex.axis=cex,cex.lab=cex,freq=freq)
  if(!is.null(legends)){
  legend(leg.loc,legend=legends,
         pch=15,col=c(rgb(0,0,1,1/4),rgb(1,0,0,1/4),rgb(0,1,0,1/4)),cex=cex)
  }
  box()
  
}

generate_group_list = function(G,group_size,indices,overlap=FALSE){
  group_list = list()
  if(overlap){
    for(g in 1:(G-1)){
      group_list[[g]] = sample(indices,group_size)
    }
    all_indices = unique(unlist(group_list))
    if(length(all_indices)==length(indices)){
      group_list[[G]] = sample(indices,group_size)
    }
    else{
      group_list[[G]] = setdiff(indices,all_indices)
    }
  }
  else{
    p = length(indices)
    group=sample(1:G,p,replace=TRUE)
    for(g in 1:G){
      group_list[[g]] = indices[which(group==g)]
    }
    
  }
  return(group_list)
}

fit_nbc = function(Y,x){
  idx1=which(Y==1)
  idx0=which(Y==0)
  mu1 = apply(x[idx1,],2,mean)
  mu0 = apply(x[idx0,],2,mean)
  sd1 = apply(x[idx1,],2,sd)
  sd0 = apply(x[idx0,],2,sd)
  return(list(mu1=mu1,mu0=mu0,sd1=sd1,sd0=sd0))
}

predict_nbc= function(newx,paras){
  log_odds = sapply_pb(1:nrow(newx), function(i){
    as.numeric(sum((log(paras$sd0)-log(paras$sd1))+0.5*(((newx[i,]-paras$mu0)/paras$sd0)^2-((newx[i,]-paras$mu1)/paras$sd1)^2)))
  })
  return(1/(1+exp(-log_odds)))
}


cv_roc = function(Y,Ypred_test,test_idx_list){
  K = length(test_idx_list)
  holp_Ypred_test = rep(NA,length=length(Y))
  for(k in 1:K)
    holp_Ypred_test[test_idx_list[[k]]] = Ypred_test[[paste(k)]]
  
  holp_thres = c(-Inf,sort(unique(holp_Ypred_test)),Inf)
  
  holp_roc = sapply(1:length(holp_thres), function(i) {
    tab = table(factor(as.numeric(holp_Ypred_test > holp_thres[i]),levels=c(1,0)),factor(Y,levels=c(1,0)))
    return(c(tab[1,2]/(tab[2,2]+tab[1,2]),tab[1,1]/(tab[2,1]+tab[1,1])))
  })
  od = order(holp_roc[1,])
  fpr = holp_roc[1,od]
  sens = holp_roc[2,od]
  k1 = 1
  k2 = 2
  idx = k1
  while(k2<=ncol(holp_roc)){
    if(fpr[k2]>=fpr[k1] && sens[k2]>=sens[k1]){
      idx = c(idx, k2)
      k1 = k2
    }
    k2 = k2 + 1
  }
  
  fpr = fpr[idx]
  sens = sens[idx]
  return(list(fpr=fpr,sens=sens,auc=sum(0.5*diff(fpr)*(sens[-length(sens)]+sens[-1]))))
}


PCA.correction = function(X,alpha=0.8){
  n = nrow(X)
  p = ncol(X)
  if(n<p){
    XtX = X%*%t(X)
    temp = eigen(XtX)
    eigen_vectors = t(X)%*%temp$vectors/matrix(sqrt(temp$values),byrow=TRUE,ncol=n,nrow=p)
    eigen_values = temp$values
  } else{
    tXX = t(X)%*%X
    temp1 = eigen(tXX)
    eigen_vectors = temp1$vectors
    eigen_values = temp$values
  }
    
  cum_sum = cumsum(eigen_values)
  ratio = cum_sum/cum_sum[length(cum_sum)]
  idx_cut = which.max(ratio>alpha)
  
  fit = lm(t(X)~0+eigen_vectors[,1:idx_cut])
  Z = t(coef(fit))
  newX = t(residuals(fit))
  return(list(Z=Z,X=newX))
}