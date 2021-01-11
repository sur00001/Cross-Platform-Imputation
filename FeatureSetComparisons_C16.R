###################################################################################
#This code compares all the feature sets - unfiltered sites (all sets in the chromosome),
#the most correlated sites (ex: top 10%) and PCs
#This is for cluster 16 (you can specify any cluster)
#This script only runs tuned penalized regressions
#######################################################################################

#libraries
library(glmnet)
library(data.table)
library(purrr)

####################### Functions ############################
#Tuned penalized regression
penregT<-function(i){
  print(i)
  #Store tuned parameters and coefficients for each regression
  coeff=c(); T.alphas=c(); T.lambdas=c(); y_pred=c()
  if(!any(is.na(y_train[,i]))){
    alpha=-1;cvm = 10000000;lambda = -1
    for(a in seq(0, 1, 0.1)) {
      cv <- cv.glmnet(x_train, y_train[,i], alpha=a, type.measure = "mse", nfolds = 10)
      if(min(cv$cvm)<cvm) {
        cvm <- min(cv$cvm)
        alpha <- a
        lambda <- cv$lambda.min
        fit=glmnet(x_train,y_train[,i],alpha=a,lambda=cv$lambda.min)
        y_pred=predict(fit,newx=x_test); rownames(y_pred)=NULL;
        y_pred=y_pred[,1]
        coeff=coef(fit)[,1]
      }}
    T.alphas=alpha
    T.lambdas=lambda
    return(list(al=T.alphas,lam=T.lambdas,coef=coeff,y_pred=y_pred))}}

#Fixed penalized regression (specify alpha and lambda)
a=0; l=.01
penregf<-function(i){
  #Store tuned parameters and coefficients for each regression
  coeff=c();y_pred=c()
  fit=glmnet(x_train,y_train[,i],alpha=a,lambda=l)
  y_pred=predict(fit,newx=x_test); rownames(y_pred)=NULL
  y_pred=y_pred[,1] #Want the y_pred vector only 
  coeff=coef(fit)[,1]
  return(list(coef=coeff,y_pred=y_pred))}

#Function to calculate pairwise correlations
cor.probe=function(i){
  cor.l=unlist(
    lapply(1:length(x_train[1,]),function(j){abs(cor(x_train[,j],y_train[,i]))}))
  return(cor.l)}

#Function for tuned penalized regression using top correlated sites
penregC<-function(i){
  print(i)
  #Store tuned parameters and coefficients for each regression
  coeff=c(); T.alphas=c(); T.lambdas=c(); y_pred=c()
  if(!any(is.na(y_train[,i]))){
    alpha=-1;cvm = 10000000;lambda = -1
    top.cors=cor.probe(i) #top correlations
    sites=which(top.cors>=quantile(top.cors, (1-Tquan))) #top correlated sites
    for(a in seq(0, 1, 0.1)) {
      cv <- cv.glmnet(x_train[,sites], y_train[,i], alpha=a, type.measure = "mse", nfolds = 10)
      if(min(cv$cvm)<cvm) {
        cvm <- min(cv$cvm)
        alpha <- a
        lambda <- cv$lambda.min
        fit=glmnet(x_train[,sites],y_train[,i],alpha=a,lambda=cv$lambda.min)
        y_pred=predict(fit,newx=x_test[,sites]); rownames(y_pred)=NULL;y_pred=y_pred[,1]
        coeff=coef(fit)[,1]}}
    T.alphas=alpha;T.lambdas=lambda
    return(list(al=T.alphas,lam=T.lambdas,coef=coeff,
                y_pred=y_pred,tcors=top.cors,x.sites=sites))}}


#Tuned regressions for different feature sets 
################## Unfiltered features ###################

#Open x_test and x_train datasets of unfiltered sites
load("C:/Users/Aparajita/Documents/Unf.testtrain22.rda")
x_test=Unf.x_test
x_train=Unf.x_train

#Open y_test and y_train datasets from cluster 16 
load("C:/Users/Aparajita/Documents/PScluster16.22.rda")
y_test=c16.y_test
y_train=c16.y_train

#Run regression which uses unfiltered feature sites
start_time <- Sys.time()
t.c16 =lapply(1:ncol(y_test),penregT) #tuned regression for cluster 16
end_time=Sys.time()
t=end_time-start_time;t

#Save results
save(glm.y_test=y_test,glm.y_train=y_train,glm.x_train=x_train,glm.x_test=x_test,
     reg.list= t.c16,sim.time=t,file=paste("c16Unf.rda"))

##################### Top correlated feature sites ######################

#Open x_test and x_train datasets
load("C:/Users/Aparajita/Documents/Unf.testtrain22.rda")
x_test=Unf.x_test
x_train=Unf.x_train

#Open y_test and y_train datasets from cluster 16 
load("C:/Users/Aparajita/Documents/PScluster16.22.rda")
y_test=c16.y_test
y_train=c16.y_train

Tquan=.20 #Specified 20% as the top quantile 

#Run regression which uses the top correlated unfiltered sites
start_time <- Sys.time()
cor20.c16 =lapply(1:ncol(y_test),penregC)
end_time=Sys.time()
t=end_time-start_time;t

#Save results
save(glm.y_test=y_test,glm.y_train=y_train,glm.x_train=x_train,glm.x_test=x_test,
     reg.list= cor20.c16,sim.time=t,Tquantile=Tquan,file=paste("c16Cor20.22.rda"))

######################## Sites with PCs from PCA ###############################

#Open x_test and x_train datasets that include PCs 
load("C:/Users/Aparajita/Documents/UnfilteredTPC.22.rda")

#Open y_test and y_train datasets from cluster 16 
load("C:/Users/Aparajita/Documents/PScluster16.22.rda")
y_test=c16.y_test
y_train=c16.y_train

#Run regression which uses unfiltered sites + PCs
start_time <- Sys.time()
pc.c16 =lapply(1:ncol(y_test),penregT)
end_time=Sys.time()
t=end_time-start_time;t

save(glm.y_test=y_test,glm.y_train=y_train,glm.x_train=x_train,glm.x_test=x_test,
     reg.list= pc.c16,sim.time=t,file=paste("c16PC.rda"))

############################# Extract results #############################

############ Functions #############
#alphas
alphas = function(res.list){alphas=map_dbl(res.list, "al", .default = NA)}
#lambdas
lamdas= function(res.list){map_dbl(res.list, "lam", .default = NA)}

#glm resuts
glm.res = function(res.list){
  preds=lapply(res.list, "[[" , "y_pred" )
  y_pred = matrix(NA,nrow=nrow(y_test),ncol=ncol(y_test))
  j=1
  for (i in 1:ncol(y_test)){
    y_pred[,i]=preds[[j]]
    j=j+1}
  glm.res = summarize_res(y_pred,y_test)
  return(glm.res)
}

#coefficients
coeffs= function(res.list){lapply(res.list , "[[" , "coef" )}

############ All results ###############
#The files loaded are the files saved from the regression ran above

#Unfiltered 
load("C:/Users/Aparajita/Documents/C16Unf.rda")
unf.al = alphas(t.c16)
unf.lambda = lamdas(t.c16)
unf.c16.res = glm.res(t.c16)

#PC
load("C:/Users/Aparajita/Documents/c16PC.rda")
pc.al = alphas(pc.c16)
pc.lambda = lamdas(pc.c16)
pc.c16.res = glm.res(pc.c16)

#Top 20% correlated sites 
load("C:/Users/Aparajita/Documents/c16Cor10.rda")
cor20.al=alphas(cor20.c16)
cor20.lambda = lamdas(cor20.c16)
cor20.c16.res = glm.res(cor20.c16)

#Save all results
save(unf.c16.res,pc.c16.res,cor20.c16.res,file="C16features.rda")




