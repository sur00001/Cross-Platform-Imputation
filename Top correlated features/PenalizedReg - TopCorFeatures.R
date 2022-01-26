# This script runs penalized regression with features that are the top 20% correlated sites for each EPIC site  

#Libraries
library(data.table)
library(glmnet)
library(sys)

#Open other methods
libpath = "~/Missingness_Bioinformatics/PR,KN_Code and Data"
source(file.path(libpath,"screening_lib_mac.R")) 
source(file.path(libpath,"kernel_pred_lib.R"))

#start_time <- Sys.time()
####################Load and setup data for chromosome######################
#Declare parameters for scenario
#1=no PCA, 2=PCA on genome, 3=PCA on chromosome
#pca.preds=TRUE mean that the PCA components are the only covariates
PCA=1
chromo=22 
pca.preds=FALSE

#Open data
#testdatpath = "/home/julianw_queue/sur00001/"
#load(file=file.path(testdatpath,"betachr22_rename.RData"))

ch=paste("betachr",chromo,"_rename.RData",sep = "")
load(paste("~/Missingness_Bioinformatics/PR,KN_Code and Data/",ch,sep=""))

#Impute mean for missing data (450K) to run PCA
d.450 <- t(beta_450) 
for(i in 1:ncol(d.450)){d.450[is.na(d.450[,i]), i] <- mean(d.450[,i], na.rm = TRUE)}
d.450=t(d.450) #rows: covariates, cols: samples

##########################PCA##############################
if (PCA==3){
  nd=paste("Ch") #name of data where PCs came from
  pca <- prcomp(d.450, scale=TRUE) 
  npc <- 50	# number of PCs to save
  pc.top <- t(data.frame(unclass(pca$rotation)[, 1:npc])) #50 PCs with values for each of the 391 samples
}

if (pca.preds==TRUE){preds=pc.top;npca=paste("PC")} #npca=name of covariates
if (pca.preds==FALSE){
  if(PCA==1){preds=d.450; npca=paste("HM450");nd=""}
  else{preds=rbind(pc.top,d.450);npca=paste("PC450")}
}
dim(preds) #cols: 391

##############################Add EPIC sites - Impute missing data##################################

pca.df=rbind(preds,beta_850)
pca.df <- t(pca.df) 
dim(pca.df) #1:391, 1:10184
for(i in 1:ncol(pca.df)){
  pca.df[is.na(pca.df[,i]), i] <- mean(pca.df[,i], na.rm = TRUE)}
sum(is.na(pca.df))
dim(pca.df) #391 rows

############################Training and Testing datasets####################

#Training-Test Split (now every 20 sites, not the whole chromosome)
#need 0.2*391 patients for test dataset, 0.8*391 for training dataset
set.seed(245)
n <- nrow(pca.df);n #391
train_rows <- sample(seq(n), size = .8 * n)
train <- pca.df[ train_rows, ]
test  <- pca.df[-train_rows, ]

#transpose train dataset
train <- t(train)

#separate 450K sites and beta850
tr.pc_450 = train[which(rownames(train)%in%rownames(preds)),]
tr850 = train[which(rownames(train)%in%rownames(beta_850)),]

#construct x_train and y_train
x_train = t(tr.pc_450)
y_train = t(tr850)

#transpose test dataset
test <- t(test)
dim(test) 

#separate 450K sites and beta850
tst.pc_450 = test[which(rownames(test)%in%rownames(preds)),]
tst850 = test[which(rownames(test)%in%rownames(beta_850)),]

#construct x_test, y_test
x_test = t(tst.pc_450)
y_test = t(tst850)

#Check dimensions
dim(x_test);dim(y_test);
dim(x_train);dim(y_train)


############# Tuned PR fit ################
#Correlation function
# cor.probe=function(i){
#   cor.l=unlist(lapply(1:length(x_train[1,]),
#   function(j){abs(cor(x_train[,j],y_train[,i]))}))
# sites=which(cor.l>=quantile(cor.l, 0.99))
# return(sites)}

cor.probe=function(i){
  cor.l=unlist(lapply(1:length(x_train[1,]),
  function(j){abs(cor(x_train[,j],y_train[,i]))}))
return(cor.l)}

#Tuning alpha AND lambda
penregC<-function(i){
  print(i)
  #Store tuned parameters and coefficients for each regression
  coeff=c(); T.alphas=c(); T.lambdas=c(); y_pred=c()
  if(!any(is.na(y_train[,i]))){
    alpha=-1;cvm = 10000000;lambda = -1
    top.cors=cor.probe(i) #top correlations
    sites=which(top.cors>=quantile(top.cors, 0.80)) #top correlated sites
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
          y_pred=y_pred,tcors=top.cors[sites],x.sites=sites))}}

penregC10<-function(i){
  print(i)
  #Store tuned parameters and coefficients for each regression
  coeff=c(); T.alphas=c(); T.lambdas=c(); y_pred=c()
  if(!any(is.na(y_train[,i]))){
    alpha=-1;cvm = 10000000;lambda = -1
    top.cors=cor.probe(i) #top correlations
    sites=which(top.cors>=(tail(sort(top.cors),10)[1]))#top correlated sites
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


############################# Serial Sim ######################################
#set.seed(1996)
#r.sites=sort(sample.int(10117, 200))
load("C:/Users/17034/Documents/Missingness_Bioinformatics/PR,KN_Code and Data/CorUnf1v2.22.rda")

#Top 20% of sites
Tquan=20
start_time <- Sys.time()
#cfull.sim =lapply(1:ncol(y_test),penregC)
c20.sim =lapply(r.sites,penregC)
end_time=Sys.time()
t=end_time-start_time;t

#l=lapply(1:10,penreg)
save(glm.y_test=y_test,glm.y_train=y_train,glm.x_train=x_train,glm.x_test=x_test,r.sites=r.sites,
     reg.list= c20.sim,sim.time=t,Tquantile=Tquan,file=paste("CorUnf20.22.rda"))

#Top 10 sites
start_time <- Sys.time()
c10.U =lapply(r.sites,penregC10)
end_time=Sys.time()
t=end_time-start_time;t

#l=lapply(1:10,penreg)
save(glm.y_test=y_test,glm.y_train=y_train,glm.x_train=x_train,glm.x_test=x_test,
     reg.list= c10.U,sim.time=t,file=paste("Cor10Unf.rda"))

# ########################## Mclapply Parallel Sim ###################################
library("parallel"); RNGkind("L'Ecuyer-CMRG")
ncores = 24
set.seed(1996)
r.sites=sort(sample.int(10117, 200)) #selecting 200 random sites

start_time <- Sys.time()
#tuned_Unfiltered = mclapply(r.sites,penreg, mc.cores = ncores)
save(glm.y_test=y_test,glm.y_train=y_train,glm.x_train=x_train,glm.x_test=x_test,
     reg.list= tuned_Unfiltered,file=paste("UnfilteredTtest.22.rda"))
# # 
end_time=Sys.time()
end_time-start_time

###################### parlapply Parallel Sim ##############################
ncores = 3
#set.seed(1996)
#r.sites=sort(sample.int(10117, 10)) #selecting random sites
r.sites=c(1:5)

start_time <- Sys.time()

cl <- makeCluster(getOption("cl.cores", ncores))
par.setup <- parLapply(cl, r.sites,
                       function(r.sites) {
                         require(glmnet)
                         require(data.table)})
clusterExport(cl, c('penreg', 'y_train','y_test','r.sites','x_test','x_train') )
p.sim=parLapply(cl, r.sites, penreg)
stopCluster(cl)

end_time=Sys.time()
end_time-start_time
