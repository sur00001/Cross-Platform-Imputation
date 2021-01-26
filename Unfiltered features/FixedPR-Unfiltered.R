###################################################################################
#This code is for fixed penalized regression (fixed lambda vals) 
#with all 450K sites in chromosome 22 as features (unfiltered features)
#Note: this function was written for PCA but PCA=1 means no PCA is run


#Libraries
library(data.table)
library(glmnet)
library(sys)
start_time <- Sys.time()

####################Load and setup data for chromosome######################
#Declare parameters for scenario
#1=no PCA, 2=PCA on genome, 3=PCA on chromosome
#pca.preds=TRUE mean that the PCA components are the only covariates
PCA=1
chromo=22
pca.preds=FALSE

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

############# Fixed PR fit ################
libpath=datpath="~/Missingness_Bioinformatics/PR,KN_Code and Data"
source(file.path(libpath,"screening_lib_mac.R")) 
source(file.path(libpath,"kernel_pred_lib.R"))

# a=0; l=.5
# penregf<-function(i){
#   #Store tuned parameters and coefficients for each regression
#   coeff=c();y_pred=c()
#   fit=glmnet(x_train,y_train[,i],alpha=a,lambda=l)
#   y_pred=predict(fit,newx=x_test); rownames(y_pred)=NULL
#   y_pred=y_pred[,1] #Want the y_pred vector only 
#   coeff=coef(fit)[,1]
#   return(list(coef=coeff,y_pred=y_pred))}

for (l in c(.025,.075)){ #Trying multiple lambdas
  for (a in c(0,.5,1)){ #do lasso and elastic net and ridge
    start_time=Sys.time()
    coeff.list=list()
    y_pred = matrix(NA,nrow=nrow(x_test),ncol=ncol(y_test))
    for(i in r.sites){
      print(i)
      if(!any(is.na(y_train[,i]))){ 
        fit=glmnet(x_train,y_train[,i],alpha=a, lambda = l)
        y_pred[,i]=predict(fit,newx=x_test)
        coeff.list[[i]]=coef(fit)[,1]}}
    
    sum.res=summarize_res(y_pred,y_test)
    m=(sum(is.na(sum.res$R2)))/ncol(y_test) #how many missing R^2s
    end_time=Sys.time();t=end_time-start_time
    
    if (a==0){ridge_fp=sum.res;ridge_pred=y_pred;ridge_coef=coeff.list;r.t=t;r.mis=m} else 
      if (a==.5){el_fp=sum.res;el_pred=y_pred;el_coef=coeff.list;el.t=t;el.mis=m} else 
      {l_fp=sum.res;l_pred=y_pred;l_coef=coeff.list;l.t=t;l.mis=m}}
  save(ridge_fp,ridge_pred,ridge_coef,r.t,el_fp,el_pred,el_coef,el.t,l_fp,l_pred,l_coef,
       l.t,r.mis,el.mis,l.mis,lambda=l,file=paste("FUnf.",l,".rda",sep=""))}

# 
# ########################## Mclapply Parallel Sim (for MACs) ###################################
# library("parallel"); RNGkind("L'Ecuyer-CMRG")
# ncores = 100
# set.seed(1996)
# r.sites=sort(sample.int(10117, 1000)) #selecting 1000 random sites
# 
# # fixed_Unfiltered = mclapply(r.sites,penregf, mc.cores = ncores)
# save(glm.y_test=y_test,glm.y_train=y_train,glm.x_train=x_train,glm.x_test=x_test,
# reg.list= l,time=time,file=paste("UnfilteredF.22.rda")) 
# 
# end_time=Sys.time() 
# time=end_time-start_time






