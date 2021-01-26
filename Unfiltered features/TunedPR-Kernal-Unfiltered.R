###################################################################################
#This code runs a tuned penalized regression and a kernal all 450K chromosome 22 
#sites as features (unfiltered sites)
#Note: this code was written for PCA but setting PCA=1 means no PCA

#Libraries
library(glmnet)
library(sys)

#Opening R file of methods for pen.reg and kernal 
libpath = "~/Missingness_Bioinformatics/PR,KN_Code and Data"
source(file.path(libpath,"screening_lib_mac.R")) 
source(file.path(libpath,"kernel_pred_lib.R"))

####################Load and setup data for chromosome######################
#Declare parameters for scenario
#1=no PCA, 2=PCA on genome, 3=PCA on chromosome
#pca.preds=TRUE mean that the PCA components are the only covariates
PCA=1
chromo=22
pca.preds=FALSE

#Open data
ch=paste("betachr",chromo,"_rename.RData",sep = "")
load(paste("~/Missingness_Bioinformatics/PR,KN_Code and Data/",ch,sep=""))

#Impute mean for missing data (450K) 
d.450 <- t(beta_450) 
for(i in 1:ncol(d.450)){d.450[is.na(d.450[,i]), i] <- mean(d.450[,i], na.rm = TRUE)}
d.450=t(d.450) #rows: covariates, cols: samples

##########################PCA##############################
#Note: "nd" or "npca" are labels I specified when this function ran histograms (you can ignore these)

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
  pca.df[is.na(pca.df[,i]), i] <- mean(pca.df[,i], na.rm = TRUE)
}

sum(is.na(pca.df))

############################Training and Testing datasets####################
#Training-Test Split (20-80 split)
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
#Tuning alpha AND lambda
penreg<-function(i){
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
    return(list(al=T.alphas,lam=T.lambdas,coef=coeff,y_pred=y_pred))}
}

######################### GLMNET with Unfiltered sites ###############################
#Just choose 200 random sites to predict (to reduce runtime)
#Note: I use the same random 200 sites for all prediction models
load("C:/Users/Aparajita/Documents/Missingness_Bioinformatics/200randomsites.rda")

start_time <- Sys.time()
unf.glm =lapply(r.sites,penreg)
end_time=Sys.time()
t=end_time-start_time;t

save(glm.y_test=y_test,glm.y_train=y_train,glm.x_train=x_train,glm.x_test=x_test,r.sites=r.sites,
     reg.list= unf.glm,sim.time=t,file=paste("UnfilteredT.22.rda"))

#Get GLM results
source(file.path(libpath,"GetGLMResults.R"))
unfPR.res = glm.res(unf.glm,r=TRUE,r.sites)

######################### Kernal #############################
#Note: Sigest() is a function to give a range of sigma parameters for kernal 

###### For entire chromosome #####
#Tune sigma 
sigest(beta_450)

#Tune sigma
#Trying sigma values within sigest range that maximize median R^2
s = r = r.new = 0 #Initialize parameters
for (i in seq(.0005,.139,.0005)){
  rbfker <- rbfdot(sigma = i)
  y_pred <- kernel_pred(x_train,y_train,x_test,ker=rbfker,alpha=0.99)
  ker_res = summarize_res(y_pred,y_test)
  r.new = median(ker_res$R2)
  if (r.new > r){
    r = r.new
    s=i
    }
}

s=.001
rbfker <- rbfdot(sigma = s)
system.time(y_pred <- kernel_pred(x_train,y_train,x_test,ker=rbfker,alpha=0.99))
Unf22.ker = summarize_res(y_pred,y_test); Unf22.ker$tab 

###### For 200 random sites #####
load("C:/Users/Aparajita/Documents/Missingness_Bioinformatics/200randomsites.rda")
r.y_train=y_train[,r.sites]
r.y_test=y_test[,r.sites]
s=.0014
rbfker <- rbfdot(sigma = s)
system.time(y_pred <- kernel_pred(x_train,r.y_train,x_test,ker=rbfker,alpha=0.99))
unfker.res = summarize_res(y_pred,r.y_test);unfker.res$tab 

