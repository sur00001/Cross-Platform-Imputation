##############################################################################
#This code runs a tuned penalized regression and a kernal for cluster 10
#You can specify the cluster
##############################################################################

library(glmnet)
library(data.table)

#Import the cluster data
testdatpath="~/Missingness_Bioinformatics/For Aparajita/For Aparajita/Select CpG sites/Clustering/grid search CV/Data"
load(file=file.path(testdatpath,"ch22_c10.RData")) #Cluster 10

#check train dataset structure
d <- c10
dim(d) #733 392 - 733 markers in cluster, 391 samples 
d <- d[,-392] #remove last column with cluster no
dim(d) #733 391
d <- t(d)

#Training-Test Split (20:80)
set.seed(245)
n <- nrow(d) #391
train_rows <- sample(seq(n), size = .8 * n)
train <- d[ train_rows, ]
test  <- d[-train_rows, ]
str(train) #1:312, 1:733
str(test) #1:79, 1:733

#load beta450 and beta850 from original for train dataset
#Import dataset
load("~/Missingness_Bioinformatics/For Aparajita/For Aparajita/Data/betachr22_rename.RData")

#Impute mean for missing data (450K)
d.450 <- t(beta_450) 
for(i in 1:ncol(d.450)){d.450[is.na(d.450[,i]), i] <- mean(d.450[,i], na.rm = TRUE)}
beta_450=t(d.450) #rows: covariates, cols: samples

#Impute mean for missing data (850K) 
d.850 <- t(beta_850) 
for(i in 1:ncol(d.850)){d.850[is.na(d.850[,i]), i] <- mean(d.850[,i], na.rm = TRUE)}
beta_850=t(d.850) #rows: covariates, cols: samples

#check original beta_450 and beta_850 structure 
str(beta_450) #1:8132, 1:391
str(beta_850) #1:10117, 1:391

#transpose train dataset
train <- t(train)
str(train) #1:733, 1:312

#separate beta450 and beta850
tr450 = train[which(rownames(train)%in%rownames(beta_450)),]
tr850 = train[which(rownames(train)%in%rownames(beta_850)),]

#check structure of tr450 and tr850
str(tr450) #1:339, 1:312
str(tr850) #1:394, 1:312

#construct x_train and y_train
x_train = t(tr450)
y_train = t(tr850)

#check structures
str(x_train) #1:312, 1:339 
str(y_train) #1:312, 1:394

#check test dataset structure
str(test) #1:79, 1:733

#transpose test dataset
test <- t(test)
str(test) #1:733, 1:79

#separate beta450 and beta850
tst450 = test[which(rownames(test)%in%rownames(beta_450)),]
tst850 = test[which(rownames(test)%in%rownames(beta_850)),]

#check structure of tst450 and tst850
str(tst450) #1:339, 1:79
str(tst850) #1:394, 1:79

#construct x_test, y_test
x_test = t(tst450)
y_test = t(tst850)

#check structure of x_test, y_test
str(x_test) #1:79, 1:339
str(y_test) #1:79, 1:394

#Source files to run kernal
libpath=datpath="~/Missingness_Bioinformatics/PR,KN_Code and Data"
source(file.path(libpath,"screening_lib_mac.R")) 
source(file.path(libpath,"kernel_pred_lib.R"))

##########################Tuned PR fit########################################
y_pred = matrix(NA,nrow=nrow(x_test),ncol=ncol(y_test))

#Store tuned parameters and coefficients for each regression
coeff.list=list()
T.alphas=c()
T.lambdas=c()

#Tuning alpha AND lambda
for(i in ncol(y_test)){ #whole chromosome
  #print(i)
  if(!any(is.na(y_train[,i]))){
    alpha=-1;cvm = 10000000;lambda = -1
    for(a in c()) {
      cv <- cv.glmnet(x_train, y_train[,i], alpha=a, type.measure = "mse", nfolds = 10)
      if(min(cv$cvm)<cvm) {
        cvm <- min(cv$cvm)
        alpha <- a
        lambda <- cv$lambda.min
        fit=glmnet(x_train,y_train[,i],alpha=a,lambda=cv$lambda.min)
        y_pred[,i]=predict(fit,newx=x_test)
        coeff.list[[i]]=coef(fit)[,1]
        T.alphas[i]=alpha
        T.lambdas[i]=lambda}}}}

#save all components to view later
save(glmy_pred=y_pred,glmy_test=y_test,
     glmy_train=y_train,glmx_train=x_train,glmx_test=x_test,glmcoef=coeff.list,
     glm.alpha=T.alphas,glm.lambda=T.lambdas,
     file=paste("KMeanC2.22.rda"))

###############################Kernal###################################
#Kernal 
rbfker <- rbfdot(sigma = 0.015)
system.time(ker_pred <- kernel_pred(x_train,y_train,x_test,ker=rbfker,alpha=0.99))
ker_res = summarize_res(ker_pred,y_test)
ker_res$tab  #0.017384059

#Penalized regression
load("C:/Users/Aparajita/Documents/KMeanC1.22.rda")
glm_res=summarize_res(y_pred,y_test)
res=glm_res$R2[1:140]
