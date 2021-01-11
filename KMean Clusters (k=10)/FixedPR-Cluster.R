##################################################################################
#This file performs a fixed penalized regression (lambda is not tuned). The penalized
#regression depends on alpha (alpha=0,.5,1 is ridge, elastic net and lasso). 
#Note: any alpha between 0 and 1 is elastic net. There are different clusters
#in the data (see the Data folder) and this script reads in cluster 1 ("c1").
#The script will be the same for all clusters so you can import in whatever cluster
#you would like!
#################################################################################

library(glmnet)
library(data.table)
library(sys)

#Import training dataset
testdatpath = "~/Missingness_Bioinformatics/For Aparajita/For Aparajita/Select CpG sites/Clustering/grid search CV/Data"
load(file=file.path(testdatpath,"ch22_c1.RData"))

#check train dataset structure
d <- c1
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

#Load beta450 and beta850 from original for train dataset
#load(file=file.path(testdatpath,"betachr22_rename.RData"))
load("~/Missingness_Bioinformatics/For Aparajita/For Aparajita/Data/betachr22_rename.RData")

#Impute mean for missing data (450K) 
d.450 <- t(beta_450) 
for(i in 1:ncol(d.450)){d.450[is.na(d.450[,i]), i] <- mean(d.450[,i], na.rm = TRUE)}
beta_450=t(d.450)

#Impute mean for missing data (850K) 
d.850 <- t(beta_850) 
for(i in 1:ncol(d.850)){d.850[is.na(d.850[,i]), i] <- mean(d.850[,i], na.rm = TRUE)}
beta_850=t(d.850) 

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

#l is lambda and a is alpha 
#l=0 since I am only interested in lambda=0 but other values can be assessed
for (l in c(0)){ 
  for (a in c(0,.5,1)){ #do lasso and elastic net and ridge
    print(a)
    start.time=Sys.time()
    coeff.list=list()
    y_pred = matrix(NA,nrow=nrow(x_test),ncol=ncol(y_test))
    for(i in 1:ncol(y_test)){
      print(i)
      if(!any(is.na(y_train[,i]))){ 
        fit=glmnet(x_train,y_train[,i],alpha=a, lambda = l)
        y_pred[,i]=predict(fit,newx=x_test)
        coeff.list[[i]]=coef(fit)[,1]}}
    
    sum.res=summarize_res(y_pred,y_test)
    m=(sum(is.na(sum.res$R2)))/ncol(y_test)
    end_time=Sys.time();t=end_time-start_time
    
    #Save the coefficients used to predict each EPIC site, t=simulation time
    #m = # of missing R^2s
    if (a==0){ridge_fp=sum.res;ridge_pred=y_pred;ridge_coef=coeff.list;r.t=t;r.mis=m} else 
      if (a==.5){el_fp=sum.res;el_pred=y_pred;el_coef=coeff.list;el.t=t;el.mis=m} else 
      {l_fp=sum.res;l_pred=y_pred;l_coef=coeff.list;l.t=t;l.mis=m}}
  save(ridge_fp,ridge_pred,ridge_coef,r.t,el_fp,el_pred,el_coef,el.t,l_fp,l_pred,l_coef,
       l.t,r.mis,el.mis,l.mis,lambda=l,file=paste("KMeanFC1.22",l,".rda",sep=""))}
