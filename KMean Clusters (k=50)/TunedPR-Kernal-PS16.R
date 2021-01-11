#################################################################################
#This code runs the tuned penalized and kernal regression for cluster 16
#of the 50 K-Means clusters 
#Note: We preselected 8 of the 50 KMeans clusters. These 8 clusters had the highest
#pairwise correlations between the EPIC and 450K sites 

#Libraries
library(data.table)
library(glmnet)
library(purrr)

#Set file paths
testdatpath = "~/Missingness_Bioinformatics/Pre-selected Cluster/"
libpath="~/Missingness_Bioinformatics/PR,KN_Code and Data"

#Open methods
source(file.path(libpath,"screening_lib_mac.R"))
source(file.path(libpath,"kernel_pred_lib.R"))
source(file.path(libpath,"GetGLMResults.R"))

#Divide testing and training
pheno = fread(input=file.path(libpath,"pheno_bmi_QC.csv"))
uni_plate = unique(pheno$sample_plate)
test_id = which(pheno$sample_plate==uni_plate[1])
train_id = which(pheno$sample_plate!=uni_plate[1])

#Get cluster 16 data
#Need to load beta_450 and beta_850 data (Ch22, cluster16)
load(file=file.path(testdatpath,"ch22_cluster16_450.RData"))
load(file=file.path(testdatpath,"ch22_cluster16_850.RData"))

#create training set
x_train = t(cluster16_450[,train_id,drop=FALSE])
y_train = t(cluster16_850[,train_id,drop=FALSE])

#check properties of training set
ncol(x_train) #62
nrow(x_train) #326
ncol(y_train) #65
nrow(y_train) #326

#confirm no intersect
#a <- intersect(colnames(x_train), colnames(y_train))
#a #NA, no intersect

#create test set for cross-validation
x_test = t(cluster16_450[,test_id,drop=FALSE])
y_test = t(cluster16_850[,test_id,drop=FALSE])

#check properties of test set
ncol(x_test) #
nrow(x_test) #65
ncol(y_test) #
nrow(y_test) #65

################################## Correlations ######################################

plot(density(cor(x_train,y_train,method="pearson")))
summary(cor(x_train,y_train,method="pearson"))

##################### GLMNET ############################
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

#Tuned penalized regression
start_time=Sys.time()
PSc16.sim=lapply(1:ncol(y_test),penreg)
end_time <- Sys.time()
t=end_time-start_time;t

#Save file to get results
save(glm.y_test=y_test,glm.y_train=y_train,glm.x_train=x_train,glm.x_test=x_test,
     reg.list= PSc16.sim,sim.time=t,file=paste("PScluster16.22.rda"))

#Summary of results (R^2 and MSE)
PR16.res = glm.res(PSc16.sim)

################## Kernal ########################
#Note: I manually found the sigma that maximized the R^2

rbfker <- rbfdot(sigma = 0.4)
system.time(ker_pred <- kernel_pred(x_train,y_train,x_test,ker=rbfker,alpha=0.99))
ker16_res = summarize_res(ker_pred,y_test)
