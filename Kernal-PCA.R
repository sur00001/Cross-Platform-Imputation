#Conducting PCR on Chromosome 22 and using the top 50 PCs in Kernal method
#This version includes the top 50 PCAs WITH the other covariates


#Opening R file of methods for pen.reg and kernal 
libpath = "~/Missingness_Bioinformatics/PR,KN_Code and Data"
source(file.path(libpath,"screening_lib_mac.R")) 
source(file.path(libpath,"kernel_pred_lib.R"))

#Get the 200 random EPIC sites that will be predicted
load("C:/Users/Aparajita/Documents/Missingness_Bioinformatics/200randomsites.rda")

#Run kernal
kernal<-function(chromo,PCA,pca.preds=FALSE,sig=.0015){
  #1=no PCA, 2=PCA on genome, 3=PCA on chromosome
  #pca.preds=TRUE mean that the PCA components are the only covariates
  
  ####################Load and setup data for chromosome######################
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
    #Looking at variation explained by PCs
    pca.var=pca$sdev^2 / sum(pca$sdev^2) #top 50 components explain 43.8% of variation
    pcs.sumvar=sum(pca.var[1:50]) 
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
  
  #Impute missing data 
  for(i in 1:ncol(pca.df)){
    pca.df[is.na(pca.df[,i]), i] <- mean(pca.df[,i], na.rm = TRUE)
    #Linearly scale the PC components when kernal is being used 
    pca.df[,i]=(pca.df[,i]-mean(pca.df[,i]))/(max(pca.df[,i])-min(pca.df[,i]))}
  
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
  
  #############Kernal fit################
  y_pred = matrix(NA,nrow=nrow(x_test),ncol=ncol(y_test))
  
  rbfker <- rbfdot(sigma = sig)
  system.time(y_pred <- kernel_pred(x_train,y_train[,r.sites],x_test,ker=rbfker,alpha=0.99))
  ker_res = summarize_res(y_pred,y_test[,r.sites])
  ker_res$tab
  
  # hist(ker_res$R2,main=paste("Kernal R-squareds - PCA from Ch",chromo),xlab="R2")
  # hist(ker_res$mse,main=paste("Kernal R-squareds - PCA from Ch",chromo),xlab="MSE")
  # min(ker_res$R2);max(ker_res$R2);median(ker_res$R2)
  # save(ker.21_res=ker_res,y21_pred=y_pred,y21_test=y_test,
  #      y21_train=y_train,x21_train=x_train,x21_test=x_test,
  #      file=paste("KerPCA_","Chr",chromo,"_",npca,nd,scale,".RData",sep=""))
  return(ker_res)
}

####################### Kernal Result ###########################
pca_ker.res=kernal(chromo=22,PCA=3,pca.preds=FALSE)










