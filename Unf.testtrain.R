###################################################################################
#This code saves the training and testing sets for the unfiltered feature set

######################## Open directory and files ############################
libpath = "~/Missingness_Bioinformatics/PR,KN_Code and Data"

#Getting unique sample plates
pheno = fread(input=file.path(libpath,"pheno_bmi_QC.csv"))  
uni_plate = unique(pheno$sample_plate)

#opening R file of methods for pen.reg and kernal 
source(file.path(libpath,"screening_lib_mac.R")) 
source(file.path(libpath,"kernel_pred_lib.R"))

############################### Open data and test/train split #########################
#Test and training data ids
test_id = which(pheno$sample_plate==uni_plate[1]) #Tsai_Sample_212
train_id = which(pheno$sample_plate!=uni_plate[1]) #Rest of the Tsai_Samples

#Open data
load("~/Missingness_Bioinformatics/PR,KN_Code and Data/betachr22_rename.RData")

#Impute missing
d.450 <- t(beta_450) 
for(i in 1:ncol(d.450)){d.450[is.na(d.450[,i]), i] <- mean(d.450[,i], na.rm = TRUE)}
d.450=t(d.450); beta_450=d.450

d.850 <- t(beta_850) 
for(i in 1:ncol(d.850)){d.850[is.na(d.850[,i]), i] <- mean(d.850[,i], na.rm = TRUE)}
d.850=t(d.850);beta_850=d.850

Unf.x_train = t(beta_450[,train_id,drop=FALSE])
Unf.y_train = t(beta_850[,train_id,drop=FALSE])

Unf.x_test = t(beta_450[,test_id,drop=FALSE])
Unf.y_test = t(beta_850[,test_id,drop=FALSE])

dim(Unf.x_test);dim(Unf.x_train);dim(Unf.y_test);dim(Unf.y_train)
sum(is.na(Unf.y_test));sum(is.na(Unf.y_train));sum(is.na(Unf.x_train));sum(is.na(Unf.x_test))

save(Unf.x_train,Unf.x_test,Unf.y_train,Unf.y_test,file=paste("Unf.testtrain22.rda"))
