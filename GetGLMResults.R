###############################################################################
#This code has the functions to find the glm results


################### Functions to find the glm results #################
#alphas
alphas = function(res.list){alphas=map_dbl(res.list, "al", .default = NA)}

#lambdas
lamdas= function(res.list){map_dbl(res.list, "lam", .default = NA)}

#glm resuts
glm.res = function(res.list,r=FALSE,r.sites){ #r=FALSE whole chromosome
  preds=lapply(res.list, "[[" , "y_pred" )
  y_pred = matrix(NA,nrow=nrow(y_test),ncol=ncol(y_test))
  if (r==TRUE){ #random sites are used 
    j=1
    for (i in r.sites){
      y_pred[,i]=preds[[j]]
      j=j+1}
    glm.res = summarize_res(y_pred[,r.sites],y_test[,r.sites])
  } else{
  j=1
  for (i in 1:ncol(y_test)){
    y_pred[,i]=preds[[j]]
    j=j+1}
  glm.res = summarize_res(y_pred,y_test)
  }

  return(glm.res)
}

#coefficients
coeffs= function(res.list){lapply(res.list , "[[" , "coef" )}