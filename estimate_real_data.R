library(crctmle)
library(glmnet)
estim_psix = function(List_matrix, K, funcname = c("logit"), v, realtrue = FALSE){
  
  if(class(v) == "character"){
    v = which(colnames(List_matrix) == v) - K
  }
  
  #removing all rows with only 0's
  List_matrix = List_matrix[which(rowSums(List_matrix[,1:K]) > 0),]
  #N = number of observed or captured units
  l = ncol(List_matrix) - K
  N = nrow(List_matrix)
  
  colnames(List_matrix) = c(paste("L", 1:K, sep = ''), paste("x", 1:l, sep = ''))
  
  set1 = sample(1:N, ceiling(N)/2, replace = FALSE)
  List1 = as.data.frame(List_matrix[-set1,])
  List2 = as.data.frame(List_matrix[set1,])

  xmat = as.matrix(List_matrix[,-c(1:K)])
  xmat1 = as.matrix(List2[,-c(1:K)])
  #if(missing(i)){
    i = 1
  #}
  #if(missing(j)){
    j = 2
  #}
  
  plugin = matrix(nrow = N, ncol = length(funcname))
  colnames(plugin) = funcname
  drl = plugin
  for(func in funcname){
    qhat = get(paste0("qhat_", func))(List1, List2, K = 2, i, j)
    q12 = qhat$q12
    q1 = qhat$q1
    q2 = pmax(qhat$q2, q12/q1)
    
    ####################################
    yi = List2[,paste("L", i, sep = '')]
    yj = List2[,paste("L", j, sep = '')]
    ########### usual estimates
    gammainv_hat = q1*q2/q12
    
    if(length(v) == 1){
      plugin[, colnames(plugin)==func] = predict(smooth.spline(xmat1[,v], gammainv_hat), xmat[,v])$y
    }else{
      plugin[, colnames(plugin)==func] = predict(cv.glmnet(xmat1[,v], gammainv_hat), xmat[,v], type="response", s="lambda.min")
    }
    pseudo = plugin[set1, colnames(plugin)==func]*(yj/q2 + yi/q1 - yi*yj/q12)
    
    if(length(v) == 1){
      drl[,func] = predict(smooth.spline(xmat1[,v],pseudo), xmat[,v])$y
    }else{
      drl[,func] = predict(cv.glmnet(xmat1[,v], pseudo), xmat[,v], type="response", s="lambda.min")
    }
  }
  if(realtrue == TRUE){
    p1 = unlist(apply(xmat1, 1, pi1))
    p2 = unlist(apply(xmat1, 1, pi2))
    
    q1 = p1/(1 - (1-p1)*(1-p2))
    q2 = p2/(1 - (1-p1)*(1-p2))
    q12 = p2*p1/(1 - (1-p1)*(1-p2))
    
    ####### oracle estimates
    plugin.or = q1*q2/q12
    
    if(length(v) == 1){
      plugin.or = predict(smooth.spline(xmat1[,v], plugin.or), xmat[,v])$y
    }else{
      plugin.or = predict(cv.glmnet(xmat1[,v], plugin.or), xmat[,v], type="response", s="lambda.min")
    }
    
    pseudo.or = plugin.or[set1,]*(yj/q2 + yi/q1 - yi*yj/q12)
    if(length(v) == 1){
      drl.or = predict(smooth.spline(xmat1[,v],pseudo.or), xmat[,v])$y
    }else{
      drl.or = predict(cv.glmnet(xmat1[,v],pseudo.or), xmat[,v], type="response", s="lambda.min")
    }
  }else{
    plugin.or = NULL
    drl.or = NULL
  }
  return(list(vvec = xmat[,v], plugin = plugin, drl = drl, plugin.or = plugin.or,
              drl.or = drl.or))
}
