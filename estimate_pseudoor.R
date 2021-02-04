pow = 1
iter = 30
library(glmnet)
estim_pseudoor = function(List_matrix, n, K, omega, alpha, v = v, plot = FALSE){
  
  if(class(v) == "character"){
    v = which(colnames(List_matrix) == v) - K
  }
  delta = 0
  #removing all rows with only 0's
  List_matrix = List_matrix[which(rowSums(List_matrix[,1:K]) > 0),]
  colnames(List_matrix) = c(paste("L", 1:K, sep = ''), paste("x", 1:(ncol(List_matrix) - K), sep = ''))
  #N = number of observed or captured units
  l = ncol(List_matrix) - K
  N = nrow(List_matrix)
  
  eta = (n/l)^alpha
  
  set1 = sample(1:2, N, prob = c(1/2, 1/2), replace = TRUE)
  xmat = as.matrix(List_matrix[,-c(1:K)])
  
  #calculating true nuisance functions
  p1 = unlist(apply(xmat, 1, pi1))
  p2 = unlist(apply(xmat, 1, pi2))
  
  q1_0 = p1/(1 - (1-p1)*(1-p2))
  q2_0 = p2/(1 - (1-p1)*(1-p2))
  q12_0 = p1*p2/(1 - (1-p1)*(1-p2))
  ####################################
  i=1;j=2
  
  yi = List_matrix[,paste("L", i, sep = '')]
  yj = List_matrix[,paste("L", j, sep = '')]
  
  ####### oracle estimates
  q1 = q1_0
  q2 = q2_0
  q12 = q12_0
  
  gammainv.or = q1*q2/q12
  if(length(v) == 1){
    #plugin.or = predict(smooth.spline(xmat[set1 == 1,v], gammainv.or[set1 == 1]), xmat[,v])$y
    xtrain = data.frame(v1 = xmat[,v], gammainv.or = gammainv.or)
    plugin.or = mgcv::predict.gam(mgcv::gam(gammainv.or ~ s(v1), data = xtrain[set1 == 1,]), newdata = xtrain)
  }else{
    plugin.or = predict(cv.glmnet(xmat[set1 == 1,v], gammainv.or[set1 == 1]), xmat[,v], type="response", s="lambda.1se")
  }
  
  pseudo.or = gammainv.or*(yj/q2 + yi/q1 - yi*yj/q12)
  if(length(v) == 1){
    #drl.or = predict(smooth.spline(xmat[set1 == 1,v], pseudo.or[set1 == 1]), xmat[,v])$y
    xtrain = data.frame(v1 = xmat[,v], pseudo.or = pseudo.or)
    drl.or = mgcv::predict.gam(mgcv::gam(pseudo.or ~ s(v1), data = xtrain[set1 == 1,]), newdata = xtrain)
  }else{
    drl.or = predict(cv.glmnet(xmat[set1 == 1,v], pseudo.or[set1 == 1]), xmat[,v], type="response", s="lambda.min")
  }
  # plot(plugin.or, drl.or)
  # plot(gammainv.or, drl.or)
  #  X = seq(min(xmat[,1]), max(xmat[,1]), length.out = 50)
  #  Y = seq(min(xmat[,2]), max(xmat[,2]), length.out = 50)
  #  ggg = outer(X, Y, Vectorize(function(x1,x2){1 - (1-pi1(c(x1,x2)))*(1-pi2(c(x1,x2)))}))
  #  plot_ly(x = ~X, y = ~Y, z = ~1/ggg, type = "surface")
  #  ggplot(as.data.frame(xmat), aes(x = x1, y = x2)) + stat_density2d_filled()
   
  ########### usual estimates
  epsiln = matrix(rnorm(3*length(q12_0), 1/eta, omega/eta), ncol = 3)

  q12 = expit(logit(q12_0) + epsiln[,3])
  q1 = pmax(pmin(expit(logit(q1_0) + epsiln[,1]), 1), q12)
  q2 = pmax(pmin(expit(logit(q2_0) + epsiln[,2]), 1 + q12 - q1), q12/q1)
  
  gammainv_hat = q1*q2/q12
  if(length(v) == 1){
    #plugin = predict(smooth.spline(xmat[set1 == 1,v], gammainv_hat[set1 == 1]), xmat[,v])$y
    xtrain = data.frame(v1 = xmat[,v], gammainv_hat = gammainv_hat)
    plugin = mgcv::predict.gam(mgcv::gam(gammainv_hat ~ s(v1), data = xtrain[set1 == 1,]), newdata = xtrain)
  }else{
    plugin = predict(cv.glmnet(xmat[set1 == 1,v], gammainv_hat[set1 == 1]), xmat[,v], type="response", s="lambda.min")
  }
  
  pseudo = gammainv_hat*(yj/q2 + yi/q1 - yi*yj/q12)
  if(length(v) == 1){
    #drl = predict(smooth.spline(xmat[set1 == 1,v], pseudo[set1 == 1]), xmat[,v])$y
    xtrain = data.frame(v1 = xmat[,v], pseudo = pseudo)
    drl = mgcv::predict.gam(mgcv::gam(pseudo ~ s(v1), data = xtrain[set1 == 1,]), newdata = xtrain)
  }else{
    drl = predict(cv.glmnet(xmat[set1 == 1,v], pseudo[set1 == 1]), xmat[,v], type="response", s="lambda.min")
  }

  if(plot == TRUE){
    ymin = min(drl, plugin, drl.or, 1/(1-(1-p1)*(1-p2)))
    ymax = max(drl, plugin, drl.or, 1/(1-(1-p1)*(1-p2)))
    plot(plugin.or, 1/(1-(1-p1)*(1-p2)), col = "yellow", pch = 19, cex = 0.4, ylim = c(ymin, ymax))
    points(plugin.or, drl, col = "red", pch = 19, cex = 0.4)
    points(plugin.or, plugin, col = "green", pch = 19, cex = 0.4)
    points(plugin.or, drl.or, pch = 19, cex = 0.4)
  }
  return(list(plugin = plugin, plugin.or = plugin.or, gammainv.or = gammainv.or,
              pseudo = pseudo, pseudo.or = pseudo.or,
              drl = drl, drl.or = drl.or, vvec = xmat[,v], set1 = set1))
}

ymin = min(drl, plugin, drl.or, 1/(1-(1-p1)*(1-p2)))
ymax = max(drl, plugin, drl.or, 1/(1-(1-p1)*(1-p2)))
par(mfrow = c(2,2), mar = c(2,2,2,0)) 
plot(plugin.or, 1/(1-(1-p1)*(1-p2)), col = "yellow", pch = 19, cex = 0.4, ylim = c(ymin, ymax))
points(plugin.or, drl, col = "red", pch = 19, cex = 0.4)
points(plugin.or, plugin, col = "green", pch = 19, cex = 0.4)
points(plugin.or, drl.or, pch = 19, cex = 0.4)

plot(plugin.or, q1, pch = 20, cex = 0.4, ylim = range(q1))
points(plugin.or, q1_0, pch = 19, col ="orange", cex = 0.6)
plot(plugin.or, q2, pch = 20, cex = 0.4, ylim = range(q2))
points(plugin.or, q2_0, pch = 19, col ="pink", cex = 0.6)
plot(plugin.or, q12, pch = 20, cex = 0.4, ylim = range(q12))
points(plugin.or, q12_0, pch = 19, col ="goldenrod4", cex = 0.6)
#dev.off()
