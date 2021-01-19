pow = 1
iter = 30

estim_pseudoor = function(List_matrix, n, K, omega, alpha){
  
  delta = 0
  eta = n^alpha
  #removing all rows with only 0's
  List_matrix = List_matrix[which(rowSums(List_matrix[,1:K]) > 0),]
  colnames(List_matrix) = c(paste("L", 1:K, sep = ''), paste("x", 1:(ncol(List_matrix) - K), sep = ''))
  #N = number of observed or captured units
  l = ncol(List_matrix) - K
  N = nrow(List_matrix)
  
  set1 = sample(1:N, ceiling(N)/2, replace = FALSE)
  List1 = as.data.frame(List_matrix[-set1,])
  List2 = as.data.frame(List_matrix[set1,])
  xmat = as.matrix(List2[,-c(1:K)])
  xmat1= as.matrix(List1[,-c(1:K)])
  p1 = unlist(apply(xmat, 1, pi1))
  p2 = unlist(apply(xmat, 1, pi2))
  
  q10_0 = p1*(1 - p2)/(1 - (1-p1)*(1-p2))
  q02_0 = p2*(1 - p1)/(1 - (1-p1)*(1-p2))
  q12_0 = p2*p1/(1 - (1-p1)*(1-p2))
  
  ####################################
  i=1;j=2
  
  fi = 1
  yi = List2[,paste("L", i, sep = '')]
  yj = List2[,paste("L", j, sep = '')]
  
  epsiln = matrix(rnorm(3*nrow(xmat), 1/eta, omega/eta), ncol = 3)
#  epsiln = rnorm(3, 0, omega)
  q12 = expit(logit(q12_0) + epsiln[,3])
  q1 = pmin(q12 + expit(logit(q10_0) + epsiln[,1]), 1)
  q2 = pmin(q12 + expit(logit(q02_0) + epsiln[,2]), 1 + q12 - q1)
  
  gammainv_hat = q1*q2/q12
  plugin = predict(smooth.spline(xmat[,ncol(xmat)], gammainv_hat),xmat1[,ncol(xmat)])$y
  
  phihat = #pmin(pmax(
    gammainv_hat*(yj/q2 + yi/q1 - yi*yj/q12)#, -100), 100)
  
  pseudo = phihat
  
  drl = predict(smooth.spline(xmat[,ncol(xmat)],pseudo),xmat1[,ncol(xmat)])$y
  #plot(drl, gammainv_hat, col = yi + 2*yj, pch = 19)
  #abline(0, 1)
  #plot(pseudo, gammainv_hat, col = yi + 2*yj, pch = 19)
  #abline(0, 1)
  #legend("bottomright", legend = c("10", "02", "12"), col = 1:3, pch = 19)
  
  q1 = q12_0 + q10_0
  q2 = q12_0 + q02_0
  q12 = q12_0
  gammainv_hat.or = q1*q2/q12
  plugin.or = gammainv_hat.or
  
  phihat.or = gammainv_hat.or*(yj/q2 + yi/q1 - yi*yj/q12) - plugin.or
  
  pseudo.or = plugin.or + phihat.or
  
  drl.or = predict(smooth.spline(xmat[,ncol(xmat)],pseudo.or),xmat1[,ncol(xmat)])$y
  
  return(list(plugin = plugin, plugin.or = plugin.or,
              pseudo = pseudo, pseudo.or = pseudo.or,
              drl = drl, drl.or = drl.or))
}
