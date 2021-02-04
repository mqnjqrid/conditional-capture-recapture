expit = function(x) {
  exp(x)/(1 + exp(x))
}
logit = function(x) {
  log(x/(1 - x))
}
#ep_vec = c(-3.906, -3.414, -3.009, -2.664, -2.31, -1.974, -1.545, -0.9765, -0.48)
#for l = 1, capture probabilities are 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9
ep_vec = c(-4.2, -3.906, -3.414, -3.009, -2.664, -2.31, -1.974, -1.545, -0.9765, -0.48, 0)
#l = 3: 0.3, 0.4
ep = ep_vec[it]
pi1 = function(x) {
  expit( ep + sum(c(0.4)*x))
}
pi2 = function(x) {
  expit( ep + sum(c(0.3)*x))
}
dat_p = function(n, l){
  x = matrix(runif(n*l, 0, 1) + 2, nrow = n, ncol = l)
  y1 = unlist(apply(x, 1, function(xi) {sample(c(0, 1), 1, replace = TRUE, prob = c( 1 - pi1(xi), pi1(xi)))}))
  y2 = unlist(apply(x, 1, function(xi) {sample(c(0, 1), 1, replace = TRUE, prob = c( 1 - pi2(xi), pi2(xi)))}))
  xp = x + matrix(rnorm(n*l, 0, 0.5), ncol = l)
  List_matrix = cbind(y1, y2, x)
  List_matrix_xstar = cbind(y1, y2, xp)
  
  p1 = unlist(apply(x, 1, pi1))
  p2 = unlist(apply(x, 1, pi2))
  
  q1 = p1/(1 - (1-p1)*(1-p2))
  q2 = p2/(1 - (1-p1)*(1-p2))
  q12 = p1*p2/(1 - (1-p1)*(1-p2))
  return(list(List_matrix = List_matrix, List_matrix_xstar = List_matrix_xstar,
              theta0 = 1 -  mean(apply(x, 1, function(xx){return((1 - pi1(xx))*(1 - pi2(xx)))})),
              psi0 = 1/(1 -  mean(apply(x, 1, function(xx){return((1 - pi1(xx))*(1 - pi2(xx)))})))
              ))
}

print(dat_p(1000, l)$theta0)
# Qnphi = mean(sapply(1:1, function(i) { 
#   x = matrix(
#     rnorm(n*l, 0, 1),
#     nrow = n, ncol = l)
#   y1 = unlist(apply(x, 1, function(xi) {sample(c(0, 1), 1, replace = TRUE, prob = c( 1 - pi1(xi), pi1(xi)))}))
#   y2 = unlist(apply(x, 1, function(xi) {sample(c(0, 1), 1, replace = TRUE, prob = c( 1 - pi2(xi), pi2(xi)))}))
#   mean((q1*q2/q12 *(y1/q1 + y2/q2 - y1*y2/q12 - 1))[pmax(y1, y2) > 0])
# }))