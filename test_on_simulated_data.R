library(ggplot2)
library(reshape2)
library(plyr)
library(ggpubr)
library(gridExtra)
library(Jmisc)
it = 9
n0 = 10000
l = 3
source("C:/Users/manja/Dropbox/conditional_capture_recapture/codes/indep_cov_simulation.R")
source("C:/Users/manja/Dropbox/conditional_capture_recapture/codes/estimate_real_data.R")
realtrue = TRUE
simuldraw = 50

result1 = numeric(0)
result2 = numeric(0)
funcname = c("logit"
             #, "sl"
             )
for(s in 1:simuldraw){
  print(s)
  
  datap = dat_p(n0, l)
  datap$theta0
  List_matrix = datap$List_matrix
  v = 1:l
  
  est_val = estim_psix(List_matrix, K = K, funcname = funcname, v = v, realtrue = realtrue)
  
  result1 <- rbind(result1, cbind(est_val$vvec, est_val$plugin, est_val$plugin.or, est_val$drl.or))
  result2 <- rbind(result2, cbind(est_val$vvec, est_val$drl, est_val$plugin.or, est_val$drl.or))
  
}
colnames(result1)[c((1 + length(v) +length(funcname)):ncol(result1))] = c("plugin.or", "drl.or")
colnames(result2)[c((1 + length(v) +length(funcname)):ncol(result2))] = c("plugin.or", "drl.or")

result1 = result1[order(result1[,1]),]
result1 = as.data.frame(result1)
result1 = melt(result1, measure.vars = funcname)
result1$method = "plugin"

result2 = result2[order(result2[,1]),]
result2 = as.data.frame(result2)
result2 = melt(result2, measure.vars = funcname)
result2$method = "drl"

result <- rbind(result1, result2)
ggplot(result, aes(x = x2, y = value, color = method)) +
  geom_smooth(se = FALSE) +
  scale_color_manual(values = c("plugin" = "red", "drl" = "goldenrod")) +
  geom_smooth(aes(y = plugin.or), se = FALSE, method = "glm", color = "black") +
  geom_smooth(aes(y = drl.or), method = "glm", color = "brown") +
  facet_wrap(~variable)

