# Using all covariates gives good estimates for drl. Only one column gives bad for all
library(ggplot2)
library(reshape2)
library(plyr)
library(ggpubr)
library(gridExtra)
library(Jmisc)
it = 8
n0 = 8000
K = 2
l = 1
source("C:/Users/manja/Dropbox/conditional_capture_recapture/codes/indep_cov_simulation.R")
source("C:/Users/manja/Dropbox/conditional_capture_recapture/codes/estimate_pseudoor.R")
realtrue = TRUE
simuldraw = 100

result1 = numeric(0)

set.seed(1)

alpha = 0.25
for(n0 in c(1000, 5000, 10000, 20000, 30000, 50000, 100000)){
  for(omega in c(1)){print(c(n0, omega))
    for(s in 1:simuldraw){
      
      datap = dat_p(n0, l)
      datap$theta0
      List_matrix = datap$List_matrix
      v = 1:l
      
      #est_val = estim_psix(List_matrix, K = K, funcname = funcname, v = v, realtrue = realtrue)
      est_val = estim_pseudoor(List_matrix, n = n0, K = K, v = v, alpha = alpha, omega = omega, plot = FALSE)
      
      result1 <- rbind(result1, cbind(est_val$vvec, est_val$plugin, est_val$drl, est_val$plugin.or, est_val$drl.or, est_val$gammainv.or, alpha, omega, n0)[est_val$set1 == 2,])
      #result2 <- rbind(result2, cbind(est_val$vvec, est_val$plugin.or, est_val$drl.or, est_val$gammainv.or, alpha, omega))
      
    }
  }
}
colnames(result1)[(length(v)+1):ncol(result1)] = c("plugin", "drl", "plugin.or", "drl.or", "gammainv.or", "alpha", "omega", "n0")
#colnames(result2)[(length(v)+1):ncol(result2)] = c("value", "plugin.or", "drl.or", "gammainv.or", "alpha", "omega")

result1 = as.data.frame(result1)
result2 = melt(result1, measure.vars = c("plugin", "drl", "drl.or"), variable.name = "method")

p1 = unlist(apply(as.matrix(result2[,1:l]), 1, pi1))
p2 = unlist(apply(as.matrix(result2[,1:l]), 1, pi2))
result2$gammainv0 = 1/(1 - (1-p1)*(1-p2))
result <- result2
colnames(result)[1] = "x1"

resultmse <- ddply(result, ~alpha + omega + method + n0, summarise,
                   mse = mean((value - plugin.or)^2),
                   mseg = mean((value - gammainv.or)^2),
                   mse0 = mean((value - gammainv0)^2))
resultmse <- melt(resultmse, id.vars = c("alpha", "omega", "method", "n0"))
theta0 = dat_p(n0, l)$theta0

#save(theta0, alpha, result1, resultmse, file = paste0("C:/Users/manja/Dropbox/conditional_capture_recapture/codes/data_K2_l1_psi0", round(theta0,2)*100, "_n_omega_mse_alpha0", round(alpha*100,0), ".Rdata"))
g1 <- ggplot(resultmse[resultmse$variable == "mse",], aes(x = n0, y = n0/l*value, color = method)) +
  geom_line() +
  geom_point() +
  ylab("n0*MSE") +
  labs(title = substitute(paste("MSE as a function of n for ", alpha, " = ", alp, ", ", psi, " = ", theta0, sep = ''), list(alp = alpha, theta0 = round(theta0, 2)))) +
  theme(legend.position = "bottom") +
  #facet_wrap(~variable) +
  scale_color_manual(values = c("plugin" = "red", "drl" = "goldenrod", "drl.or" = "grey"))
#  geom_line(aes(y = msedrl), color = "gray", linetype = "dashed")

pdf(paste0("C:/Users/manja/Dropbox/conditional_capture_recapture/codes/images/crc_simulated_barplots/lineplot_K2_l1_psi0", round(theta0,2)*100, "_n_omega_mse_alpha0", round(alpha*100,0), ".pdf"), width = 4.5, height = 4.5, onefile = FALSE) #height = 375 and 415
g1
dev.off()
