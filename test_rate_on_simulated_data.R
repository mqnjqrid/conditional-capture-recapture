# Using all covariates gives good estimates for drl. Only one column gives bad for all
library(ggplot2)
library(reshape2)
library(plyr)
library(ggpubr)
library(gridExtra)
library(Jmisc)
it = 6
n0 = 30000
K = 2
l = 1
source("C:/Users/manja/Dropbox/conditional_capture_recapture/codes/indep_cov_simulation.R")
source("C:/Users/manja/Dropbox/conditional_capture_recapture/codes/estimate_pseudoor.R")
realtrue = TRUE
simuldraw = 100

result1 = numeric(0)

set.seed(1)

for(alpha in seq(0.1, 0.5, by = 0.05)){
  for(omega in c(1)){print(c(alpha, omega))
    for(s in 1:simuldraw){

      datap = dat_p(n0, l)
      datap$theta0
      List_matrix = datap$List_matrix
      v = 1:l
      
      #est_val = estim_psix(List_matrix, K = K, funcname = funcname, v = v, realtrue = realtrue)
      est_val = estim_pseudoor(List_matrix, n = n0, K = K, v = v, alpha = alpha, omega = omega, plot = FALSE)
      
      result1 <- rbind(result1, cbind(est_val$vvec, est_val$plugin, est_val$drl, est_val$plugin.or, est_val$drl.or, est_val$gammainv.or, alpha, omega)[est_val$set1 == 2,])
      #result2 <- rbind(result2, cbind(est_val$vvec, est_val$plugin.or, est_val$drl.or, est_val$gammainv.or, alpha, omega))
      
    }
  }
}
colnames(result1)[(length(v)+1):ncol(result1)] = c("plugin", "drl", "plugin.or", "drl.or", "gammainv.or", "alpha", "omega")
#colnames(result2)[(length(v)+1):ncol(result2)] = c("value", "plugin.or", "drl.or", "gammainv.or", "alpha", "omega")

result1 = as.data.frame(result1)
result2 = melt(result1, measure.vars = c("plugin", "drl", "drl.or"), variable.name = "method")

p1 = unlist(apply(as.matrix(result2[,1:l]), 1, pi1))
p2 = unlist(apply(as.matrix(result2[,1:l]), 1, pi2))
result2$gammainv0 = 1/(1 - (1-p1)*(1-p2))
result <- result2
colnames(result)[1] = "x1"

resultmse <- ddply(result, ~alpha + omega + method, summarise,
                   mse = n0/l*mean((value - plugin.or)^2),
                   mseg = n0/l*mean((value - gammainv.or)^2),
                   mse0 = n0/l*mean((value - gammainv0)^2))
resultmse <- melt(resultmse, id.vars = c("alpha", "omega", "method"))
theta0 = dat_p(n0, l)$theta0

save(theta0, result1, result, resultmse, file = paste0("C:/Users/manja/Dropbox/conditional_capture_recapture/codes/data_K2_l1_psi0", round(theta0,2)*100, "_alpha_omega_mse_n", n0, ".Rdata"))

g1 <- ggplot(resultmse[resultmse$variable == "mse",], aes(x = alpha, y = value, color = method)) +
  geom_line() +
  geom_point() +
  ylab("n*MSE") +
  xlab(bquote(alpha)) +
  labs(title = substitute(paste("MSE as a function of error rate, ", psi, " = ", theta0, sep = ''), list(theta0 = round(theta0, 2)))) +
  theme(legend.position = "bottom") +
  #facet_wrap(~variable) +
  scale_color_manual(values = c("plugin" = "red", "drl" = "goldenrod", "drl.or" = "grey"))
#  geom_line(aes(y = msedrl), color = "gray", linetype = "dashed")

pdf(paste0("C:/Users/manja/Dropbox/conditional_capture_recapture/codes/images/crc_simulated_barplots/lineplot_K2_l1_psi0", round(theta0,2)*100, "_alpha_omega_mse_n", n0, ".pdf"), width = 4.5, height = 4.5, onefile = FALSE) #height = 375 and 415
g1
dev.off()


ggplot(result, aes(x = x1, y = value, color = method)) +
  geom_point() +
  facet_wrap(~method) +#, labeller = label_both) +
  scale_color_manual(values = c("plugin" = "red", "drl" = "goldenrod", "drl.or" = "grey")) +
  #geom_line(aes(x = v, y = plugin.or), color = "green")
  geom_smooth(aes(x = x1, y = plugin.or), method = "glm", color = "black") +
  geom_smooth(aes(x = x1, y = gammainv0), method = "glm", color = "green")
