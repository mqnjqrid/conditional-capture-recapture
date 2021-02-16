# Using all covariates gives good estimates for drl. Only one column gives bad for all
library(ggplot2)
library(reshape2)
library(plyr)
library(ggpubr)
library(gridExtra)
library(Jmisc)

n0 = 2000*3
K = 2
l = 1
source("C:/Users/manja/Dropbox/conditional_capture_recapture/codes/estimate_pseudoor.R")
realtrue = TRUE
simuldraw = 16

result1 = numeric(0)

set.seed(1)

alpha = 0.25
for(it in c(1:11)){
  source("C:/Users/manja/Dropbox/conditional_capture_recapture/codes/indep_cov_simulation.R")
  
  for(omega in c(1)){print(c(it, omega))
    for(s in 1:simuldraw){
      
      datap = dat_p(n0, l)
      theta0 = datap$theta0
      List_matrix = datap$List_matrix
      v = 1
      
      #est_val = estim_psix(List_matrix, K = K, funcname = funcname, v = v, realtrue = realtrue)
      est_val = estim_pseudoor(List_matrix, n = n0, K = K, v = v, alpha = alpha, omega = omega, plot = FALSE)
      
      result1 <- rbind(result1, cbind(est_val$vvec, est_val$plugin, est_val$drl, est_val$plugin.or, est_val$drl.or, est_val$gammainv.or, est_val$drlk, est_val$drlk.or, alpha, omega, n0, round(theta0, 2))[est_val$set1 == 2,])
      #result2 <- rbind(result2, cbind(est_val$vvec, est_val$plugin.or, est_val$drl.or, est_val$gammainv.or, alpha, omega))
      
    }
  }
}
colnames(result1)[(length(v)+1):ncol(result1)] = c("plugin", "drl", "plugin.or", "drl.or", "gammainv.or", "drlk", "drlk.or", "alpha", "omega", "n0", "theta0")
#colnames(result2)[(length(v)+1):ncol(result2)] = c("value", "plugin.or", "drl.or", "gammainv.or", "alpha", "omega")

result1 = as.data.frame(result1)
result2 = melt(result1, measure.vars = c("plugin", "drl", "drl.or", "drlk", "drlk.or"), variable.name = "method")

p1 = unlist(apply(as.matrix(result2[,1:l]), 1, pi1))
p2 = unlist(apply(as.matrix(result2[,1:l]), 1, pi2))
result2$gammainv0 = 1/(1 - (1-p1)*(1-p2))
result <- result2
colnames(result)[1] = "x1"

resultmse <- ddply(result, ~alpha + omega + method + n0 + theta0, summarise,
                   mse = mean((value - plugin.or)^2),
                   mseg = mean((value - gammainv.or)^2),
                   mse0 = mean((value - gammainv0)^2))
resultmse <- melt(resultmse, id.vars = c("alpha", "omega", "method", "n0", "theta0"))

#save(result1, result2, resultmse, n0, l, file = paste0("C:/Users/manja/Dropbox/conditional_capture_recapture/codes/images/crc_simulated_barplots/data_K2_l1_n", n0, "_omega_mse_alpha0", round(alpha*100,0), ".Rdata"))

#### OLD location of Rdata files
#load("C:/Users/manja/Dropbox/conditional_capture_recapture/codes/data_K2_l1_n20000_omega_mse_alpha025.Rdata")
#### NEW location of Rdata files
#load("C:/Users/manja/OneDrive/Documents/codes_Rdata_from_conditional_crc_folder_codes/data_K2_l1_n20000_omega_mse_alpha025.Rdata")
g1 <- ggplot(resultmse[resultmse$variable == "mse",], aes(x = theta0, y = n0/l*value, color = method)) +
  geom_line() +
  geom_point() +
  ylab("n0*MSE") +
  xlab(bquote(psi~" (expected number of observations)")) +
  labs(title = substitute(paste("MSE as a function of ", psi, " for ", alpha, " = ", alp, " and n = ", n0, sep = ''), list(alp = alpha, n0 = n0))) +
  scale_x_continuous(breaks = seq(0,1, by = 0.1), labels = paste0(seq(0, 1, by = 0.1), ' (', n0*seq(0, 1, by = 0.1), ')')) +
  #scale_color_manual(values = c("plugin" = "red", "drl" = "goldenrod", "drl.or" = "grey")) +
  theme(legend.position = "bottom", text = element_text(size = 12))
  #facet_wrap(~variable) +
#  geom_line(aes(y = msedrl), color = "gray", linetype = "dashed")


g2 <- ggplot(result[abs(result$theta0 - 0.84) < 0.05,], aes(x = x1, y = value, color = method)) +
  geom_point(size = 0.3) +
  facet_wrap(~theta0) +
  labs(title = substitute(paste("MSE as a function of n for ", alpha, " = ", alp, sep = ''), list(alp = alpha))) +
  theme(legend.position = "bottom", text = element_text(size = 12)) +
  scale_color_manual(values = c("plugin" = "red", "drl" = "goldenrod", "drl.or" = "grey")) +
  #facet_wrap(~variable) +
  geom_smooth(aes(x = x1, y = plugin.or), method = "glm", color = "black")

#pdf(paste0("C:/Users/manja/Dropbox/conditional_capture_recapture/codes/images/crc_simulated_barplots/lineplot_K2_l1_n", n0, "_omega_mse_alpha0", round(alpha*100,0), ".pdf"), width = 9, height = 4.5, onefile = FALSE) #height = 375 and 415
g1 + scale_y_log10()
#dev.off()

