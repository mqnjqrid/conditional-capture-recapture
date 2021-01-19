library(ggplot2)
library(reshape2)
library(plyr)
library(ggpubr)
library(gridExtra)
library(Jmisc)
source("C:/Users/manja/Dropbox/conditional_capture_recapture/codes/indep_cov_simulation.R")
source("C:/Users/manja/Dropbox/conditional_capture_recapture/codes/estimate_pseudoor.R")
simuldraw = 100
n = 1000;l = 2
alpha_vec = c(0.25)
omega_vec = c(0.5, 1)
datorg = numeric(0)
mseorg = numeric(0)
theta0 =  dat_p(10000, l)$theta0
psi0 = 1/theta0
theta0
for(n in c(5000, 10000, 15000)){
  print(n)
  for(alp in 1:length(alpha_vec)) {
    alpha = alpha_vec[alp]
    for(sig in 1:length(omega_vec)){print(c(alp, sig))
      omega = omega_vec[sig]
      K = 2;l = 2
    
      for (s in 1:simuldraw){
       # print(s)
      
        datap = dat_p(n, l)
        List_matrix = datap$List_matrix

        est_val = estim_pseudoor(List_matrix, n = n, K = 2, omega = omega, alpha = alpha)
        
        datmat = cbind(s, est_val$plugin, est_val$plugin.or, est_val$pseudo, est_val$pseudo.or,
                       est_val$drl, est_val$drl.or, alpha, omega, n)
        colnames(datmat) = c("iter", "plugin", "plugin.or", "pseudo", "pseudo.or", "drl", "drl.or", "alpha", "omega", 'n')
        datorg = rbind(datorg, datmat)
      
        mse = c(colMeans((datmat[,2:7] - datmat[,3])^2), alpha, omega, n)
        names(mse) = c("plugin", "plugin.or", "pseudo", "pseudo.or", "drl", "drl.or", "alpha", "omega", "n")
        mseorg = rbind(mseorg, mse)
      
      }
    }
  }
}
datorg = data.frame(datorg)
mseorg = data.frame(mseorg)
mseorg[,1:6] = sqrt(mseorg[,1:6])
mseorgagg = aggregate(cbind(plugin , plugin.or, pseudo, pseudo.or, drl, drl.or) ~ alpha + omega + n, data = mseorg, median)
datorgagg = aggregate(cbind(plugin , plugin.or, pseudo, pseudo.or, drl, drl.or) ~ iter + alpha + omega + n, data = datorg, mean)
datorgagg = melt(datorgagg, id.vars = c("iter", "alpha", "omega", "n"))
datorgagg$value = abs(datorgagg$value - psi0)
datorgagg1 = aggregate(value ~ variable + alpha + omega + n, data = datorgagg, mean)

mseorgagg = melt(mseorgagg, id.vars = c("alpha", "omega", "n"))

#save(theta0, psi0, mseorgagg, file = paste("C:/Users/manja/Dropbox/conditional_capture_recapture/codes/data_theta0", round(10*theta0), "0.Rdata", sep = ''))
############## manually calculating coverage due to error on line 211-212

#load(file = paste("C:/Users/manja/Dropbox/conditional_capture_recapture/codes/data_theta060.Rdata", sep = ''))
#theta0 = 0.6;psio = 1/theta0

tsize = 12

v1 = ggplot(data = mseorgagg[mseorgagg$variable %in% c("plugin", "drl", "drl.or"),]) +
  geom_bar(mapping =  aes(x = n, y = value, fill = variable), stat = "identity", color = "black", position = "dodge") +
  # ylim(c(0, 1)) + #yupperborder)) +
  theme(text = element_text(size = tsize), legend.position = "bottom") +
  labs(title = substitute(paste("MSE of ", gamma, '(v), ', theta, ' = ', var), list(var = round(theta0, 1))), x = NULL, y = NULL) + 
  facet_grid(omega ~ alpha, labeller = label_both, scales = "free_y")
#scale_fill_manual("Estimation method", values=c("red", "#E69F00", "#56B4E9", "gray"))
v1

v2 = ggplot(data = mseorgagg[mseorgagg$variable %in% c("plugin", "drl", "drl.or"),]) +
  geom_line(mapping =  aes(x = n, y = value, color = variable), stat = "identity") +
  # ylim(c(0, 1)) + #yupperborder)) +
  theme(text = element_text(size = tsize), axis.title.x=element_blank(), legend.position = "bottom") +
  labs(title = substitute(paste("MSE of ", gamma, '(v), ', theta, ' = ', var), list(var = round(theta0, 1))), x = NULL, y = NULL) + 
  facet_grid(omega ~ alpha, labeller = label_both, scales = "free_y")
v2

ggplot(data = datorgagg1) +
  geom_bar(mapping =  aes(x = n, y = value, fill = variable), stat = "identity", color = "black", position = "dodge") +
  # ylim(c(0, 1)) + #yupperborder)) +
  theme(text = element_text(size = tsize), axis.text.x = element_blank(), axis.ticks.x=element_blank(), axis.title.x=element_blank(), legend.position = "bottom") +
  labs(title = substitute(paste("MSE of ", gamma, '(v), ', theta, ' = ', var), list(var = round(theta0, 1))), x = NULL, y = NULL) + 
  facet_grid(omega ~ alpha, labeller = label_both, scales = "free_y")

#pdf("C:/Users/manja/Dropbox/capture_recapture/codes/images/crc_simulated_barplots/barplot_K2_l1_theta030_alpha_omega12_bias_rmse_cvrg.pdf", width = 1100/60, height = 415/60, onefile = FALSE) #height = 375 and 415
#png(paste("C:/Users/manja/Dropbox/conditional_capture_recapture/codes/images/crc_simulated_barplots/barplot_K2_l1_theta0", round(10*theta0), "0_alpha_omega12_mse_n.png", sep = ''), width = 500, height = 300) #height = 375 and 415
v1
#dev.off()
#png(paste("C:/Users/manja/Dropbox/conditional_capture_recapture/codes/images/crc_simulated_barplots/lineplot_K2_l1_theta0", round(10*theta0), "0_alpha_omega12_mse_n.png", sep = ''), width = 400, height = 240) #height = 375 and 415
v2
#dev.off()

