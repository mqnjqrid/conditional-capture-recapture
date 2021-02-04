library(ggplot2)
library(ggpubr)
library(robustbase)

#devtools::install("C:/Users/manja/Onedrive/Documents/crctmle")
library(crctmle)
source("C:/Users/manja/Dropbox/capture_recapture/codes/indep_cov_Tilling_simulation.R")
source("C:/Users/manja/Dropbox/capture_recapture/codes/indep_cov_bias_sqmse_functions.R")
source("C:/Users/manja/Dropbox/capture_recapture/codes/superlearner_function.R")
source("C:/Users/manja/Dropbox/capture_recapture/codes/psi_indep_tmle_crossfitting_corrected.R")
source("C:/Users/manja/Dropbox/capture_recapture/codes/Peru_codes/peru_1.R")
K = 2
nfolds = 5
eps = 0.005

funcname = c("logit"
             #                "Mlogit",
             #,"NP"#,
             #                ,"Gam",
             #  ,"sl"
)

list2 = pmax(x[,16], x[,17])
misage = rep(1, length(age))
misage[is.na(age)] = 0
age_nona = age
age_nona[is.na(age_nona)] = 0
iddpto = round(ui/10000)
idprov = round(ui/100)

datap0 = cbind(x[,c(10)], list2, x$AnoHechos, age_nona, perpe1, misage, x[,c(2, 6)], ui, iddpto, idprov, strata)
datap0 = merge(datap0, cordarea, by.x = "strata", by.y = "strata", all.x = TRUE)
datap0 = merge(datap0, depcoord, by.x = "iddpto", by.y = "IDDPTO", all.x = TRUE)

datap0 = datap0[,-which(names(datap0) %in% c("ui", "iddpto", "idprov", "strata"))]
datap0$Sexo[which(is.na(datap0$Sexo))] = 'N'
datap0$misage = factor(datap0$misage)
datap0$Sexo = factor(datap0$Sexo)
datap0$perpe1 = factor(datap0$perpe1)
datap0$Situacion = factor(datap0$Situacion)
colnames(datap0)[3] = "AnoHechos"
datap0 = data.frame(datap0)
datap0 = datap0[datap0$AnoHechos >= 1980,]
datap0 = na.omit(datap0)

delcol = which(names(datap0) %in% c("perpe1", "misage", "Sexo", "Situacion", "ubina", "depna", "provna"))
datap_modelmatrix = cbind(datap0[,-delcol], model.matrix(~misage + Sexo + Situacion + perpe1, datap0))
List_matrix = List_matrix[List_matrix$Sexo %in% c('F', 'M', 'I'), ]

List_matrix = datap_modelmatrix
est_val = estim_psix(List_matrix, K = K, funcname = funcname, v = v, realtrue = TRUE)

result1 <- cbind(est_val$vvec, est_val$plugin)
result2 <- cbind(est_val$vvec, est_val$drl)

colnames(result1)[1] = c("v")
colnames(result2)[1] = c("v")

result1 = result1[order(result1[,1]),]
result1 = as.data.frame(result1)
result1 = melt(result1, id.vars = "v")
result1$method = "plugin"

result2 = result2[order(result2[,1]),]
result2 = as.data.frame(result2)
result2 = melt(result2, id.vars = "v")
result2$method = "drl"

result <- rbind(result1, result2)
orgsize = data.frame(table(x$v))
colnames(orgsize) = c(v, "n0")
result = merge(result, orgsize, all.x = "TRUE")
ggplot(result, aes(x = v, y = value, color = variable)) +
  geom_line() +
  geom_line(aes(x = v, y = n0), color = "green") +
  xlab("AnoHechos")


#save(psimat, psiestim, psibarplot, file = "C:/Users/manja/Dropbox/capture_recapture/codes/Peru_codes/killing_barplot2_eps0005_latlong_de_strata.Rdata")

#load("C:/Users/manja/Dropbox/capture_recapture/codes/Peru_codes/killing_barplot2_eps0005_latlong_de_strata.Rdata")
