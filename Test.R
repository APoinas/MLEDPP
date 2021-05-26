library(foreach)
library(doSNOW)
library(parallel)

###############Gauss Check##################
rho0 = 100
alpha0 = 0.03
S_length = 1
nb_core = 7
Simu = readRDS(paste("datasets/2000_Gauss_",toString(rho0),"_",toString(alpha0*100),"e-2_[0_",toString(S_length),"].dat",sep=""))
n = 50

#Setting up the progress bar
pb = txtProgressBar(max = n, style = 3)
progress = function(n) setTxtProgressBar(pb, n)
opts = list(progress = progress)

#Loading the values of the estimates of alpha
T2 = readRDS(paste("results/Result_Gauss_",toString(alpha0*100),"e-2_[0_",toString(S_length),"].dat",sep=""))
old_result = T2$MLE_corrected

#Computing the corresponding Fisher information
cl=makeCluster(nb_core)
clusterExport(cl, c('pairdist', 'area'))
registerDoSNOW(cl)
templist=foreach(a=1:n, .combine="cbind", .options.snow = opts) %dopar% {
  Fish = Fisher_Info(Simu[[a]], "Gauss", old_result[a], edgecorr=TRUE)
  c(Fish[1,1], Fish[1,2], Fish[2,2])}
stopCluster(cl)
F1 = unname(templist[1,])
F2 = unname(templist[2,])
F3 = unname(templist[3,])

Check_List = readRDS(paste("results/Fisher/Result_Gauss_",toString(alpha0*100),"e-2_[0_",toString(S_length),"].dat",sep=""))
print(c(mean((F1-Check_List$Fisher_1_1[1:n])^2),mean((F2-Check_List$Fisher_1_2[1:n])^2),mean((F3-Check_List$Fisher_2_2[1:n])^2)))

###############Cauchy Check##################
rho0 = 100
alpha0 = 0.02
S_length = 1
nb_core = 7
Simu = readRDS(paste("datasets/500_Cauchy_",toString(rho0),"_",toString(alpha0*1000),"e-3_[0_",toString(S_length),"].dat",sep=""))
n = 50

#Setting up the progress bar
pb = txtProgressBar(max = n, style = 3)
progress = function(n) setTxtProgressBar(pb, n)
opts = list(progress = progress)

#Loading the values of the estimates of alpha
T2 = readRDS(paste("results/Result_Cauchy_",toString(alpha0*1000),"e-3_[0_",toString(S_length),"].dat",sep=""))
old_result = T2$MLE_corrected

#Computing the corresponding Fisher information
cl=makeCluster(nb_core)
clusterExport(cl, c('pairdist', 'area'))
registerDoSNOW(cl)
templist=foreach(a=1:n, .combine="cbind", .options.snow = opts) %dopar% {
  Fish = Fisher_Info(Simu[[a]], "Cauchy", old_result[a], edgecorr=TRUE)
  c(Fish[1,1], Fish[1,2], Fish[2,2])}
stopCluster(cl)
F1 = unname(templist[1,])
F2 = unname(templist[2,])
F3 = unname(templist[3,])

Check_List = readRDS(paste("results/Fisher/Result_Cauchy_",toString(alpha0*1000),"e-3_[0_",toString(S_length),"].dat",sep=""))
print(c(mean((F1-Check_List$Fisher_1_1[1:n])^2),mean((F2-Check_List$Fisher_1_2[1:n])^2),mean((F3-Check_List$Fisher_2_2[1:n])^2)))

###############Bessel Check##################
rho0 = 100
alpha0 = 0.03
S_length = 1
nb_core = 7
Simu = readRDS(paste("datasets/2000_Bessel_",toString(rho0),"_",toString(alpha0*100),"e-2_[0_",toString(S_length),"].dat",sep=""))
n = 50

#Setting up the progress bar
pb = txtProgressBar(max = n, style = 3)
progress = function(n) setTxtProgressBar(pb, n)
opts = list(progress = progress)

#Loading the values of the estimates of alpha
T2 = readRDS(paste("results/Result_Bessel_",toString(alpha0*100),"e-2_[0_",toString(S_length),"].dat",sep=""))
old_result = T2$MLE_corrected

#Computing the corresponding Fisher information
cl=makeCluster(nb_core)
clusterExport(cl, c('pairdist', 'area'))
registerDoSNOW(cl)
templist=foreach(a=1:n, .combine="cbind", .options.snow = opts) %dopar% {
  Fish = Fisher_Info(Simu[[a]], "Bessel", old_result[a], edgecorr=TRUE)
  c(Fish[1,1], Fish[1,2], Fish[2,2])}
stopCluster(cl)
F1 = unname(templist[1,])
F2 = unname(templist[2,])
F3 = unname(templist[3,])

Check_List = readRDS(paste("results/Fisher/Result_Bessel_",toString(alpha0*100),"e-2_[0_",toString(S_length),"].dat",sep=""))
print(c(mean((F1-Check_List$Fisher_1_1[1:n])^2),mean((F2-Check_List$Fisher_1_2[1:n])^2),mean((F3-Check_List$Fisher_2_2[1:n])^2)))
