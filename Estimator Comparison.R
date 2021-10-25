library(spatstat)
library(stats)

#For parallelization
library(foreach)
library(doSNOW)
library(parallel)

source("MLE_DPP.R") #Load the DPP function, make sure it is in the working directory.

#######################################################################################################
############################### Examples of using the function MLEDPP #################################
#######################################################################################################

###################################
##### On a rectangular window #####
###################################

### Gauss ###
rho0 = 100
alpha0 = 0.03
S_length = 2
S = simulate(dppGauss(lambda=rho0, alpha=alpha0, d=2), 1, W=owin(c(0,S_length), c(0,S_length)))
plot(S)

MLEDPP(S, "Gauss", edgecorr=FALSE)
MLEDPP(S, "Gauss", edgecorr=TRUE)

### Bessel ###
rho0 = 100
alpha0 = 0.03
S_length = 2
S = simulate(dppBessel(lambda=rho0, alpha=alpha0, sigma=0, d=2), 1, W=owin(c(0,S_length), c(0,S_length)))
plot(S)

MLEDPP(S, "Bessel", edgecorr=TRUE)

### Cauchy ###
rho0 = 100
alpha0 = 0.02
S_length = 2
S = simulate(dppCauchy(lambda=rho0, alpha=alpha0, nu=0.5, d=2), 1, W=owin(c(0,S_length), c(0,S_length)))
plot(S)

MLEDPP(S, "Cauchy", edgecorr=TRUE)

### Whittle-Matern ###
rho0 = 100
alpha0 = 0.015
S_length = 2
S = simulate(dppMatern(lambda=rho0, alpha=alpha0, nu=2, d=2), 1, W=owin(c(0,S_length), c(0,S_length)))
plot(S)

MLEDPP(S, "WM", edgecorr=TRUE, sigma=2)

#################################
##### On an R-shaped window #####
#################################

rho0 = 100
alpha0 = 0.05
Simu = readRDS(paste("datasets/100_Gauss_", toString(alpha0*100), "e-2_letterR.dat", sep=""))
S = Simu[[1]]
plot(S)

MLEDPP(S, "Gauss", edgecorr=FALSE)
MLEDPP(S, "Gauss", edgecorr=TRUE)

###############################################################################################
############################### Comparison between estimators #################################
###############################################################################################

##### Gauss #####
rho0 = 100
alpha0 = 0.05
S_length = 1
nb_core = 7
n = 500
Simu = readRDS(paste("datasets/2000_Gauss_", toString(rho0), "_", toString(alpha0*100), "e-2_[0_", toString(S_length), "].dat", sep=""))

#MCE with the pcf
T_pcf = unlist(unname(mclapply(Simu[1:n], function(x) dppm(x~1, dppGauss, method="mincon", statistic="pcf", statargs=list(divisor='d'), rmin=0.01, q=1/2)$fitted$fixedpar$alpha, mc.cores=nb_core)))

#MCE with K
T_K = unlist(unname(mclapply(Simu[1:n], function(x) dppm(x~1, dppGauss, method="mincon", statistic="K", rmin=0.01, q=1/2)$fitted$fixedpar$alpha, mc.cores=nb_core)))

#Setting up the progress bar
pb = txtProgressBar(max = n, style = 3)
progress = function(n) setTxtProgressBar(pb, n)
opts = list(progress = progress)

#MLE (Without tore correction)
cl = makeCluster(nb_core)
clusterExport(cl, c('pairdist', 'is.ppp', 'area', 'dppGauss'))
registerDoSNOW(cl)
templist = foreach(a=1:n, .combine="cbind", .options.snow = opts) %dopar% {
  MLEDPP(Simu[[a]], "Gauss", edgecorr=FALSE)$fixedpar$alpha}
stopCluster(cl)
T_MLE = unname(templist[1,])

#MLE (With tore correction)
cl = makeCluster(nb_core)
clusterExport(cl, c('pairdist', 'is.ppp', 'area', 'dppGauss'))
registerDoSNOW(cl)
templist = foreach(a=1:n, .combine="cbind", .options.snow = opts) %dopar% {
  MLEDPP(Simu[[a]], "Gauss", edgecorr=TRUE)$fixedpar$alpha}
stopCluster(cl)
T_MLE_corr = unname(templist[1,])

#Boxplot of all statistics
boxplot(T_MLE_corr-alpha0, T_MLE-alpha0, T_pcf-alpha0, T_K-alpha0, names=c("MLE w/ edge corr.", "MLE w/o edge corr.", "MCE w/ the pcf", "MCE w/ K"))
abline(h=0, col="red", lty=2)
#Mean square errors
c(10^4*mean((T_MLE_corr-alpha0)^2), 10^4*mean((T_MLE-alpha0)^2), 10^4*mean((T_pcf-alpha0)^2), 10^4*mean((T_K-alpha0)^2))

#Saving data
all_data = data.frame(MLE_not_corrected = T_MLE, MLE_corrected = T_MLE_corr , MCE_with_g = T_pcf , MCE_with_K = T_K)
saveRDS(all_data ,paste("results/Result_Gauss_", toString(alpha0*100), "e-2_[0_", toString(S_length), "].dat", sep=""))

##### Bessel #####
rho0 = 100
alpha0 = 0.05
S_length = 2
nb_core = 7
n = 500
Simu = readRDS(paste("datasets/2000_Bessel_",toString(rho0),"_",toString(alpha0*100),"e-2_[0_",toString(S_length),"].dat",sep=""))

#MCE with the pcf
T_pcf = unlist(unname(mclapply(Simu[1:n], function(x) dppm(x~1, dppBessel(sigma=0), method="mincon", statistic="pcf", statargs=list(divisor='d'), rmin=0.01, q=1/2)$fitted$fixedpar$alpha, mc.cores=nb_core)))

#MCE with K
T_K = unlist(unname(mclapply(Simu[1:n], function(x) dppm(x~1, dppBessel(sigma=0), method="mincon", statistic="K", rmin=0.01, q=1/2)$fitted$fixedpar$alpha, mc.cores=nb_core)))

#Setting up the progress bar
pb = txtProgressBar(max = n, style = 3)
progress = function(n) setTxtProgressBar(pb, n)
opts = list(progress = progress)

#MLE (Without tore correction)
cl = makeCluster(nb_core)
clusterExport(cl, c('pairdist', 'is.ppp', 'area', 'dppBessel'))
registerDoSNOW(cl)
templist = foreach(a=1:n, .combine="cbind", .options.snow = opts) %dopar% {
  MLEDPP(Simu[[a]], "Bessel", edgecorr=FALSE)$fixedpar$alpha}
stopCluster(cl)
T_MLE = unname(templist[1,])

#MLE (With tore correction)
cl = makeCluster(nb_core)
clusterExport(cl, c('pairdist', 'is.ppp', 'area', 'dppBessel'))
registerDoSNOW(cl)
templist = foreach(a=1:n, .combine="cbind", .options.snow = opts) %dopar% {
  MLEDPP(Simu[[a]], "Bessel", edgecorr=TRUE)$fixedpar$alpha}
stopCluster(cl)
T_MLE_corr = unname(templist[1,])

#Boxplot of all statistics
boxplot(T_MLE_corr-alpha0, T_MLE-alpha0, T_pcf-alpha0, T_K-alpha0, names=c("MLE w/ edge corr.", "MLE w/o edge corr.", "MCE w/ the pcf", "MCE w/ K"))
abline(h=0, col="red", lty=2)
#Mean square errors
c(10^4*mean((T_MLE_corr-alpha0)^2), 10^4*mean((T_MLE-alpha0)^2), 10^4*mean((T_pcf-alpha0)^2), 10^4*mean((T_K-alpha0)^2))

#Saving data
all_data = data.frame(MLE_not_corrected = T_MLE, MLE_corrected = T_MLE_corr , MCE_with_g = T_pcf , MCE_with_K = T_K)
saveRDS(all_data ,paste("results/Result_Bessel_", toString(alpha0*100), "e-2_[0_", toString(S_length), "].dat", sep=""))

##### Cauchy #####
rho0 = 100
alpha0 = 0.02
S_length = 1
nb_core = 7
n = 500
Simu = readRDS(paste("datasets/500_Cauchy_",toString(rho0),"_",toString(alpha0*1000),"e-3_[0_",toString(S_length),"].dat",sep=""))

#MCE with the pcf
T_pcf = unlist(unname(mclapply(Simu[1:n], function(x) dppm(x~1, dppCauchy(nu=0.5), method="mincon", statistic="pcf", statargs=list(divisor='d'), rmin=0.01, q=1/2)$fitted$fixedpar$alpha, mc.cores=nb_core)))

#MCE with K
T_K = unlist(unname(mclapply(Simu[1:n], function(x) dppm(x~1, dppCauchy(nu=0.5), method="mincon", statistic="K", rmin=0.01, q=1/2)$fitted$fixedpar$alpha, mc.cores=nb_core)))

#Setting up the progress bar
pb = txtProgressBar(max = n, style = 3)
progress = function(n) setTxtProgressBar(pb, n)
opts = list(progress = progress)

#MLE (Without tore correction)
cl = makeCluster(nb_core)
clusterExport(cl, c('pairdist', 'is.ppp', 'area', 'dppCauchy'))
registerDoSNOW(cl)
templist = foreach(a=1:n, .combine="cbind", .options.snow = opts) %dopar% {
  MLEDPP(Simu[[a]], "Cauchy", edgecorr=FALSE)$fixedpar$alpha}
stopCluster(cl)
T_MLE = unname(templist[1,])

#MLE (With tore correction)
cl = makeCluster(nb_core)
clusterExport(cl, c('pairdist', 'is.ppp', 'area', 'dppCauchy'))
registerDoSNOW(cl)
templist = foreach(a=1:n, .combine="cbind", .options.snow = opts) %dopar% {
  MLEDPP(Simu[[a]], "Cauchy", edgecorr=TRUE)$fixedpar$alpha}
stopCluster(cl)
T_MLE_corr = unname(templist[1,])

#Boxplot of all statistics
boxplot(T_MLE_corr-alpha0, T_MLE-alpha0, T_pcf-alpha0, T_K-alpha0, names=c("MLE w/ edge corr.", "MLE w/o edge corr.", "MCE w/ the pcf", "MCE w/ K"))
abline(h=0, col="red", lty=2)
#Mean square errors
c(10^4*mean((T_MLE_corr-alpha0)^2), 10^4*mean((T_MLE-alpha0)^2), 10^4*mean((T_pcf-alpha0)^2), 10^4*mean((T_K-alpha0)^2))

#Saving data
all_data = data.frame(MLE_not_corrected = T_MLE, MLE_corrected = T_MLE_corr , MCE_with_g = T_pcf , MCE_with_K = T_K)
saveRDS(all_data ,paste("results/Result_Cauchy_", toString(alpha0*1000), "e-3_[0_", toString(S_length), "].dat", sep=""))

##### Whittle-Matérn #####
rho0 = 100
alpha0 = 0.015
S_length = 1
nb_core = 7
n = 500
Simu = readRDS(paste("datasets/500_Matern_sigma=2_",toString(rho0),"_",toString(alpha0*1000),"e-3_[0_",toString(S_length),"].dat",sep=""))

#MCE with the pcf
T_pcf = unlist(unname(mclapply(Simu[1:n], function(x) dppm(x~1, dppMatern(nu=2), method="mincon", statistic="pcf", statargs=list(divisor='d'), rmin=0.01, q=1/2)$fitted$fixedpar$alpha, mc.cores=nb_core)))

#MCE with K
T_K = unlist(unname(mclapply(Simu[1:n], function(x) dppm(x~1, dppMatern(nu=2), method="mincon", statistic="K", rmin=0.01, q=1/2)$fitted$fixedpar$alpha, mc.cores=nb_core)))

#Setting up the progress bar
pb = txtProgressBar(max = n, style = 3)
progress = function(n) setTxtProgressBar(pb, n)
opts = list(progress = progress)

#MLE (Without tore correction)
cl = makeCluster(nb_core)
clusterExport(cl, c('pairdist', 'is.ppp', 'area', 'dppMatern'))
registerDoSNOW(cl)
templist = foreach(a=1:n, .combine="cbind", .options.snow = opts) %dopar% {
  MLEDPP(Simu[[a]], "WM", edgecorr=FALSE, sigma=2)$fixedpar$alpha}
stopCluster(cl)
T_MLE = unname(templist[1,])

#MLE (With tore correction)
cl = makeCluster(nb_core)
clusterExport(cl, c('pairdist', 'is.ppp', 'area', 'dppMatern'))
registerDoSNOW(cl)
templist = foreach(a=1:n, .combine="cbind", .options.snow = opts) %dopar% {
  MLEDPP(Simu[[a]], "WM", edgecorr=TRUE, sigma=2)$fixedpar$alpha}
stopCluster(cl)
T_MLE_corr = unname(templist[1,])

#Boxplot of all statistics
boxplot(T_MLE_corr-alpha0, T_MLE-alpha0, T_pcf-alpha0, T_K-alpha0, names=c("MLE w/ edge corr.", "MLE w/o edge corr.", "MCE w/ the pcf", "MCE w/ K"))
abline(h=0, col="red", lty=2)
#Mean square errors
c(10^4*mean((T_MLE_corr-alpha0)^2), 10^4*mean((T_MLE-alpha0)^2), 10^4*mean((T_pcf-alpha0)^2), 10^4*mean((T_K-alpha0)^2))

#Saving data
all_data = data.frame(MLE_not_corrected = T_MLE, MLE_corrected = T_MLE_corr , MCE_with_g = T_pcf , MCE_with_K = T_K)
saveRDS(all_data ,paste("results/Result_Matern_", toString(alpha0*1000), "e-3_[0_", toString(S_length), "].dat", sep=""))

##### On an R-shaped window #####
rho0 = 100
alpha0 = 0.05
nb_core = 7
n = 500
Simu = readRDS(paste("datasets/500_Gauss_", toString(alpha0*100), "e-2_letterR.dat", sep=""))

#MCE with the pcf
T_pcf = unlist(unname(mclapply(Simu[1:n], function(x) dppm(x~1, dppGauss, method="mincon", statistic="pcf", rmin=0.01, q=1/2)$fitted$fixedpar$alpha, mc.cores=nb_core)))

#MCE with K
T_K = unlist(unname(mclapply(Simu[1:n], function(x) dppm(x~1, dppGauss, method="mincon", statistic="K", rmin=0.01, q=1/2)$fitted$fixedpar$alpha, mc.cores=nb_core)))

#Setting up the progress bar
pb = txtProgressBar(max = n, style = 3)
progress = function(n) setTxtProgressBar(pb, n)
opts = list(progress = progress)

#MLE (Without tore correction)
cl = makeCluster(nb_core)
clusterExport(cl, c('pairdist', 'is.ppp', 'area', 'dppGauss'))
registerDoSNOW(cl)
templist = foreach(a=1:n, .combine="cbind", .options.snow = opts) %dopar% {
  MLEDPP(Simu[[a]], "Gauss", edgecorr=FALSE)$fixedpar$alpha}
stopCluster(cl)
T_MLE = unname(templist[1,])

#MLE (With tore correction)
cl = makeCluster(nb_core)
clusterExport(cl, c('pairdist', 'is.ppp', 'area', 'dppGauss', 'bdist.points'))
registerDoSNOW(cl)
templist = foreach(a=1:n, .combine="cbind", .options.snow = opts) %dopar% {
  MLEDPP(Simu[[a]], "Gauss", edgecorr=TRUE)$fixedpar$alpha}
stopCluster(cl)
T_MLE_corr = unname(templist[1,])

#Boxplot of all statistics
boxplot(T_MLE_corr-alpha0, T_MLE-alpha0, T_pcf-alpha0, T_K-alpha0, names=c("MLE w/ edge corr.", "MLE w/o edge corr.", "MCE w/ the pcf", "MCE w/ K"))
abline(h=0, col="red", lty=2)
#Mean square errors
c(10^4*mean((T_MLE_corr-alpha0)^2), 10^4*mean((T_MLE-alpha0)^2), 10^4*mean((T_pcf-alpha0)^2), 10^4*mean((T_K-alpha0)^2))

#Saving data
all_data = data.frame(MLE_not_corrected = T_MLE, MLE_corrected = T_MLE_corr , MCE_with_g = T_pcf , MCE_with_K = T_K)
saveRDS(all_data ,paste("results/Result_Gauss_", toString(alpha0*100), "e-2_letterR.dat", sep=""))

##############################################################################################################################
############################################## Plotting all results ##########################################################
##############################################################################################################################

### Gauss ###
par(mfrow=c(3,3), mar=c(1, 2.5, 2, 1))
cmain = 2.5
caxis = 1.5
lwd = 1.5
for (wind in c(1,2,3)){
  for (alpha in c(1,3,5)){
    titletext = bquote(alpha*"*="*.(alpha/100)*", W=[0,"*.(wind)*"]²"~"")
    T = readRDS(paste("results/Result_Gauss_", alpha, "e-2_[0_", wind,"].dat", sep=""))
    boxplot(T-alpha/100,yaxt="n", cex.main=cmain, main=titletext, xaxt="n", lwd=lwd)
    cat('--------------------------------------------------------\n')
    cat(paste("Mean square Errror (x10^4) for alpha=", alpha/100, " and window=[0,", wind, "]^2\n", sep=""))
    print(10^4 * colMeans((T-alpha/100)^2))
    abline(h=0, col="red", lty=2, lwd=lwd)
    axis(2, cex.axis=caxis)}}

### Bessel ###
par(mfrow=c(3,3), mar=c(1, 2.5, 2, 1))
cmain = 2.5
caxis = 1.5
lwd = 1.5
for (wind in c(1,2,3)){
  for (alpha in c(1,3,5)){
    titletext = bquote(alpha*"*="*.(alpha/100)*", W=[0,"*.(wind)*"]²"~"")
    T = readRDS(paste("results/Result_Bessel_", alpha, "e-2_[0_", wind,"].dat", sep=""))
    boxplot(T-alpha/100,yaxt="n", cex.main=cmain, main=titletext, xaxt="n", lwd=lwd)
    cat('--------------------------------------------------------\n')
    cat(paste("Mean square Errror (x10^4) for alpha=", alpha/100, " and window=[0,", wind, "]^2\n", sep=""))
    print(10^4 * colMeans((T-alpha/100)^2))
    abline(h=0, col="red", lty=2, lwd=lwd)
    axis(2, cex.axis=caxis)}}

#### Cauchy ###
par(mfrow=c(3,3), mar=c(1, 2.5, 2, 1))
cmain = 2.5
caxis = 1.5
lwd = 1.5
for (wind in c(1,2,3)){
  for (alpha in c(5,20,35)){
    titletext = bquote(alpha*"*="*.(alpha/1000)*", W=[0,"*.(wind)*"]²"~"")
    T = readRDS(paste("results/Result_Cauchy_", alpha, "e-3_[0_", wind,"].dat", sep=""))
    boxplot(T-alpha/1000,yaxt="n", cex.main=cmain, main=titletext, xaxt="n", lwd=lwd)
    cat('--------------------------------------------------------\n')
    cat(paste("Mean square Errror (x10^4) for alpha=", alpha/1000, " and window=[0,", wind, "]^2\n", sep=""))
    print(10^4 * colMeans((T-alpha/1000)^2))
    abline(h=0, col="red", lty=2, lwd=lwd)
    axis(2, cex.axis=caxis)}}

### Whittle-Matérn ###
par(mfrow=c(1,2), mar=c(1, 2.5, 2, 1))
cmain = 2.5
caxis = 1.5
lwd = 1.5
alpha = 15
for (wind in c(1,2)){
  titletext = bquote(alpha*"*="*.(alpha/1000)*", W=[0,"*.(wind)*"]²"~"")
  T = readRDS(paste("results/Result_Matern_", alpha, "e-3_[0_", wind,"].dat", sep=""))
  boxplot(T-alpha/1000,yaxt="n", cex.main=cmain, main=titletext, xaxt="n", lwd=lwd)
  cat('--------------------------------------------------------\n')
  cat(paste("Mean square Errror (x10^4) for alpha=", alpha/1000, " and window=[0,", wind, "]^2\n", sep=""))
  print(10^4 * colMeans((T-alpha/1000)^2))
  abline(h=0, col="red", lty=2, lwd=lwd)
  axis(2, cex.axis=caxis)}

### Gauss (R-shaped window) ###
par(mfrow=c(1,3), mar=c(1, 2.5, 2, 1))
cmain = 2.5
caxis = 1.5
lwd = 1.5
for (alpha in c(1,3,5)){
  titletext = bquote(alpha*"*="*.(alpha/100)*""~"")
  T = readRDS(paste("results/Result_Gauss_", alpha, "e-2_letterR.dat", sep=""))
  boxplot(T-alpha/100,yaxt="n", cex.main=cmain, main=titletext, xaxt="n", lwd=lwd)
  cat('--------------------------------------------------------\n')
  cat(paste("Mean square Errror (x10^4) for alpha=", alpha/100, "\n", sep=""))
  print(10^4 * colMeans((T-alpha/100)^2))
  abline(h=0, col="red", lty=2, lwd=lwd)
  axis(2, cex.axis=caxis)}