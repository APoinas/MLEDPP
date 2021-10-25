library(spatstat)
library(stats)

#For parallelization
library(foreach)
library(doSNOW)
library(parallel)

source("MLE_DPP.R") #Load the DPP function, make sure it is in the working directory.

############################### Examples of using the function Fisher_Info #################################
#Gauss
rho0 = 100
alpha0 = 0.03
S_length = 2
S = simulate(dppGauss(lambda=rho0, alpha=alpha0, d=2), 1, W=owin(c(0,S_length), c(0,S_length)))
plot(S)

alpha_est = MLEDPP(S, "Gauss", edgecorr=TRUE)$fixedpar$alpha
I_fisher = Fisher_Info(S, "Gauss", alpha_est, edgecorr=TRUE)
Range_IC = 1.96 * sqrt(1/(I_fisher[2,2]-((I_fisher[1,2]^2)/I_fisher[1,1])))
print(paste("95% confidence Interval for the estimation of alpha:[", round(alpha_est - Range_IC, 3), " , ", round(alpha_est + Range_IC, 3), "]"))

#Bessel
rho0 = 100
alpha0 = 0.03
S_length = 2
S = simulate(dppBessel(lambda=rho0, alpha=alpha0, sigma=0, d=2), 1, W=owin(c(0,S_length), c(0,S_length)))
plot(S)

alpha_est = MLEDPP(S, "Bessel", edgecorr=TRUE)$fixedpar$alpha
I_fisher = Fisher_Info(S, "Bessel", alpha_est, edgecorr=TRUE)
Range_IC = 1.96 * sqrt(1/(I_fisher[2,2]-((I_fisher[1,2]^2)/I_fisher[1,1])))
print(paste("95% confidence Interval for the estimation of alpha:[", round(alpha_est - Range_IC, 3), " , ", round(alpha_est + Range_IC, 3), "]"))

#Cauchy
rho0 = 100
alpha0 = 0.02
S_length = 2
S = simulate(dppCauchy(lambda=rho0, alpha=alpha0, nu=0.5, d=2), 1, W=owin(c(0,S_length), c(0,S_length)))
plot(S)

alpha_est = MLEDPP(S, "Cauchy", edgecorr=TRUE)$fixedpar$alpha
I_fisher = Fisher_Info(S, "Cauchy", alpha_est, edgecorr=TRUE)
Range_IC = 1.96 * sqrt(1/(I_fisher[2,2]-((I_fisher[1,2]^2)/I_fisher[1,1])))
print(paste("95% confidence Interval for the estimation of alpha:[", round(alpha_est - Range_IC, 3), " , ", round(alpha_est + Range_IC, 3), "]"))

############################### Testing the correctness of the confidence interval #################################
#####Gauss#####
rho0 = 100
alpha0 = 0.03
S_length = 1
nb_core = 7
Simu = readRDS(paste("datasets/2000_Gauss_",toString(rho0),"_",toString(alpha0*100),"e-2_[0_",toString(S_length),"].dat",sep=""))
n = 500

#Setting up the progress bar
pb = txtProgressBar(max = n, style = 3)
progress = function(n) setTxtProgressBar(pb, n)
opts = list(progress = progress)

#Loading the values of the estimates of alpha
T2 = readRDS(paste("results/Result_Gauss_",toString(alpha0*100),"e-2_[0_",toString(S_length),"].dat",sep=""))
old_result = T2$MLE_corrected[1:n]

#Computing the corresponding Fisher information
cl=makeCluster(nb_core)
clusterExport(cl, c('pairdist','area'))
registerDoSNOW(cl)
templist=foreach(a=1:n, .combine="cbind", .options.snow = opts) %dopar% {
  S = Simu[[a]]
  alpha_est = old_result[a]
  Fish = Fisher_Info(S, "Gauss", alpha_est, edgecorr=TRUE)
  c(Fish[1,1], Fish[1,2], Fish[2,2])}
stopCluster(cl)
F1 = unname(templist[1,])
F2 = unname(templist[2,])
F3 = unname(templist[3,])

#Saving the obtained data
all_data = data.frame(rho_est=Intens, alpha_est=T, Fisher_1_1=F1, Fisher_1_2=F2, Fisher_2_2=F3)
saveRDS(all_data ,paste("results/Fisher/Result_Gauss_",toString(alpha0*100),"e-2_[0_",toString(S_length),"].dat",sep=""))

#Number of correct IC for the estimate of alpha using the full information matrix
print(paste(100*sum(F3<=0)/n, "% of confidence interval computations failed"))
alpha_est_list = old_result[F3>0]
F1 = F1[F3>0]
F2 = F2[F3>0]
F3 = F3[F3>0]
range = 1.96 * sqrt(1/(F3-((F2^2)/F1)))
print(paste(100*sum(alpha0 <= alpha_est_list+range & alpha0 >= alpha_est_list-range)/n, "% of the confidence intervals contains the true value of alpha"))

#####Bessel#####
rho0 = 100
alpha0 = 0.03
S_length = 1
nb_core = 7
Simu = readRDS(paste("datasets/2000_Bessel_",toString(rho0),"_",toString(alpha0*100),"e-2_[0_",toString(S_length),"].dat",sep=""))
n = 500

#Setting up the progress bar
pb = txtProgressBar(max = n, style = 3)
progress = function(n) setTxtProgressBar(pb, n)
opts = list(progress = progress)

#Loading the values of the estimates of alpha
T2 = readRDS(paste("results/Result_Bessel_",toString(alpha0*100),"e-2_[0_",toString(S_length),"].dat",sep=""))
old_result = T2$MLE_corrected[1:n]

#Computing the corresponding Fisher information
cl=makeCluster(nb_core)
clusterExport(cl, c('pairdist','area'))
registerDoSNOW(cl)
templist=foreach(a=1:n, .combine="cbind", .options.snow = opts) %dopar% {
  S = Simu[[a]]
  alpha_est = old_result[a]
  Fish = Fisher_Info(S, "Bessel", alpha_est, edgecorr=TRUE)
  c(Fish[1,1], Fish[1,2], Fish[2,2])}
stopCluster(cl)
F1 = unname(templist[1,])
F2 = unname(templist[2,])
F3 = unname(templist[3,])

#Saving the obtained data
all_data = data.frame(rho_est=Intens, alpha_est=T, Fisher_1_1=F1, Fisher_1_2=F2, Fisher_2_2=F3)
saveRDS(all_data ,paste("results/Fisher/Result_Bessel_",toString(alpha0*100),"e-2_[0_",toString(S_length),"].dat",sep=""))

#Number of correct IC for the estimate of alpha using the full information matrix
print(paste(100*sum(F3<=0)/n, "% of confidence interval computations failed"))
alpha_est_list = old_result[F3>0]
F1 = F1[F3>0]
F2 = F2[F3>0]
F3 = F3[F3>0]
range = 1.96 * sqrt(1/(F3-((F2^2)/F1)))
print(paste(100*sum(alpha0 <= alpha_est_list+range & alpha0 >= alpha_est_list-range)/n, "% of the confidence intervals contains the true value of alpha"))

##################
##### Cauchy #####
##################
rho0 = 100
alpha0 = 0.02
S_length = 1
nb_core = 7
Simu = readRDS(paste("datasets/500_Cauchy_",toString(rho0),"_",toString(alpha0*1000),"e-3_[0_",toString(S_length),"].dat",sep=""))
n = 500

#Setting up the progress bar
pb = txtProgressBar(max = n, style = 3)
progress = function(n) setTxtProgressBar(pb, n)
opts = list(progress = progress)

#Loading the values of the estimates of alpha
T2 = readRDS(paste("results/Result_Cauchy_",toString(alpha0*1000),"e-3_[0_",toString(S_length),"].dat",sep=""))
old_result = T2$MLE_corrected[1:n]

#Computing the corresponding Fisher information
cl=makeCluster(nb_core)
clusterExport(cl, c('pairdist','area'))
registerDoSNOW(cl)
templist=foreach(a=1:n, .combine="cbind", .options.snow = opts) %dopar% {
  S = Simu[[a]]
  alpha_est = old_result[a]
  Fish = Fisher_Info(S, "Cauchy", alpha_est, edgecorr=TRUE)
  c(Fish[1,1], Fish[1,2], Fish[2,2])}
stopCluster(cl)
F1 = unname(templist[1,])
F2 = unname(templist[2,])
F3 = unname(templist[3,])

all_data = data.frame(rho_est=Intens, alpha_est=T, Fisher_1_1=F1, Fisher_1_2=F2, Fisher_2_2=F3)
saveRDS(all_data ,paste("results/Fisher/Result_Cauchy_",toString(alpha0*1000),"e-3_[0_",toString(S_length),"].dat",sep=""))

#Number of correct IC for the estimate of alpha using the full information matrix
print(paste(100*sum(F3<=0)/n, "% of confidence interval computations failed"))
alpha_est_list = old_result[F3>0]
F1 = F1[F3>0]
F2 = F2[F3>0]
F3 = F3[F3>0]
range = 1.96 * sqrt(1/(F3-((F2^2)/F1)))
print(paste(100*sum(alpha0 <= alpha_est_list+range & alpha0 >= alpha_est_list-range)/n, "% of the confidence intervals contains the true value of alpha"))



############################################## Printing the results ##########################################################
