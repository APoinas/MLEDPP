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
w = Window(Simu[[1]])
vol = area(w)

#Setting up the progress bar
pb = txtProgressBar(max = n, style = 3)
progress = function(n) setTxtProgressBar(pb, n)
opts = list(progress = progress)

#Loading the values of the estimates of alpha
T2 = readRDS(paste("Version propre pour github/results/Result_Gauss_",toString(alpha0*100),"e-2_[0_",toString(S_length),"].dat",sep=""))
old_result = T2$MLE_corrected

#Computing the corresponding Fisher information
cl=makeCluster(nb_core)
clusterExport(cl, c('pairdist'))
registerDoSNOW(cl)
templist=foreach(a=1:n, .combine="cbind", .options.snow = opts) %dopar% {
  S = Simu[[a]]
  M = pairdist(S, periodic=TRUE)
  res = old_result[a]
  Fish = Fisher_Info(S$n/vol, res, M, 50)
  c(res, S$n/vol, Fish[1,1], Fish[1,2], Fish[2,2])}
stopCluster(cl)
T = unname(templist[1,])
Intens = unname(templist[2,])
F1 = unname(templist[3,])
F2 = unname(templist[4,])
F3 = unname(templist[5,])

all_data = data.frame(rho_est=Intens, alpha_est=T, Fisher_1_1=F1, Fisher_1_2=F2, Fisher_2_2=F3)
saveRDS(all_data ,paste("Version propre pour github/results/Fisher/Result_Gauss_",toString(alpha0*100),"e-2_[0_",toString(S_length),"].dat",sep=""))

#Number of correct IC for the estimate of alpha using the partial information matrix
range = 1.96/sqrt(F3)
sum(alpha0 <= T+range & alpha0 >= T-range)/n

boxplot(1/sqrt(F3))
abline(h=sd(T), col="red", lty=2)

#1, 0.01, 88.2%
#1, 0.03, 89.6%
#1, 0.05, 86%
#2, 0.01, 88.6%
#2, 0.03, 94.2%
#2, 0.05, 90.8%
#3, 0.01, 93.2%
#3, 0.03, 93.2%
#3, 0.05, 92.2%

#Number of correct IC for the estimate of alpha using the full information matrix
range = 1.96 * sqrt(1/(F3-((F2^2)/F1)))
sum(alpha0 <= T+range & alpha0 >= T-range)/n

boxplot(sqrt(1/(F3-((F2^2)/F1))))
abline(h=sd(T), col="red", lty=2)

############################################## Printing the results ##########################################################
