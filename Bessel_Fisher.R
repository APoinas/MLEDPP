library(spatstat)
library(tictoc)
library(stats)
library(foreach)
library(doSNOW)
library(parallel)

#####Functions#####

#L0
L0=function(rho,alpha,M){
  N = (besselJ(2*M/alpha,1)/(M/alpha))*rho/(1-rho*pi*alpha^2)
  diag(N) = rho/(1-rho*pi*alpha^2)
  return(N)
}

#Derivative of L0 with respect to alpha
DalpL0=function(rho,alpha,M){
  N = (besselJ(2*M/alpha, 1)/(M/alpha))* 2*pi*alpha*rho^2 /(1-rho*pi*alpha^2)^2 + besselJ(2*M/alpha, 2) * 2*rho /(alpha * (1-rho*pi*alpha^2))
  diag(N) = 2*pi*alpha*rho^2/(1-rho*pi*alpha^2)^2
  return(N)
}

#2nd derivative of L0 with respect to alpha
D2alpL0=function(rho,alpha,M){
  N = besselJ(2*M/alpha, 1)* (2*pi*alpha*rho^2 * (1+ 3*pi*rho*alpha^2)/ (M * (1-rho*pi*alpha^2)^3) - 4 * rho * M /(alpha^3 * (1-pi*rho*alpha^2))) + besselJ(2*M/alpha, 2) * 2*rho*(1+3*pi*rho*alpha^2) /(alpha^2 * (1-rho*pi*alpha^2)^2)
  diag(N) = 2*pi*rho^2 * (1+ 3*pi*rho*alpha^2)/(1-rho*pi*alpha^2)^3
  return(N)
}

#Derivative of L0 with respect to rho
DrhoL0=function(rho,alpha,M){
  N = (besselJ(2*M/alpha,1)/(M/alpha))/(1-rho*pi*alpha^2)^2
  diag(N) = 1/(1-rho*pi*alpha^2)^2
  return(N)
}

#2nd derivative of L0 with respect to rho
D2rhoL0=function(rho,alpha,M){
  N = (besselJ(2*M/alpha,1)/(M/alpha)) * 2*pi*alpha^2/(1-rho*pi*alpha^2)^3
  diag(N) = 2*pi*alpha^2 /(1-rho*pi*alpha^2)^3
  return(N)
}

#Derivative of L0 with respect to alpha and rho
DalpDrhoL0=function(rho,alpha,M){
  N = (besselJ(2*M/alpha, 1)/(M/alpha))* 4*pi*alpha*rho /(1-rho*pi*alpha^2)^3 + besselJ(2*M/alpha, 2) * 2 /(alpha * (1-rho*pi*alpha^2)^2)
  diag(N) = 4*pi*alpha*rho /(1-rho*pi*alpha^2)^3
  return(N)
}

fun = function(h){
  return ( (DalpInteg(100,0.03+h) - DalpInteg(100,0.03))/h )
}
fun(0.00000000001)
D2alpInteg(100,0.03)

fun2 = function(h){
  return ( (DalpInteg(100+h,0.03) - DalpInteg(100,0.03))/h )
}
fun2(0.000000001)
DalpDrhoInteg(100,0.03)

#Integral part of the log-likelihood
Integ = function(rho,alpha){
  return (log(1-rho*pi*alpha^2)/(pi*alpha^2))
}

#Derivative of the integral part of the log-likelihood with respect to alpha
DalpInteg = function(rho,alpha){
  return ((rho*pi*alpha^2/(1-rho*pi*alpha^2) + log(1-rho*pi*alpha^2))* -2/(pi*alpha^3))
}

#2nd derivative of the integral part of the log-likelihood with respect to alpha
D2alpInteg = function(rho,alpha){
  return (2*rho*(3-5*pi*rho*alpha^2)/(alpha^2 * (1-rho*pi*alpha^2)^2) + log(1-rho*pi*alpha^2)* 6/(pi*alpha^4))
}

#Derivative of the integral part of the log-likelihood with respect to rho
DrhoInteg = function(rho,alpha){
  return (-1/(1-rho*pi*alpha^2))
}

#2nd derivative of the integral part of the log-likelihood with respect to rho
D2rhoInteg = function(rho,alpha){
  return (-pi*alpha^2/(1-rho*pi*alpha^2)^2)
}

#Derivative of the integral part of the log-likelihood with respect to rho and alpha
DalpDrhoInteg = function(rho,alpha){
  return (-2*rho*pi*alpha/(1-rho*pi*alpha^2)^2)
}

#Fisher information matrix
Fisher_Info = function(rho,alpha){
  if (rho*pi*alpha^2>1) {return(NA)}
  N = L0(rho,alpha,M)
  DalpN = DalpL0(rho,alpha,M)
  DrhoN = DrhoL0(rho,alpha,M)
  D2alpN = D2alpL0(rho,alpha,M)
  D2rhoN = D2rhoL0(rho,alpha,M)
  DalpDrhoN = DalpDrhoL0(rho,alpha,M)
  invN = solve(N)
  tempmat_alpha = DalpN %*% invN
  tempmat_rho = DrhoN %*% invN
  a = vol * D2rhoInteg(rho,alpha) + sum(D2rhoN * invN) - sum(tempmat_rho * t(tempmat_rho))
  b = vol * DalpDrhoInteg(rho,alpha) + sum(DalpDrhoN * invN) - sum(tempmat_alpha * t(tempmat_rho))
  c = vol * D2alpInteg(rho,alpha) + sum(D2alpN * invN) - sum(tempmat_alpha * t(tempmat_alpha))
  print(c(N[1,1], DalpN[1,1], DrhoN[1,1], D2alpN[1,1], D2rhoN[1,1], DalpDrhoN[1,1]))
  return(array(c(-a,-b,-b,-c), dim=c(2,2)))
}

##############################################################################################
rho0 = 100
alpha0 = 0.05
S_length = 3
nb_core = 7
Simu = readRDS(paste("datasets/2000_Bessel_",toString(rho0),"_",toString(alpha0*100),"e-2_[0_",toString(S_length),"].dat",sep=""))
n = 500
w = Window(Simu[[1]])
vol = area(w)

#Setting up the progress bar
pb = txtProgressBar(max = n, style = 3)
progress = function(n) setTxtProgressBar(pb, n)
opts = list(progress = progress)

#Loading the values of the estimates of alpha
#T2 = readRDS(paste("Version propre pour github/results/Result_Bessel_",toString(alpha0*100),"e-2_[0_",toString(S_length),"].dat",sep=""))
T2 = readRDS(paste("../Git_MLEDPP/MLEDPP/results/Result_Bessel_",toString(alpha0*100),"e-2_[0_",toString(S_length),"].dat",sep=""))
old_result = T2$MLE_corrected

#Computing the corresponding Fisher information
cl=makeCluster(nb_core)
clusterExport(cl, c('pairdist'))
registerDoSNOW(cl)
templist=foreach(a=1:n, .combine="cbind", .options.snow = opts) %dopar% {
  S = Simu[[a]]
  M = pairdist(S, periodic=TRUE)
  res = old_result[a]
  Fish = Fisher_Info(S$n/vol, res)
  c(res, S$n/vol, Fish[1,1], Fish[1,2], Fish[2,2])}
stopCluster(cl)
T = unname(templist[1,])
Intens = unname(templist[2,])
F1 = unname(templist[3,])
F2 = unname(templist[4,])
F3 = unname(templist[5,])

all_data = data.frame(rho_est=Intens, alpha_est=T, Fisher_1_1=F1, Fisher_1_2=F2, Fisher_2_2=F3)
saveRDS(all_data ,paste("Version propre pour github/results/Fisher/Result_Bessel_",toString(alpha0*100),"e-2_[0_",toString(S_length),"].dat",sep=""))

sum(F3<=0)/n
sum(pi*Intens*T^2>0.99999)/n

#1, 0.01, 0% / 0%
#1, 0.03, 0.2% / 0.2%
#1, 0.05, 16.2% / 16.2%
#2, 0.01, 0% / 0%
#2, 0.03, 0% / 0%
#2, 0.05, 7.6% / 7.6%
#3, 0.01, 0% / 0%
#3, 0.03, 0% / 0%
#3, 0.05, 2.2% / 2.2%

T = T[F3>0]
Intens = Intens[F3>0]
F1 = F1[F3>0]
F2 = F2[F3>0]
F3 = F3[F3>0]

#Number of correct IC for the estimate of alpha using the partial information matrix
range = 1.96/sqrt(F3)
sum(alpha0 <= T+range & alpha0 >= T-range)/n

boxplot(1/sqrt(F3))
abline(h=sd(T), col="red", lty=2)

#1, 0.01, 66%
#1, 0.03, 76.6%
#1, 0.05, 56%
#2, 0.01, 77.2%
#2, 0.03, 78.8%
#2, 0.05, 82.2%
#3, 0.01, 81%
#3, 0.03, 74.4%
#3, 0.05, 12%

#Number of correct IC for the estimate of alpha using the full information matrix
range = 1.96 * sqrt(1/(F3-((F2^2)/F1)))
sum(alpha0 <= T+range & alpha0 >= T-range)/n

boxplot(sqrt(1/(F3-((F2^2)/F1))))
abline(h=sd(T), col="red", lty=2)

#1, 0.01, 66%
#1, 0.03, 76.6%
#1, 0.05, 56%
#2, 0.01, 77.2%
#2, 0.03, 78.8%
#2, 0.05, 82.2%
#3, 0.01, 81%
#3, 0.03, 74.4%
#3, 0.05, 12%


#Number of correct IC for the estimate of rho using the partial information matrix
range = 1.96 / sqrt(F1)
sum(rho0 <= Intens+range & rho0 >= Intens-range)/n

boxplot(sqrt(1/(F1)))
abline(h=sd(Intens), col="red", lty=2)

#Number of correct IC for the estimate of rho using the full information matrix
range = 1.96 * sqrt(1/(F1-((F2^2)/F3)))
sum(rho0 <= Intens+range & rho0 >= Intens-range)/n

boxplot(sqrt(1/(F1-((F2^2)/F3))))
abline(h=sd(Intens), col="red", lty=2)
