library(spatstat)
library(tictoc)
library(stats)
library(foreach)
library(doSNOW)
library(parallel)

#####Functions#####

#L0
L0=function(rho,alpha,r,Nsum){
  d = dim(r)[1]
  bigM = array(rep(0, d*d), dim=c(d,d))
  for (i in 1:Nsum){
    bigM = bigM + (2*rho*pi*alpha^2)^(i-1)/(i^2*(1+(r/(alpha*i))^2)^(3/2))
  }
  return (rho * bigM)
}

#Derivative of L0 with respect to alpha
DalpL0=function(rho,alpha,r,Nsum){
  d = dim(r)[1]
  bigM = array(rep(0, d*d), dim=c(d,d))
  for (i in 1:Nsum){
    bigM = bigM + (2*rho*pi)^(i-1) * alpha^(2*i-3) * (2*i+1-3/(1 + (r / (i * alpha))^2))    /(i^2*(1+(r/(alpha*i))^2)^(3/2))
  }
  return (rho * bigM)
}

#2nd derivative of L0 with respect to alpha
D2alpL0=function(rho,alpha,r,Nsum){
  d = dim(r)[1]
  bigM = array(rep(0, d*d), dim=c(d,d))
  for (i in 1:Nsum){
    bigM = bigM + ((2*rho*pi)^(i-1)/(i^2)) * ( ( (2*i-2) * alpha^(2*i-4) * (2*i-3/(1 + (r / (i * alpha))^2))    /((1+(r/(alpha*i))^2)^(3/2)) ) + ( 3*(r/i)^2 * alpha^(2*i-6) * (2*i-5/(1 + (r / (i * alpha))^2))    /((1+(r/(alpha*i))^2)^(5/2)) ) )
  }
  return (rho * bigM)
}

#Derivative of L0 with respect to rho
DrhoL0=function(rho,alpha,r,Nsum){
  d = dim(r)[1]
  bigM = array(rep(0, d*d), dim=c(d,d))
  for (i in 1:Nsum){
    bigM = bigM + (2*rho*pi*alpha^2)^(i-1)/(i*(1+(r/(alpha*i))^2)^(3/2))
  }
  return (bigM)
}

#2nd derivative of L0 with respect to rho
D2rhoL0=function(rho,alpha,r,Nsum){
  d = dim(r)[1]
  bigM = array(rep(0, d*d), dim=c(d,d))
  for (i in 1:Nsum){
    bigM = bigM + (2*pi*alpha^2)^(i-1) * rho^(i-2) * (i-1)/(i*(1+(r/(alpha*i))^2)^(3/2))
  }
  return (bigM)
}

#Derivative of L0 with respect to alpha and rho
DalpDrhoL0=function(rho,alpha,r,Nsum){
  d = dim(r)[1]
  bigM = array(rep(0, d*d), dim=c(d,d))
  for (i in 1:Nsum){
    bigM = bigM + (2*rho*pi)^(i-1) * alpha^(2*i-3) * (2*i-2+3/(1 + i^2 * alpha^2 / r^2))    /(i*(1+(r/(alpha*i))^2)^(3/2))
  }
  return (bigM)
}

#Integral part of the log-likelihood
Integ = function(rho,alpha,Nsum){
  N = 1:Nsum
  return(-rho * sum( (2*rho*pi*alpha^2)^(N-1) / (N^3)))
}

#Derivative of the integral part of the log-likelihood with respect to alpha
DalpInteg = function(rho,alpha,Nsum){
  N = 1:Nsum
  return(- 2 * rho * sum( (N-1) * (2*rho*pi)^(N-1) * alpha^(2*N-3) / (N^3)))
}

#2nd derivative of the integral part of the log-likelihood with respect to alpha
D2alpInteg = function(rho,alpha,Nsum){
  N = 1:Nsum
  return(- 2 * rho * sum( (N-1) * (2*N-3) * (2*rho*pi)^(N-1) * alpha^(2*N-4) / (N^3)))
}

#Derivative of the integral part of the log-likelihood with respect to rho
DrhoInteg = function(rho,alpha,Nsum){
  N = 1:Nsum
  return(- sum( (2*rho*pi*alpha^2)^(N-1) / (N^2)))
}

#2nd derivative of the integral part of the log-likelihood with respect to rho
D2rhoInteg = function(rho,alpha,Nsum){
  N = 1:Nsum
  return(- sum( (N-1) * (2*pi*alpha^2)^(N-1) * rho^(N-2) / (N^2)))
}

#Derivative of the integral part of the log-likelihood with respect to rho and alpha
DalpDrhoInteg = function(rho,alpha,Nsum){
  N = 1:Nsum
  return(- 2 * sum( (N-1) * (2*rho*pi)^(N-1) * alpha^(2*N-3) / (N^2)))
}

#Fisher information matrix
Fisher_Info = function(rho,alpha,Nsum){
  if (2*rho*pi*alpha^2>1) {return(NA)}
  N = L0(rho,alpha,M,Nsum)
  DalpN = DalpL0(rho,alpha,M,Nsum)
  DrhoN = DrhoL0(rho,alpha,M,Nsum)
  D2alpN = D2alpL0(rho,alpha,M,Nsum)
  D2rhoN = D2rhoL0(rho,alpha,M,Nsum)
  DalpDrhoN = DalpDrhoL0(rho,alpha,M,Nsum)
  invN = solve(N)
  tempmat_alpha = DalpN %*% invN
  tempmat_rho = DrhoN %*% invN
  a = vol * D2rhoInteg(rho,alpha,Nsum) + sum(D2rhoN * invN) - sum(tempmat_rho * t(tempmat_rho))
  b = vol * DalpDrhoInteg(rho,alpha,Nsum) + sum(DalpDrhoN * invN) - sum(tempmat_alpha * t(tempmat_rho))
  c = vol * D2alpInteg(rho,alpha,Nsum) + sum(D2alpN * invN) - sum(tempmat_alpha * t(tempmat_alpha))
  print(c(N[1,1], DalpN[1,1], DrhoN[1,1], D2alpN[1,1], D2rhoN[1,1], DalpDrhoN[1,1]))
  return(array(c(-a,-b,-b,-c), dim=c(2,2)))
}

##############################################################################################
rho0 = 100
alpha0 = 0.035
S_length = 3
nb_core = 7
Simu = readRDS(paste("datasets/500_Cauchy_",toString(rho0),"_",toString(alpha0*1000),"e-3_[0_",toString(S_length),"].dat",sep=""))
Nsum = 50
n = 500
w = Window(Simu[[1]])
vol = area(w)

#Setting up the progress bar
pb = txtProgressBar(max = n, style = 3)
progress = function(n) setTxtProgressBar(pb, n)
opts = list(progress = progress)

#Loading the values of the estimates of alpha
T2 = readRDS(paste("Version propre pour github/results/Result_Cauchy_",toString(alpha0*1000),"e-3_[0_",toString(S_length),"].dat",sep=""))
old_result = T2$MLE_corrected

#Computing the corresponding Fisher information
cl=makeCluster(nb_core)
clusterExport(cl, c('pairdist'))
registerDoSNOW(cl)
templist=foreach(a=1:n, .combine="cbind", .options.snow = opts) %dopar% {
  S = Simu[[a]]
  M = pairdist(S, periodic=TRUE)
  res = old_result[a]
  Fish = Fisher_Info(S$n/vol, res, 50)
  c(res, S$n/vol, Fish[1,1], Fish[1,2], Fish[2,2])}
stopCluster(cl)
T = unname(templist[1,])
Intens = unname(templist[2,])
F1 = unname(templist[3,])
F2 = unname(templist[4,])
F3 = unname(templist[5,])

all_data = data.frame(rho_est=Intens, alpha_est=T, Fisher_1_1=F1, Fisher_1_2=F2, Fisher_2_2=F3)
saveRDS(all_data ,paste("Version propre pour github/results/Fisher/Result_Cauchy_",toString(alpha0*1000),"e-3_[0_",toString(S_length),"].dat",sep=""))

sum(F3<=0)/n
sum(2*pi*Intens*T^2>0.99)/n

#1, 0.005, 2% / 2.6%
#1, 0.020, 3.6% / 6.2%
#1, 0.035, 25.4% / 46.2 %
#2, 0.005, 0% / 0%
#2, 0.020, 0% / 0%
#2, 0.035, 0% / 4.4%
#3, 0.005, 0% / 0%
#3, 0.020, 0% / 0%
#3, 0.035, 0% / 0.2%

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

#1, 0.005, 89.8%
#1, 0.020, 88.4%
#1, 0.035, 72%
#2, 0.005, 92.4%
#2, 0.020, 92.2%
#2, 0.035, 83.6%
#3, 0.005, 91.4%
#3, 0.020, 95%
#3, 0.035, 87.2%

#Number of correct IC for the estimate of alpha using the full information matrix
range = 1.96 * sqrt(1/(F3-((F2^2)/F1)))
sum(alpha0 <= T+range & alpha0 >= T-range)/n

boxplot(sqrt(1/(F3-((F2^2)/F1))))
abline(h=sd(T), col="red", lty=2)

#1, 0.005, 89.8%
#1, 0.020, 88.4%
#1, 0.035, 72%
#2, 0.005, 92.4%
#2, 0.020, 92.2%
#2, 0.035, 84%
#3, 0.005, 91.4%
#3, 0.020, 95%
#3, 0.035, 87.2%






#Computing the corresponding Fisher information
cl=makeCluster(nb_core)
clusterExport(cl, c('pairdist'))
registerDoSNOW(cl)
templist=foreach(a=1:n, .combine="cbind", .options.snow = opts) %dopar% {
  S = Simu[[a]]
  M = pairdist(S, periodic=TRUE)
  res = old_result[a]
  Fish = Fisher_Info(S$n/vol, res, 50)
  c(res, S$n/vol, Fish[1,1], Fish[1,2], Fish[2,2])}
stopCluster(cl)
T = unname(templist[1,])
Intens = unname(templist[2,])
F1 = unname(templist[3,])
F2 = unname(templist[4,])
F3 = unname(templist[5,])

