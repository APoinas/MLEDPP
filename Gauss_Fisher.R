library(spatstat)
library(tictoc)
library(stats)
library(foreach)
library(doSNOW)
library(parallel)

#####Functions#####

#L0
L0=function(rho,alpha,r,Nsum){
  N = 1:Nsum
  if (is.null(dim(r))){return(rho*sum((rho*pi*alpha^2)^(N-1)*exp(-r^2/(alpha^2*N))/N))}
  d = dim(r)[1]
  bigM = array(rep(0, d*d), dim=c(d,d))
  for (i in N){
    bigM = bigM + (rho*pi*alpha^2)^(i-1)*exp(-r^2/(alpha^2*i))/i
  }
  return (rho * bigM)
}

#Derivative of L0 with respect to alpha
DalpL0=function(rho,alpha,r,Nsum){
  N = 1:Nsum
  if (is.null(dim(r))){return(2*rho*sum((rho*pi)^(N-1) * ((N-1) * alpha^(2*N-3) + r^2 * alpha^(2*N-5)/N) * exp(-r^2/(alpha^2*N))/N))}
  d = dim(r)[1]
  bigM = array(rep(0, d*d), dim=c(d,d))
  for (i in N){
    bigM = bigM + (rho*pi)^(i-1)*((i-1) * alpha^(2*i-3) + r^2 * alpha^(2*i-5)/i)*exp(-r^2/(alpha^2*i))/i
  }
  return (2 * rho * bigM)
}

#2nd derivative of L0 with respect to alpha
D2alpL0=function(rho,alpha,r,Nsum){
  N = 1:Nsum
  if (is.null(dim(r))){return(2*rho*sum((rho*pi)^(N-1) * ((N-1) * (2*N-3) * alpha^(2*N-4) + r^2 * (4*N-7) * alpha^(2*N-6)/N + 2*r^4 * alpha^(2*N-8)/ (N^2)) * exp(-r^2/(alpha^2*N))/N))}
  d = dim(r)[1]
  bigM = array(rep(0, d*d), dim=c(d,d))
  for (i in N){
    bigM = bigM + (rho*pi)^(i-1) * ((i-1) * (2*i-3) * alpha^(2*i-4) + r^2 * (4*i-7) * alpha^(2*i-6)/i + 2*r^4 * alpha^(2*i-8)/ (i^2))*exp(-r^2/(alpha^2*i))/i
  }
  return (2 * rho * bigM)
}

#Derivative of L0 with respect to rho
DrhoL0=function(rho,alpha,r,Nsum){
  N = 1:Nsum
  if (is.null(dim(r))){return(sum((rho*pi*alpha^2)^(N-1)*exp(-r^2/(alpha^2*N))))}
  d = dim(r)[1]
  bigM = array(rep(0, d*d), dim=c(d,d))
  for (i in N){
    bigM = bigM + (rho*pi*alpha^2)^(i-1)*exp(-r^2/(alpha^2*i))
  }
  return (bigM)
}

#2nd derivative of L0 with respect to rho
D2rhoL0=function(rho,alpha,r,Nsum){
  N = 1:Nsum
  if (is.null(dim(r))){return(sum(rho^(N-2) * (pi*alpha^2)^(N-1) * (N-1) * exp(-r^2/(alpha^2*N))))}
  d = dim(r)[1]
  bigM = array(rep(0, d*d), dim=c(d,d))
  for (i in N){
    bigM = bigM + rho^(i-2) * (pi*alpha^2)^(i-1) * (i-1) * exp(-r^2/(alpha^2*i))
  }
  return (bigM)
}

#Derivative of L0 with respect to alpha and rho
DalpDrhoL0=function(rho,alpha,r,Nsum){
  N = 1:Nsum
  if (is.null(dim(r))){return(2*sum((rho*pi)^(N-1) * ((N-1) * alpha^(2*N-3) + r^2 * alpha^(2*N-5)/N) * exp(-r^2/(alpha^2*N))))}
  d = dim(r)[1]
  bigM = array(rep(0, d*d), dim=c(d,d))
  for (i in N){
    bigM = bigM + (rho*pi)^(i-1)*((i-1) * alpha^(2*i-3) + r^2 * alpha^(2*i-5)/i)*exp(-r^2/(alpha^2*i))
  }
  return (2 * bigM)
}

#Integral part of the log-likelihood
Integ = function(rho,alpha,Nsum){
  N = 1:Nsum
  return(-rho * sum( (rho*pi*alpha^2)^(N-1) / (N^2)))
}

#Derivative of the integral part of the log-likelihood with respect to alpha
DalpInteg = function(rho,alpha,Nsum){
  N = 1:Nsum
  return(-rho * sum( (rho*pi)^(N-1) * alpha^(2*N - 3) * (2*N - 2) / (N^2)))
}

#2nd derivative of the integral part of the log-likelihood with respect to alpha
D2alpInteg = function(rho,alpha,Nsum){
  N = 1:Nsum
  return(-rho * sum( (rho*pi)^(N-1) * alpha^(2*N - 4) * (2*N - 3) * (2*N - 2) / (N^2)))
}

#Derivative of the integral part of the log-likelihood with respect to rho
DrhoInteg = function(rho,alpha,Nsum){
  return(log(1-rho*pi*alpha^2)/(rho*pi*alpha^2))
}

#2nd derivative of the integral part of the log-likelihood with respect to rho
D2rhoInteg = function(rho,alpha,Nsum){
  return(-(1/(1-rho*pi*alpha^2) + log(1-rho*pi*alpha^2)/(rho*pi*alpha^2))/rho)
}

#Derivative of the integral part of the log-likelihood with respect to rho and alpha
DalpDrhoInteg = function(rho,alpha,Nsum){
  N = 1:Nsum
  return(-2 * sum( (rho*pi)^(N-1) * alpha^(2*N - 3) * (N - 1) / N))
}

#(minus) likelihood divided by the volume
func = function(rho,alpha,M,Nsum){
  if (rho*pi*alpha^2>1) {return(NA)}
  N=L0(rho,alpha,M,Nsum)
  a=Integ(rho,alpha,Nsum)+(1/vol)*determinant(N,logarithm=TRUE)$modulus[1]
  return(a)
}

Grad_LL = function(rho,alpha,M,Nsum){
  if (rho*pi*alpha^2>1) {return(NA)}
  N = L0(rho,alpha,M,Nsum)
  DalpN = DalpL0(rho,alpha,M,Nsum)
  DrhoN = DrhoL0(rho,alpha,M,Nsum)
  invN = solve(N)
  a = vol * DrhoInteg(rho,alpha,Nsum) + sum(DrhoN * invN)
  b = vol * DalpInteg(rho,alpha,Nsum) + sum(DalpN * invN)
  return(c(a,b))
}

MM = M[1:3,1:3]
func(100,0.03,MM,50)
fun = function(h){
  return((func(100+h,0.03,M,50)-func(100,0.03,M,50))/h)
}
fun2 = function(h){
  return((func(100,0.03+h,M,50)-func(100,0.03,M,50))/h)
}
fun(0.00000000001)
fun2(0.00000000001)
Grad_LL(100,0.03,M,50)

#Fisher information matrix
Fisher_Info = function(rho,alpha,M,Nsum){
  if (rho*pi*alpha^2>1) {return(NA)}
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
  return(array(c(-a,-b,-b,-c), dim=c(2,2)))
}

#################################################################################################
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

#1, 0.01, 88.2%
#1, 0.03, 89.6%
#1, 0.05, 92%
#2, 0.01, 88.6%
#2, 0.03, 94.2%
#2, 0.05, 92.6%
#3, 0.01, 93.2%
#3, 0.03, 93.2%
#3, 0.05, 92.8%

#Number of correct IC for the estimate of rho using the partial information matrix
range = 1.96 / sqrt(F1)
sum(rho0 <= Intens+range & rho0 >= Intens-range)/n

boxplot(sqrt(1/(F1)))
abline(h=sd(Intens), col="red", lty=2)

#1, 0.01, 95%
#1, 0.03, 95.6%
#1, 0.05, 91.2%
#2, 0.01, 93.8%
#2, 0.03, 93.4%
#2, 0.05, 94.6%
#3, 0.01, 95.6%
#3, 0.03, 92%
#3, 0.05, 94.6%

#Number of correct IC for the estimate of rho using the full information matrix
range = 1.96 * sqrt(1/(F1-((F2^2)/F3)))
sum(rho0 <= Intens+range & rho0 >= Intens-range)/n

boxplot(sqrt(1/(F1-((F2^2)/F3))))
abline(h=sd(Intens), col="red", lty=2)

#1, 0.01, 95%
#1, 0.03, 95.6%
#1, 0.05, 93%
#2, 0.01, 93.8%
#2, 0.03, 93.4%
#2, 0.05, 94.8%
#3, 0.01, 95.6%
#3, 0.03, 92%
#3, 0.05, 94.6%

FF = Intens*(1 - pi * Intens * T^2 / 2)/vol
range = 1.96 * sqrt(FF)
sum(rho0 <= Intens+range & rho0 >= Intens-range)/n

boxplot(sqrt(1/(F1-((F2^2)/F3))))
abline(h=sd(Intens), col="red", lty=2)

#######################################Sandwich method###########################################
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
  Grad = Grad_LL(S$n/vol, res, M, 50)
  Fish = Fisher_Info(S$n/vol, res, M, 50)
  c(res, S$n/vol, Grad[1], Grad[2], Fish[1,1], Fish[1,2], Fish[2,2])}
stopCluster(cl)
T = unname(templist[1,])
Intens = unname(templist[2,])
G1 = unname(templist[3,])
G2 = unname(templist[4,])
F1 = unname(templist[5,])
F2 = unname(templist[6,])
F3 = unname(templist[7,])

range = c()
for (i in 1:n){
  G = array(c(G1[i],G2[i]), dim=c(2,1))
  H = array(c(F1[i],F2[i],F2[i],F3[i]), dim=c(2,2))
  Sandwich = solve(H) %*% G %*% t(G) %*% solve(H)
}

##############################################R-shaped#########################################
newM.nonrectangular_edgecorrection = function(Nmax, M, bdist){
  l = dim(M)[1]
  dim(Nmax) = c(l,l)
  tempbool = Nmax > 0.001*max(Nmax)
  rmax = max(M[tempbool])
  
  noborder = bdist>rmax
  if(sum(noborder) == 0) noborder = which.max(bdist)
  
  distpossible = M[noborder,]
  distpossible = matrix(distpossible, ncol=l)
  
  nb.replace.possible = apply(distpossible < rmax, 1, sum)
  
  index = order(bdist)
  indexmax = which(sort(bdist) > rmax)
  if(length(indexmax) == 0){imax=l}else{imax = min(indexmax)}
  
  Mnew = M
  for(i in 1:(imax-1)){
    distpossible.i = distpossible[distpossible > bdist[index[i]]]
    nb.replace.i = sample(nb.replace.possible, 1) - sum(M[index[i],] < rmax)
    if(nb.replace.i > 0){
      for(j in (i+1):min(l,i+nb.replace.i)){
        if(Mnew[index[i], index[j]] > rmax){
          Mnew[index[i], index[j]] = sample(distpossible.i, 1)
          Mnew[index[j], index[i]] = Mnew[index[i], index[j]]
        }
      }}
  }
  return(Mnew)
}

Lker = function(rho, alpha, r, Max_trunc){
  Nsum = min(Max_trunc, 1 + abs(ceiling(log(0.0001)/log(rho*pi*alpha^2))))
  d = dim(r)[1]
  bigM = array(rep(0, d*d), dim=c(d,d))
  for (i in 1:Nsum){
    bigM = bigM + (rho*pi*alpha^2)^(i-1)*exp(-r^2/(alpha^2*i))/i
  }
  return (rho * bigM)
}

###########################################
rho0 = 100
alpha0 = 0.01
nb_core = 7
Simu = readRDS(paste("../Git_MLEDPP/MLEDPP/datasets/500_Gauss_",toString(alpha0*100),"e-2_letterR.dat",sep=""))
n = 500
w = Window(Simu[[1]])
vol = area(w)

#Setting up the progress bar
pb = txtProgressBar(max = n, style = 3)
progress = function(n) setTxtProgressBar(pb, n)
opts = list(progress = progress)

#Loading the values of the estimates of alpha
T2 = readRDS(paste("../Git_MLEDPP/MLEDPP/results/Result_Gauss_",toString(alpha0*100),"e-2_letterR.dat",sep=""))
old_result = T2$MLE_corrected

#Computing the corresponding Fisher information
cl = makeCluster(nb_core)
clusterExport(cl, c('pairdist', 'bdist.points', 'Lker', 'newM.nonrectangular_edgecorrection'))
registerDoSNOW(cl)
templist = foreach(a=1:n, .combine="cbind", .options.snow = opts) %dopar% {
  S = Simu[[a]]
  bdist = bdist.points(S)
  alphamax = sqrt(vol/(pi*S$n))
  
  M = as.matrix(dist(t(rbind(S$x, S$y)), method = "euclidean", diag = TRUE, upper = TRUE, p = 2))
  Nmax = Lker(S$n/vol, alphamax, M, 50)
  M = newM.nonrectangular_edgecorrection(Nmax, M, bdist)
  
  res = old_result[a]
  Fish = Fisher_Info(S$n/vol, res, M, 50)
  c(res, S$n/vol, Fish[1,1], Fish[1,2], Fish[2,2])}
stopCluster(cl)
T = unname(templist[1,])
Intens = unname(templist[2,])
F1 = unname(templist[3,])
F2 = unname(templist[4,])
F3 = unname(templist[5,])

sum(F3<=0)/n
sum(pi*Intens*T^2>0.99)/n

#0.01, 0%
#0.03, 0%
#0.05, 2.6%

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

#0.01, 87.4%
#0.03, 90%
#0.05, 57.8%

#Number of correct IC for the estimate of alpha using the full information matrix
range = 1.96 * sqrt(1/(F3-((F2^2)/F1)))
sum(alpha0 <= T+range & alpha0 >= T-range)/n

boxplot(sqrt(1/(F3-((F2^2)/F1))))
abline(h=sd(T), col="red", lty=2)

#0.01, 87.4%
#0.03, 90%
#0.05, 63.2%



###########################################
rho0 = 100
alpha0 = 0.01
nb_core = 7
Simu = readRDS(paste("../Git_MLEDPP/MLEDPP/datasets/500_Gauss_",toString(alpha0*100),"e-2_letterR.dat",sep=""))
n = 500
w = Window(Simu[[1]])
vol = area(w)

#Setting up the progress bar
pb = txtProgressBar(max = n, style = 3)
progress = function(n) setTxtProgressBar(pb, n)
opts = list(progress = progress)

#Loading the values of the estimates of alpha
T2 = readRDS(paste("../Git_MLEDPP/MLEDPP/results/Result_Gauss_",toString(alpha0*100),"e-2_letterR.dat",sep=""))
old_result = T2$MLE_corrected

#Computing the corresponding Fisher information
cl = makeCluster(nb_core)
registerDoSNOW(cl)
templist = foreach(a=1:n, .combine="cbind", .options.snow = opts) %dopar% {
  S = Simu[[a]]
  M = as.matrix(dist(t(rbind(S$x, S$y)), method = "euclidean", diag = TRUE, upper = TRUE, p = 2))
  Fish = Fisher_Info(S$n/vol, res, M, 50)
  c(res, S$n/vol, Fish[1,1], Fish[1,2], Fish[2,2])}
stopCluster(cl)
T = unname(templist[1,])
Intens = unname(templist[2,])
F1 = unname(templist[3,])
F2 = unname(templist[4,])
F3 = unname(templist[5,])

sum(F3<=0)/n
sum(pi*Intens*T^2>0.99)/n

#0.01, 0%
#0.03, 
#0.05, 

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

#0.01, 
#0.03, 
#0.05, 

#Number of correct IC for the estimate of alpha using the full information matrix
range = 1.96 * sqrt(1/(F3-((F2^2)/F1)))
sum(alpha0 <= T+range & alpha0 >= T-range)/n

boxplot(sqrt(1/(F3-((F2^2)/F1))))
abline(h=sd(T), col="red", lty=2)

#0.01, 
#0.03,
#0.05,

