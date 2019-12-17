library(spatstat)
library(tictoc)
library(stats)
library(foreach)
library(doSNOW)
library(parallel)

#####Functions#####

Lker=function(rho,alpha,r,Nsum){
  return(rho*sum((rho*pi*alpha^2)^(1:Nsum-1)*exp(-r^2/(alpha^2*(1:Nsum)))/(1:Nsum)))
}

Lker_vec=Vectorize(Lker)

func <- function(rho,alpha,Nsum){
  if (rho*pi*alpha^2<1) {
    l=dim(M)[1]
    N=Lker_vec(rho,alpha,M,Nsum)
    dim(N)=c(l,l)
    a=integrate(function(r) 2*pi*r*log(1-rho*pi*alpha^2*exp(-(pi*alpha*r)^2)),0,Inf)$value+(1/vol)*determinant(N,logarithm=TRUE)$modulus[1]
    return(-a)}
  else {return(1000000)}
}

#####Generating Data#####

n=1000
w <- owin(c(0,3), c(0,3))
X=dppGauss(lambda=100,alpha=0.03,d=2)
cl<-makeCluster(48)
registerDoSNOW(cl)
S=foreach(a=1:n, .combine="c") %dopar% {simulate(X,2,W=w)}
stopCluster(cl)
names(S)=lapply(1:2000, function(x) paste("Simulation ",toString(x),sep=""))
saveRDS(S,"simu.dat")

#####Estimating alpha on the square#####

#Initialisation
rho0=100
alpha0=0.03
S_length=1
n=100
Simu=readRDS(paste("Gauss_DPP/2000_Gauss_",toString(rho0),"_",toString(alpha0*100),"e-2_[0_",toString(S_length),"].dat",sep=""))
Nsum=switch((100*alpha0+1)/2,20,30,50)
w=Window(Simu[[1]])
vol=area(w)
ncores=6

#MCE with g
T2=unlist(unname(mclapply(Simu[1:n], function(x) dppm(x~1, dppGauss, method="mincon",statistic="pcf", rmin=0.01, q=1/2)$fitted$fixedpar$alpha, mc.cores=ncores)))

#MCE with K
T2b=unlist(unname(mclapply(Simu[1:n], function(x) dppm(x~1, dppGauss, method="mincon",statistic="K", rmin=0.01, q=1/2)$fitted$fixedpar$alpha, mc.cores=ncores)))

#Approximated MLE without edge effect correction
cl=makeCluster(ncores)
registerDoSNOW(cl)
templist=foreach(a=1:n, .combine="cbind") %dopar% {
  S=Simu[[a]]
  M=as.matrix(dist(t(rbind(S$x,S$y)), method = "euclidean", diag = TRUE, upper = TRUE, p = 2))
  optim(par=sqrt(1/(2*pi*S$n/vol)),fn=func,Nsum=Nsum,rho=S$n/vol,lower=0.0001,upper=1/sqrt(pi*S$n/vol),method="Brent")$par}
stopCluster(cl)
T=unname(templist[1,])

#Approximated MLE with edge effect correction
cl=makeCluster(ncores)
registerDoSNOW(cl)
templist=foreach(a=1:n, .combine="cbind") %dopar% {
  S=Simu[[a]]
  l=S_length
  M=pmin(sqrt(outer(S$x-l,S$x,'-')^2+outer(S$y-l,S$y,'-')^2),sqrt(outer(S$x,S$x,'-')^2+outer(S$y-l,S$y,'-')^2),sqrt(outer(S$x+l,S$x,'-')^2+outer(S$y-l,S$y,'-')^2),sqrt(outer(S$x-l,S$x,'-')^2+outer(S$y,S$y,'-')^2),sqrt(outer(S$x,S$x,'-')^2+outer(S$y,S$y,'-')^2),sqrt(outer(S$x+l,S$x,'-')^2+outer(S$y,S$y,'-')^2),sqrt(outer(S$x-l,S$x,'-')^2+outer(S$y+l,S$y,'-')^2),sqrt(outer(S$x,S$x,'-')^2+outer(S$y+l,S$y,'-')^2),sqrt(outer(S$x+l,S$x,'-')^2+outer(S$y+l,S$y,'-')^2))
  optim(par=sqrt(1/(2*pi*S$n/vol)),fn=func,Nsum=Nsum,rho=S$n/vol,lower=0.0001,upper=1/sqrt(pi*S$n/vol),method="Brent")$par}
stopCluster(cl)
T1=unname(templist[1,])

#Presenting and saving the results
boxplot(T-alpha0,T1-alpha0,T2-alpha0,T2b-alpha0)
c(10^4*mean((T-alpha0)^2),10^4*mean((T1-alpha0)^2),10^4*mean((T2-alpha0)^2),10^4*mean((T2b-alpha0)^2))

A=rbind(T,T1,T2,T2b)
row.names(A)=c("MLE (not corrected)","MLE (corrected)","MCE with the pcf","MLE with K")
saveRDS(A,paste("Result_",toString(n),"_Gauss_",toString(rho0),"_",toString(100*alpha0),"e-2_[0_",toString(S_length),"].dat",sep=""))

#####Presenting the results nicely#####

par(mfrow=c(3,3),mar=c(1, 2.5, 2, 1))
cmain=2.5
caxis=1.5
thiccc=1.5
for (wind in c(1,2,3)){
  for (alpha in c(1,3,5)){
    filename=paste("Result_500_Gauss_100_",alpha,"e-2_[0_",wind,"].dat",sep="")
    titletext=bquote(alpha*"*="*.(alpha/100)*", W=[0,"*.(wind)*"]²"~"")
    T=t(readRDS(paste("results/",filename,sep="")))
    T1=t(readRDS(paste("results/To",filename,sep="")))
    boxplot(cbind(T1[2:501]-alpha/100,T[2:501,c(2,3,4)-(wind==3)]-alpha/100),yaxt="n",cex.main=cmain,main=titletext,xaxt="n",lwd=thiccc)
    print(c(alpha/100,wind,mean((T1[2:501]-alpha/100)^2),mean((T[2:501,2-(wind==3)]-alpha/100)^2),mean((T[2:501,3-(wind==3)]-alpha/100)^2),mean((T[2:501,4-(wind==3)]-alpha/100)^2)))
    abline(h = 0, col = "red",lty=2,lwd=thiccc)
    axis(2,cex.axis=caxis)}}

#1000 by 1000