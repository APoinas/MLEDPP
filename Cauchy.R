library(spatstat)
library(stats)
library(foreach)
library(doSNOW)
library(parallel)

#####Functions#####

Lker=function(rho,alpha,r,Nsum){
  return(rho*sum((2*rho*pi*alpha^2)^(1:Nsum-1)/((1:Nsum)^2*(1+(r/(alpha*1:Nsum))^2)^(3/2))))
}

func <- function(rho,alpha,Nsum){
  if (2*rho*pi*alpha^2<1) {
    l=dim(M)[1]
    N=diag(l)
    for (i in 1:l){
      for (j in i:l){
        N[i,j]=Lker(rho,alpha,M[i,j],Nsum)
        N[j,i]=N[i,j]
      }
    }
    a=-rho*sum((2*rho*pi*alpha^2)^(1:Nsum-1)/(1:Nsum)^3)+(1/vol)*determinant(N,logarithm=TRUE)$modulus[1]
    return(-a)}
  else {return(1000000)}
}

#####Plotting the approximated LL in an example#####

w <- owin(c(0,2), c(0,2))
X=dppCauchy(lambda=100,alpha=0.035,nu=0.5,d=2)
S=simulate(X,1,W=w)
l=2
vol=4
M=pmin(sqrt(outer(S$x-l,S$x,'-')^2+outer(S$y-l,S$y,'-')^2),sqrt(outer(S$x,S$x,'-')^2+outer(S$y-l,S$y,'-')^2),sqrt(outer(S$x+l,S$x,'-')^2+outer(S$y-l,S$y,'-')^2),sqrt(outer(S$x-l,S$x,'-')^2+outer(S$y,S$y,'-')^2),sqrt(outer(S$x,S$x,'-')^2+outer(S$y,S$y,'-')^2),sqrt(outer(S$x+l,S$x,'-')^2+outer(S$y,S$y,'-')^2),sqrt(outer(S$x-l,S$x,'-')^2+outer(S$y+l,S$y,'-')^2),sqrt(outer(S$x,S$x,'-')^2+outer(S$y+l,S$y,'-')^2),sqrt(outer(S$x+l,S$x,'-')^2+outer(S$y+l,S$y,'-')^2))
rho=S$n/vol
T=seq(0.001,1/sqrt(2*pi*S$n/vol),0.001)
L=c()
for (i in T){
  L=c(L,func(rho=rho,alpha=i,Nsum=10))
}
plot(T,-L)

#####Generating Data#####

n=250
w <- owin(c(0,2), c(0,2))
X=dppCauchy(lambda=100,alpha=0.005,nu=0.5,d=2)
cl<-makeCluster(6)
registerDoSNOW(cl)
S=foreach(a=1:n, .combine="c") %dopar% {simulate(X,2,W=w)}
stopCluster(cl)
names(S)=lapply(1:(2*n), function(x) paste("Simulation ",toString(x),sep=""))
saveRDS(S,"500_Cauchy_100_5e-3_[0_2].dat")

#####Estimating alpha#####

#Initialisation
rho0=100
alpha0=0.035
S_length=3
n=500
Simu=readRDS(paste("Cauchy_DPP/500_Cauchy_",toString(rho0),"_",toString(alpha0*1000),"e-3_[0_",toString(S_length),"].dat",sep=""))
w=Window(Simu[[1]])
vol=area(w)
Ncores=4

#MCE with g
T2=unlist(unname(mclapply(Simu[1:n], function(x) dppm(x~1, dppCauchy(nu=0.5), method="mincon",statistic="pcf", rmin=0.01, q=1/2)$fitted$fixedpar$alpha, mc.cores=Ncores)))

#MCE with K
T2b=unlist(unname(mclapply(Simu[1:n], function(x) dppm(x~1, dppCauchy(nu=0.5), method="mincon",statistic="K", rmin=0.01, q=1/2)$fitted$fixedpar$alpha, mc.cores=Ncores)))

#Approwimated MLE with edge effect correction
cl=makeCluster(Ncores)
registerDoSNOW(cl)
templist=foreach(a=1:n, .combine="cbind") %dopar% {
  S=Simu[[a]]
  l=S_length
  M=pmin(sqrt(outer(S$x-l,S$x,'-')^2+outer(S$y-l,S$y,'-')^2),sqrt(outer(S$x,S$x,'-')^2+outer(S$y-l,S$y,'-')^2),sqrt(outer(S$x+l,S$x,'-')^2+outer(S$y-l,S$y,'-')^2),sqrt(outer(S$x-l,S$x,'-')^2+outer(S$y,S$y,'-')^2),sqrt(outer(S$x,S$x,'-')^2+outer(S$y,S$y,'-')^2),sqrt(outer(S$x+l,S$x,'-')^2+outer(S$y,S$y,'-')^2),sqrt(outer(S$x-l,S$x,'-')^2+outer(S$y+l,S$y,'-')^2),sqrt(outer(S$x,S$x,'-')^2+outer(S$y+l,S$y,'-')^2),sqrt(outer(S$x+l,S$x,'-')^2+outer(S$y+l,S$y,'-')^2))
  optim(par=1/(2*sqrt(2*pi*S$n/vol)),fn=func,rho=S$n/vol,Nsum=50,lower=0.0001,upper=1/sqrt(2*pi*S$n/vol),method="Brent")$par}
stopCluster(cl)
T=unname(templist[1,])

#Presenting and saving the results
boxplot(T-alpha0,T2-alpha0,T2b-alpha0)
10^4*c(mean((T-alpha0)^2),mean((T2-alpha0)^2),mean((T2b-alpha0)^2))

A=rbind(T,T2,T2b)
row.names(A)=c("MLE Tore","MCE with g","MCE with K")
saveRDS(A,paste("ToResult_",toString(n),"_Cauchy_",toString(rho0),"_",toString(1000*alpha0),"e-3_[0_",toString(S_length),"].dat",sep=""))

#####Showcase Results#####

cmain=2.5
caxis=1.5
thiccc=1.5
par(mfrow=c(3,3),mar=c(1, 2.5, 2, 1))
for (wind in c(1,2,3)){
  for (alpha in c(5,20,35)){
    filename=paste("Results/Cauchy_DPP/Result_500_Cauchy_100_",alpha,"e-3_[0_",wind,"].dat",sep="")
    titletext=bquote(alpha*"*="*.(alpha/1000)*", W=[0,"*.(wind)*"]"~"")
    T=t(readRDS(filename))
    boxplot(cbind(as.numeric(T[,1])-alpha/1000,as.numeric(T[,2])-alpha/1000,as.numeric(T[,3])-alpha/1000),yaxt="n",cex.main=cmain,main=titletext,xaxt="n",lwd=thiccc)
    abline(h = 0, col = "red",lty=2,lwd=thiccc)
    axis(2,cex.axis=caxis)
    print(c(alpha,wind,10^4*mean((as.numeric(T[,1])-alpha/1000)^2),10^4*mean((as.numeric(T[,2])-alpha/1000)^2),10^4*mean((as.numeric(T[,3])-alpha/1000)^2)))
  }}

#1000 x 1000