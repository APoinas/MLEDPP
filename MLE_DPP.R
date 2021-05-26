library(spatstat)
library(stats)
source("LL_derivatives.R")

############################### General MLE function ################################

MLEDPP = function(ppp, DPPfamily, startpar=NULL, sigma=NULL, edgecorr=FALSE, Trunc=50){
  if (!is.ppp(ppp)){
    stop('First argument is not of the ppp type.')}
  if (!(DPPfamily %in% c("Gauss", "Bessel", "Cauchy", "Whittle-Matern", "WM"))){
    stop('Second argument should be either Gauss, Bessel, Cauchy or Whittle-Matern.')}
  if ((DPPfamily == "WM" | DPPfamily == "Whittle-Matern") & is.null(sigma)){
    stop('The sigma parameter has to be specified for the Whittle-Matern family')}
  if(!is.null(sigma) & DPPfamily %in% c("Gauss", "Bessel", "Cauchy")){stop('The sigma parameter can only be specified for the Whittle-Matern family')}
  
  Win = ppp$window #Observation window
  vol = area(Win) #Area of the observation window
  
  #Estimation of the intensity
  rho_est = ppp$n/vol
    
  #Generating the distance matrix
  if (!edgecorr){
    M = as.matrix(dist(t(rbind(ppp$x, ppp$y)), method = "euclidean", diag = TRUE, upper = TRUE, p = 2))}
  else {
    if (Win$type == "rectangle") {M=pairdist(ppp, periodic=TRUE)}
    else {  
      warning("Edge correction is experimental for non rectangular windows.")
      bdist = bdist.points(ppp)
      M = as.matrix(dist(t(rbind(ppp$x, ppp$y)), method = "euclidean", diag = TRUE, upper = TRUE, p = 2)) #M is modified below, depending on the family
    }
  }
  
  #Setting up the log-likelihood function
  if (DPPfamily == "Gauss"){
    DPPfam = dppGauss(lambda=0, alpha=0, d=2) #DPP family
    alphamax = sqrt(1/(pi*rho_est)) #Max value of alpha
    
    #Function corresponding to L_0^{rho,alpha}(r). Nsum is the truncation parameter.
    Lker = function(rho, alpha, r, Max_trunc){
      Nsum = min(Max_trunc, 1 + abs(ceiling(log(0.0001)/log(rho*pi*alpha^2))))
      d = dim(r)[1]
      bigM = array(rep(0, d*d), dim=c(d,d))
      for (i in 1:Nsum){
        bigM = bigM + (rho*pi*alpha^2)^(i-1)*exp(-r^2/(alpha^2*i))/i
      }
      return (rho * bigM)
    }
    
    
    if(edgecorr & Win$type != "rectangle"){
      Nmax = Lker(rho_est, alphamax, M, Max_trunc=Trunc)
      M = newM.nonrectangular_edgecorrection(Nmax, M, bdist)
    }
      
    #Approximate log-likelihood of Gauss-type DPPs
    Log_Likelihood = function(rho,alpha,M,vol,Max_trunc){
      if (rho*pi*alpha^2 >= 1) {return(NA)}
      N = Lker(rho,alpha,M,Max_trunc)
      temp = determinant(N,logarithm=TRUE)
      if(temp$sign == -1){return(NA)} #Returns NA if the determinant in the stochastic part of the log-likelihood is negative.
      a = integrate(function(r) 2*pi*r*log(1-rho*pi*alpha^2*exp(-(pi*alpha*r)^2)),0,Inf)$value+(1/vol)*temp$modulus[1]
      return(-a)} #The opposite of the log-likelihood is returned since the optimize function looks for minima.
    }
  
  if (DPPfamily == "Cauchy"){
    DPPfam = dppCauchy(lambda=0, alpha=0, nu=0.5, d=2) #DPP family
    alphamax = sqrt(1/(2*pi*rho_est)) #Max value of alpha
    
    #Function corresponding to L_0^{rho,alpha}(r). Nsum is the truncation parameter.
    Lker = function(rho, alpha, r, Max_trunc){
      Nsum = min(Max_trunc, 1 + abs(ceiling(log(0.0001)/log(2*rho*pi*alpha^2))))
      d = dim(r)[1]
      bigM = array(rep(0, d*d), dim=c(d,d))
      for (i in 1:Nsum){
        bigM = bigM + (2*rho*pi*alpha^2)^(i-1)/(i^2*(1+(r/(alpha*i))^2)^(3/2))
      }
      return (rho * bigM)
    }
    
    if(edgecorr & Win$type != "rectangle"){
      Nmax = Lker(rho_est, alphamax, M, Max_trunc=Trunc)
      M = newM.nonrectangular_edgecorrection(Nmax, M, bdist)
    }
    
    
    #Approximate log-likelihood of Cauchy-type DPPs
    Log_Likelihood = function(rho, alpha, M, vol, Max_trunc){
      if (2*rho*pi*alpha^2 >= 1) {return(NA)}
      Nsum = min(Max_trunc, 1 + abs(ceiling(log(0.0001)/log(2*rho*pi*alpha^2))))
      N = Lker(rho, alpha, M, Max_trunc)
      temp = determinant(N,logarithm=TRUE)
      if(temp$sign == -1){return(NA)} #Returns NA if the determinant in the stochastic part of the log-likelihood is negative.
      a=-rho*sum((2*rho*pi*alpha^2)^(1:Nsum-1)/(1:Nsum)^3)+(1/vol)*temp$modulus[1]
      return(-a) #The opposite of the log-likelihood is returned since the optimize function looks for minima.
      }
    }
  
  if (DPPfamily == "Bessel"){
    DPPfam = dppBessel(lambda=0, alpha=0, sigma=0, d=2) #DPP family
    alphamax = sqrt(1/(pi*rho_est)) #Max value of alpha
    
    if(edgecorr & Win$type != "rectangle"){
      Nmax = (besselJ(2*M/alphamax,1)/(M/alphamax))*rho_est/(1-rho_est*pi*alphamax^2)
      for (i in 1:sqrt(length(M))){
        Nmax[i,i] = rho_est/(1-rho_est*pi*alphamax^2)
      }
      M = newM.nonrectangular_edgecorrection(Nmax,M,bdist)
    }
    
    
    #Approximate log-likelihood of Bessel-type DPPs
    Log_Likelihood = function(rho,alpha,M,vol,Max_trunc=NULL){
      if (rho*pi*alpha^2 >= 1) {return(NA)}
      N = (besselJ(2*M/alpha,1)/(M/alpha))*rho/(1-rho*pi*alpha^2)
      for (i in 1:sqrt(length(M))){
        N[i,i] = rho/(1-rho*pi*alpha^2)
      }
      temp = determinant(N,logarithm=TRUE)
      if(temp$sign == -1){return(NA)} #Returns NA if the determinant in the stochastic part of the log-likelihood is negative.
      a = log(1-rho*pi*alpha^2)/(pi*alpha^2)+(1/vol)*temp$modulus[1]
      return(-a) #The opposite of the log-likelihood is returned since the optimize function looks for minima.
      }
    }
  
  if (DPPfamily == "WM" | DPPfamily == "Whittle-Matern"){
    DPPfam = dppMatern(lambda=0, alpha=0, nu=sigma, d=2) #DPP family
    alphamax = sqrt(1/(pi*rho_est*4*sigma)) #Max value of alpha
    
    #Function corresponding to L_0^{rho,alpha}(r). Nsum is the truncation parameter.
    Lker = function(rho,alpha,r,Max_trunc){
      Nsum = min(Max_trunc, 1 + abs(ceiling(log(0.0001)/log(4*rho*pi*sigma*alpha^2))))
      d = dim(r)[1]
      bigM = array(rep(0, d*d), dim=c(d,d))
      for (i in 1:Nsum){
        bigM = bigM + (4*rho*pi*alpha^2*sigma)^(i-1)/gamma(i*(sigma+1))*exp((i*(sigma+1)-1)*log(r/(2*alpha))-r/alpha)*besselK(r/alpha,(sigma+1)*i-1,expon.scaled = TRUE)
      }
      return (2 * rho * sigma * bigM)
    }
    
    if(edgecorr & Win$type != "rectangle"){
      Nmax = Lker(rho_est, alphamax, M, Max_trunc=Trunc)
      M = newM.nonrectangular_edgecorrection(Nmax, M, bdist)
    }
    
    
    #Approximate log-likelihood of Whittle-Matern-type DPPs
    Log_Likelihood = function(rho, alpha, M, vol, Max_trunc=NULL){
      if (4*rho*pi*sigma*alpha^2 >= 1) {return(NA)}
      Nsum = min(Max_trunc, 1 + abs(ceiling(log(0.0001)/log(4*rho*pi*sigma*alpha^2))))
      N = Lker(rho, alpha, M, Max_trunc)
      diagN = rho*sigma*sum((4*rho*pi*alpha^2*sigma)^(1:Nsum-1)/((sigma+1)*(1:Nsum)-1))
      diag(N) = diagN
      temp = determinant(N, logarithm=TRUE)
      if(temp$sign == -1){return(NA)} #Returns NA if the determinant in the stochastic part of the log-likelihood is negative.
      else{
        a = integrate(function(r) 2*pi*r*log(1-rho*sigma*4*pi*alpha^2/(1+4*pi^2*alpha^2*r^2)^(sigma+1)),0,Inf)$value+(1/vol)*temp$modulus[1]
        return(-a) #The opposite of the log-likelihood is returned since the optimize function looks for minima.
        }
      }
    }
  
  if (is.null(startpar)){startpar = alphamax/2} #By default, the starting parameter is half the maximum possible value  of alpha
  if (startpar > alphamax | startpar < 0){stop("Starting paramater not in range")}
  
  #Estimation of alpha by the maximum of the approximate log-likelihood
  alpha_est = optim(par=startpar, fn=Log_Likelihood, M=M, vol=vol, Max_trunc=Trunc, rho=rho_est, lower=alphamax/1000, upper=alphamax, method="Brent")$par
  
  DPPfam$fixedpar$lambda = rho_est
  DPPfam$fixedpar$alpha = alpha_est
  
  return(DPPfam)
}

##### (Experimental) Edge correction for non-rectangular windows #####

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

##### (Experimental) Variance estimation using Fisher's information matrix #####

Fisher_Info = function(ppp, DPPfamily, alpha_est, edgecorr=FALSE, Max_Trunc=50){
  if (DPPfamily == "WM" | DPPfamily == "Whittle-Matern"){
    stop('Fisher information not implemented for the Whittle-Matern family')}

  #Estimated parameters
  vol = area(ppp$window)
  rho = ppp$n / vol
  alpha = alpha_est
  
  #Generating the distance matrix
  if (!edgecorr){
    M = as.matrix(dist(t(rbind(ppp$x, ppp$y)), method = "euclidean", diag = TRUE, upper = TRUE, p = 2))}
  else {
    if (ppp$window$type == "rectangle") {M = pairdist(ppp, periodic=TRUE)}
    else {  
      stop('Edge correction for the Fisher information matrix not implemented for non-rectangular windows.')
    }
  }
  
  #Loading the various quantities used in the derivatives of the log-likelihood
  N = L0(rho, alpha, DPPfamily, M, Max_Trunc)
  DalpN = DalpL0(rho, alpha, DPPfamily, M, Max_Trunc)
  DrhoN = DrhoL0(rho, alpha, DPPfamily, M, Max_Trunc)
  D2alpN = D2alpL0(rho, alpha, DPPfamily, M, Max_Trunc)
  D2rhoN = D2rhoL0(rho, alpha, DPPfamily, M, Max_Trunc)
  DalpDrhoN = DalpDrhoL0(rho, alpha, DPPfamily, M, Max_Trunc)
  invN = solve(N)
  tempmat_alpha = DalpN %*% invN
  tempmat_rho = DrhoN %*% invN
  
  #2nd derivative of the log-likelihood with respect to rho
  D2rhoLL = vol * D2rhoInteg(rho, alpha, DPPfamily, Max_Trunc) + sum(D2rhoN * invN) - sum(tempmat_rho * t(tempmat_rho))
  #Derivative of the log-likelihood with respect to rho and alpha
  DalpDrhoLL = vol * DalpDrhoInteg(rho, alpha, DPPfamily, Max_Trunc) + sum(DalpDrhoN * invN) - sum(tempmat_alpha * t(tempmat_rho))
  #2nd derivative of the log-likelihood with respect to alpha
  D2alpLL = vol * D2alpInteg(rho, alpha, DPPfamily, Max_Trunc) + sum(D2alpN * invN) - sum(tempmat_alpha * t(tempmat_alpha))
  #Fisher information matrix
  return(array(c(-D2rhoLL,-DalpDrhoLL,-DalpDrhoLL,-D2alpLL), dim=c(2,2)))
}