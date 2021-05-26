#L0
L0=function(rho, alpha, DPPfamily, M, Nsum){
  if (DPPfamily=="Gauss"){
    d = dim(M)[1]
    bigM = array(rep(0, d*d), dim=c(d,d))
    for (i in 1:Nsum){
      bigM = bigM + (rho*pi*alpha^2)^(i-1)*exp(-M^2/(alpha^2*i))/i
    }
    return (rho * bigM)}
  
  if (DPPfamily=="Cauchy"){
    d = dim(M)[1]
    bigM = array(rep(0, d*d), dim=c(d,d))
    for (i in 1:Nsum){
      bigM = bigM + (2*rho*pi*alpha^2)^(i-1)/(i^2*(1+(M/(alpha*i))^2)^(3/2))
    }
    return (rho * bigM)}
  
  if (DPPfamily=="Bessel"){
    N = (besselJ(2*M/alpha,1)/(M/alpha))*rho/(1-rho*pi*alpha^2)
    diag(N) = rho/(1-rho*pi*alpha^2)
    return(N)}
}

#Derivative of L0 with respect to alpha
DalpL0=function(rho, alpha, DPPfamily, M, Nsum){
  if (DPPfamily=="Gauss"){
    d = dim(M)[1]
    bigM = array(rep(0, d*d), dim=c(d,d))
    for (i in 1:Nsum){
      bigM = bigM + (rho*pi)^(i-1)*((i-1) * alpha^(2*i-3) + M^2 * alpha^(2*i-5)/i)*exp(-M^2/(alpha^2*i))/i
    }
    return (2 * rho * bigM)}
  
  if (DPPfamily=="Cauchy"){  
    d = dim(M)[1]
    bigM = array(rep(0, d*d), dim=c(d,d))
    for (i in 1:Nsum){
      bigM = bigM + (2*rho*pi)^(i-1) * alpha^(2*i-3) * (2*i+1-3/(1 + (M / (i * alpha))^2))    /(i^2*(1+(M/(alpha*i))^2)^(3/2))
    }
    return (rho * bigM)}
  
  if (DPPfamily=="Bessel"){  
    N = (besselJ(2*M/alpha, 1)/(M/alpha))* 2*pi*alpha*rho^2 /(1-rho*pi*alpha^2)^2 + besselJ(2*M/alpha, 2) * 2*rho /(alpha * (1-rho*pi*alpha^2))
    diag(N) = 2*pi*alpha*rho^2/(1-rho*pi*alpha^2)^2
    return(N)}
}

#2nd derivative of L0 with respect to alpha
D2alpL0=function(rho, alpha, DPPfamily, M, Nsum){
  if (DPPfamily=="Gauss"){
    d = dim(M)[1]
    bigM = array(rep(0, d*d), dim=c(d,d))
    for (i in 1:Nsum){
      bigM = bigM + (rho*pi)^(i-1) * ((i-1) * (2*i-3) * alpha^(2*i-4) + M^2 * (4*i-7) * alpha^(2*i-6)/i + 2*M^4 * alpha^(2*i-8)/ (i^2))*exp(-M^2/(alpha^2*i))/i
    }
    return (2 * rho * bigM)}
  
  if (DPPfamily=="Cauchy"){
    d = dim(M)[1]
    bigM = array(rep(0, d*d), dim=c(d,d))
    for (i in 1:Nsum){
      bigM = bigM + ((2*rho*pi)^(i-1)/(i^2)) * ( ( (2*i-2) * alpha^(2*i-4) * (2*i-3/(1 + (M / (i * alpha))^2))    /((1+(M/(alpha*i))^2)^(3/2)) ) + ( 3*(M/i)^2 * alpha^(2*i-6) * (2*i-5/(1 + (M / (i * alpha))^2))    /((1+(M/(alpha*i))^2)^(5/2)) ) )
    }
    return (rho * bigM)}
  
  if (DPPfamily=="Bessel"){
    N = besselJ(2*M/alpha, 1)* (2*pi*alpha*rho^2 * (1+ 3*pi*rho*alpha^2)/ (M * (1-rho*pi*alpha^2)^3) - 4 * rho * M /(alpha^3 * (1-pi*rho*alpha^2))) + besselJ(2*M/alpha, 2) * 2*rho*(1+3*pi*rho*alpha^2) /(alpha^2 * (1-rho*pi*alpha^2)^2)
    diag(N) = 2*pi*rho^2 * (1+ 3*pi*rho*alpha^2)/(1-rho*pi*alpha^2)^3
    return(N)}
}

#Derivative of L0 with respect to rho
DrhoL0=function(rho, alpha, DPPfamily, M, Nsum){
  if (DPPfamily=="Gauss"){
    d = dim(M)[1]
    bigM = array(rep(0, d*d), dim=c(d,d))
    for (i in 1:Nsum){
      bigM = bigM + (rho*pi*alpha^2)^(i-1)*exp(-M^2/(alpha^2*i))
    }
    return (bigM)}
  
  if (DPPfamily=="Cauchy"){
    d = dim(M)[1]
    bigM = array(rep(0, d*d), dim=c(d,d))
    for (i in 1:Nsum){
      bigM = bigM + (2*rho*pi*alpha^2)^(i-1)/(i*(1+(M/(alpha*i))^2)^(3/2))
    }
    return (bigM)}
  
  if (DPPfamily=="Bessel"){
    N = (besselJ(2*M/alpha,1)/(M/alpha))/(1-rho*pi*alpha^2)^2
    diag(N) = 1/(1-rho*pi*alpha^2)^2
    return(N)}
}

#2nd derivative of L0 with respect to rho
D2rhoL0=function(rho, alpha, DPPfamily, M, Nsum){
  if (DPPfamily=="Gauss"){
    d = dim(M)[1]
    bigM = array(rep(0, d*d), dim=c(d,d))
    for (i in 1:Nsum){
      bigM = bigM + rho^(i-2) * (pi*alpha^2)^(i-1) * (i-1) * exp(-M^2/(alpha^2*i))
    }
    return (bigM)}
  
  if (DPPfamily=="Cauchy"){
    d = dim(M)[1]
    bigM = array(rep(0, d*d), dim=c(d,d))
    for (i in 1:Nsum){
      bigM = bigM + (2*pi*alpha^2)^(i-1) * rho^(i-2) * (i-1)/(i*(1+(M/(alpha*i))^2)^(3/2))
    }
    return (bigM)}
  
  if (DPPfamily=="Bessel"){
    N = (besselJ(2*M/alpha,1)/(M/alpha)) * 2*pi*alpha^2/(1-rho*pi*alpha^2)^3
    diag(N) = 2*pi*alpha^2 /(1-rho*pi*alpha^2)^3
    return(N)}
}

#Derivative of L0 with respect to alpha and rho
DalpDrhoL0=function(rho, alpha, DPPfamily, M, Nsum){
  if (DPPfamily=="Gauss"){
    d = dim(M)[1]
    bigM = array(rep(0, d*d), dim=c(d,d))
    for (i in 1:Nsum){
      bigM = bigM + (rho*pi)^(i-1)*((i-1) * alpha^(2*i-3) + M^2 * alpha^(2*i-5)/i)*exp(-M^2/(alpha^2*i))
    }
    return (2 * bigM)}
  
  if (DPPfamily=="Cauchy"){
    d = dim(M)[1]
    bigM = array(rep(0, d*d), dim=c(d,d))
    for (i in 1:Nsum){
      bigM = bigM + (2*rho*pi)^(i-1) * alpha^(2*i-3) * (2*i-2+3/(1 + i^2 * alpha^2 / M^2))    /(i*(1+(M/(alpha*i))^2)^(3/2))
    }
    return (bigM)}
  
  if (DPPfamily=="Bessel"){
    N = (besselJ(2*M/alpha, 1)/(M/alpha))* 4*pi*alpha*rho /(1-rho*pi*alpha^2)^3 + besselJ(2*M/alpha, 2) * 2 /(alpha * (1-rho*pi*alpha^2)^2)
    diag(N) = 4*pi*alpha*rho /(1-rho*pi*alpha^2)^3
    return(N)}
}

#Integral part of the log-likelihood
Integ = function(rho, alpha, DPPfamily, Nsum){
  if (DPPfamily=="Gauss"){
    N = 1:Nsum
    return(-rho * sum( (rho*pi*alpha^2)^(N-1) / (N^2)))}
  
  if (DPPfamily=="Cauchy"){
    N = 1:Nsum
    return(-rho * sum( (2*rho*pi*alpha^2)^(N-1) / (N^3)))}
  
  if (DPPfamily=="Bessel"){return (log(1-rho*pi*alpha^2)/(pi*alpha^2))}
}

#Derivative of the integral part of the log-likelihood with respect to alpha
DalpInteg = function(rho, alpha, DPPfamily, Nsum){
  if (DPPfamily=="Gauss"){
    N = 1:Nsum
    return(-rho * sum( (rho*pi)^(N-1) * alpha^(2*N - 3) * (2*N - 2) / (N^2)))}
  
  if (DPPfamily=="Cauchy"){
    N = 1:Nsum
    return(- 2 * rho * sum( (N-1) * (2*rho*pi)^(N-1) * alpha^(2*N-3) / (N^3)))}
  
  if (DPPfamily=="Bessel"){return ((rho*pi*alpha^2/(1-rho*pi*alpha^2) + log(1-rho*pi*alpha^2))* -2/(pi*alpha^3))}
}

#2nd derivative of the integral part of the log-likelihood with respect to alpha
D2alpInteg = function(rho, alpha, DPPfamily, Nsum){
  if (DPPfamily=="Gauss"){
    N = 1:Nsum
    return(-rho * sum( (rho*pi)^(N-1) * alpha^(2*N - 4) * (2*N - 3) * (2*N - 2) / (N^2)))}
  
  if (DPPfamily=="Cauchy"){
    N = 1:Nsum
    return(- 2 * rho * sum( (N-1) * (2*N-3) * (2*rho*pi)^(N-1) * alpha^(2*N-4) / (N^3)))}
  
  if (DPPfamily=="Bessel"){return (2*rho*(3-5*pi*rho*alpha^2)/(alpha^2 * (1-rho*pi*alpha^2)^2) + log(1-rho*pi*alpha^2)* 6/(pi*alpha^4))}
}

#Derivative of the integral part of the log-likelihood with respect to rho
DrhoInteg = function(rho, alpha, DPPfamily, Nsum){
  if (DPPfamily=="Gauss"){
    return(log(1-rho*pi*alpha^2)/(rho*pi*alpha^2))}
  
  if (DPPfamily=="Cauchy"){
    N = 1:Nsum
    return(- sum( (2*rho*pi*alpha^2)^(N-1) / (N^2)))}
  
  if (DPPfamily=="Bessel"){return (-1/(1-rho*pi*alpha^2))}
}

#2nd derivative of the integral part of the log-likelihood with respect to rho
D2rhoInteg = function(rho, alpha, DPPfamily, Nsum){
  if (DPPfamily=="Gauss"){
    return(-(1/(1-rho*pi*alpha^2) + log(1-rho*pi*alpha^2)/(rho*pi*alpha^2))/rho)}
  
  if (DPPfamily=="Cauchy"){
    N = 1:Nsum
    return(- sum( (N-1) * (2*pi*alpha^2)^(N-1) * rho^(N-2) / (N^2)))}
  
  if (DPPfamily=="Bessel"){return (-pi*alpha^2/(1-rho*pi*alpha^2)^2)}
}

#Derivative of the integral part of the log-likelihood with respect to rho and alpha
DalpDrhoInteg = function(rho, alpha, DPPfamily, Nsum){
  if (DPPfamily=="Gauss"){
    N = 1:Nsum
    return(-2 * sum( (rho*pi)^(N-1) * alpha^(2*N - 3) * (N - 1) / N))}
  
  if (DPPfamily=="Cauchy"){
    N = 1:Nsum
    return(- 2 * sum( (N-1) * (2*rho*pi)^(N-1) * alpha^(2*N-3) / (N^2)))}
  
  if (DPPfamily=="Bessel"){return (-2*rho*pi*alpha/(1-rho*pi*alpha^2)^2)}
}