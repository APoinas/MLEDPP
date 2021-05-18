#L0
L0=function(rho, alpha, M, Nsum){
  d = dim(M)[1]
  bigM = array(rep(0, d*d), dim=c(d,d))
  for (i in 1:Nsum){
    bigM = bigM + (2*rho*pi*alpha^2)^(i-1)/(i^2*(1+(M/(alpha*i))^2)^(3/2))
  }
  return (rho * bigM)
}

#Derivative of L0 with respect to alpha
DalpL0=function(rho, alpha, M, Nsum){
  d = dim(M)[1]
  bigM = array(rep(0, d*d), dim=c(d,d))
  for (i in 1:Nsum){
    bigM = bigM + (2*rho*pi)^(i-1) * alpha^(2*i-3) * (2*i+1-3/(1 + (M / (i * alpha))^2))    /(i^2*(1+(M/(alpha*i))^2)^(3/2))
  }
  return (rho * bigM)
}

#2nd derivative of L0 with respect to alpha
D2alpL0=function(rho, alpha, M, Nsum){
  d = dim(M)[1]
  bigM = array(rep(0, d*d), dim=c(d,d))
  for (i in 1:Nsum){
    bigM = bigM + ((2*rho*pi)^(i-1)/(i^2)) * ( ( (2*i-2) * alpha^(2*i-4) * (2*i-3/(1 + (M / (i * alpha))^2))    /((1+(M/(alpha*i))^2)^(3/2)) ) + ( 3*(M/i)^2 * alpha^(2*i-6) * (2*i-5/(1 + (M / (i * alpha))^2))    /((1+(M/(alpha*i))^2)^(5/2)) ) )
  }
  return (rho * bigM)
}

#Derivative of L0 with respect to rho
DrhoL0=function(rho, alpha, M, Nsum){
  d = dim(M)[1]
  bigM = array(rep(0, d*d), dim=c(d,d))
  for (i in 1:Nsum){
    bigM = bigM + (2*rho*pi*alpha^2)^(i-1)/(i*(1+(M/(alpha*i))^2)^(3/2))
  }
  return (bigM)
}

#2nd derivative of L0 with respect to rho
D2rhoL0=function(rho, alpha, M, Nsum){
  d = dim(M)[1]
  bigM = array(rep(0, d*d), dim=c(d,d))
  for (i in 1:Nsum){
    bigM = bigM + (2*pi*alpha^2)^(i-1) * rho^(i-2) * (i-1)/(i*(1+(M/(alpha*i))^2)^(3/2))
  }
  return (bigM)
}

#Derivative of L0 with respect to alpha and rho
DalpDrhoL0=function(rho, alpha, M, Nsum){
  d = dim(M)[1]
  bigM = array(rep(0, d*d), dim=c(d,d))
  for (i in 1:Nsum){
    bigM = bigM + (2*rho*pi)^(i-1) * alpha^(2*i-3) * (2*i-2+3/(1 + i^2 * alpha^2 / M^2))    /(i*(1+(M/(alpha*i))^2)^(3/2))
  }
  return (bigM)
}

#Integral part of the log-likelihood
Integ = function(rho, alpha, Nsum){
  N = 1:Nsum
  return(-rho * sum( (2*rho*pi*alpha^2)^(N-1) / (N^3)))
}

#Derivative of the integral part of the log-likelihood with respect to alpha
DalpInteg = function(rho, alpha, Nsum){
  N = 1:Nsum
  return(- 2 * rho * sum( (N-1) * (2*rho*pi)^(N-1) * alpha^(2*N-3) / (N^3)))
}

#2nd derivative of the integral part of the log-likelihood with respect to alpha
D2alpInteg = function(rho, alpha, Nsum){
  N = 1:Nsum
  return(- 2 * rho * sum( (N-1) * (2*N-3) * (2*rho*pi)^(N-1) * alpha^(2*N-4) / (N^3)))
}

#Derivative of the integral part of the log-likelihood with respect to rho
DrhoInteg = function(rho, alpha, Nsum){
  N = 1:Nsum
  return(- sum( (2*rho*pi*alpha^2)^(N-1) / (N^2)))
}

#2nd derivative of the integral part of the log-likelihood with respect to rho
D2rhoInteg = function(rho, alpha, Nsum){
  N = 1:Nsum
  return(- sum( (N-1) * (2*pi*alpha^2)^(N-1) * rho^(N-2) / (N^2)))
}

#Derivative of the integral part of the log-likelihood with respect to rho and alpha
DalpDrhoInteg = function(rho, alpha, Nsum){
  N = 1:Nsum
  return(- 2 * sum( (N-1) * (2*rho*pi)^(N-1) * alpha^(2*N-3) / (N^2)))
}