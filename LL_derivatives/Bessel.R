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
