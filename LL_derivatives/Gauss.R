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