#####Gauss#####
rho0 = 100
alpha0 = 0.03
S_length = 2
Simu = readRDS(paste("datasets/2000_Gauss_", toString(rho0), "_", toString(alpha0*100), "e-2_[0_", toString(S_length), "].dat", sep=""))
alpha_est = MLEDPP(Simu[[2]], "Gauss", edgecorr=TRUE)$fixedpar$alpha

I = Fisher_Info(Simu[[2]], "Gauss", alpha_est, edgecorr=TRUE)
range = 1.96 * sqrt(1/(I[2,2]-((I[1,2]^2)/I[1,1])))
