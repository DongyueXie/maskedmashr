R = 5
u = toeplitz(c(1,0.4,0.3,0.2,0.1))
u[4:5,] = 0
u[,4:5] = 0

u2 = toeplitz(c(1,-0.8,0,0,0))
u2[1:3,] = 0
u2[,1:3] = 0

Ulist = list(U1 = matrix(0,nrow=R,ncol=R),
             U2 = diag(R),
             U3 = u,
             U4 = u2)
datax = simDataI.ult(1000,Ulist,prior = 'uniform')

fit_maskedmash = maskedmash_wrapper(datax$Bhat)
fit_mm_fdr = maskedmashFDR(fit_maskedmash)

fit_maskedmash_onZ = mash_on_maskedZ(datax$Bhat)
fit_maskedmash_onZ_fdr = maskedmashFDR(fit_maskedmash_onZ)

fit_mash = mash_wrapper(datax$Bhat)
fit_mash_fdr = mashFDR(fit_mash)
