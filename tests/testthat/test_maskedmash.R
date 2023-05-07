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












#
# sim_study = function(Ulist,prior_weight,N,prior,df,
#                      seed=12345,nreps = 20,mc.cores = 4,npc=5,
#                      half.uniform=FALSE,mean.range=4){
#
#   set.seed(seed)
#
#   result = mclapply(1:nreps,function(i){
#     # generate data
#     simdata = simDataI.ult(N,Ulist,prior_weight,prior,df,half.uniform,mean.range)
#     #datax = mash_set_data(simdata$Bhat,simdata$Shat)
#     # fit model
#     mash.out = mash_wrapper(simdata$Bhat,npc=npc)
#
#     mashonmaskedZ.out = mash_mask(simdata$Bhat,simdata$P,p.thresh = 0.5,npc=npc)
#     mashonmaskedZ.out2 = mash_mask(simdata$Bhat,simdata$P,p.thresh = 0.4,npc=npc)
#     mashonmaskedZ.out3 = mash_mask(simdata$Bhat,simdata$P,p.thresh = 0.3,npc=npc)
#     mashonmaskedZ.out4 = mash_mask(simdata$Bhat,simdata$P,p.thresh = 0.2,npc=npc)
#     mashonmaskedZ.out5 = mash_mask(simdata$Bhat,simdata$P,p.thresh = 0.1,npc=npc)
#     mashonmaskedZ.out6 = mash_mask(simdata$Bhat,simdata$P,p.thresh = 0.05,npc=npc)
#
#     maskedmash.out = maskedmash_wrapper(simdata$Bhat,simdata$P,p.thresh = 0.5,npc=npc)
#     maskedmash.out2 = maskedmash_wrapper(simdata$Bhat,simdata$P,p.thresh=0.4,npc=npc)
#     maskedmash.out3 = maskedmash_wrapper(simdata$Bhat,simdata$P,p.thresh=0.3,npc=npc)
#     maskedmash.out4 = maskedmash_wrapper(simdata$Bhat,simdata$P,p.thresh=0.2,npc=npc)
#     maskedmash.out5 = maskedmash_wrapper(simdata$Bhat,simdata$P,p.thresh=0.1,npc=npc)
#     maskedmash.out6 = maskedmash_wrapper(simdata$Bhat,simdata$P,p.thresh=0.05,npc=npc)
#
#     list(data = simdata,
#          mash.out=mash.out,
#          mashonmaskedZ.out=mashonmaskedZ.out,
#          mashonmaskedZ.out2=mashonmaskedZ.out2,
#          mashonmaskedZ.out3=mashonmaskedZ.out3,
#          mashonmaskedZ.out4=mashonmaskedZ.out4,
#          mashonmaskedZ.out5=mashonmaskedZ.out5,
#          mashonmaskedZ.out6=mashonmaskedZ.out6,
#          maskedmash.out=maskedmash.out,
#          maskedmash.out2=maskedmash.out2,
#          maskedmash.out3=maskedmash.out3,
#          maskedmash.out4=maskedmash.out4,
#          maskedmash.out5=maskedmash.out5,
#          maskedmash.out6=maskedmash.out6)
#
#   },mc.cores = mc.cores)
#
#   result
#
# }