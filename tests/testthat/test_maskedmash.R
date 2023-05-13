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
datax = simDataI.ult(2000,Ulist,prior = 'uniform')

fit_maskedmash = maskedmash_wrapper(datax$Bhat)
fit_mm_fdr = maskedmashFDR(fit_maskedmash)
fdp(fit_mm_fdr$rej.set,which(datax$B!=0))
powr(fit_mm_fdr$rej.set,which(datax$B!=0))

fit_maskedmash_onZ = mash_on_maskedZ(datax$Bhat,verbose = T)
fit_maskedmash_onZ_fdr = maskedmashFDR(fit_maskedmash_onZ)
fdp(fit_maskedmash_onZ_fdr$rej.set,which(datax$B!=0))
powr(fit_maskedmash_onZ_fdr$rej.set,which(datax$B!=0))


fit_ash_on_maskedZ = ash_on_maskedZ(datax$Bhat,verbose = T)
fit_ash_on_maskedZ_fdr = maskedmashFDR(fit_ash_on_maskedZ)
fdp(fit_ash_on_maskedZ_fdr$rej.set,which(datax$B!=0))
powr(fit_ash_on_maskedZ_fdr$rej.set,which(datax$B!=0))

fit_mash = mash_wrapper(datax$Bhat)
fit_mash_fdr = mashFDR(fit_mash)
fdp(fit_mash_fdr,which(datax$B!=0))
powr(fit_mash_fdr,which(datax$B!=0))











library(parallel)
sim_study = function(Ulist,prior_weight,N,prior,df,
                     seed=12345,nreps = 20,mc.cores = 4,npc=5,mean.range=4){

  set.seed(seed)

  result = mclapply(1:nreps,function(i){
    print(paste('Running ',i,sep=''))
    # generate data
    simdata = simDataI.ult(N,Ulist,prior_weight,prior,df,mean.range)
    #datax = mash_set_data(simdata$Bhat,simdata$Shat)
    # fit model
    mash.out = mash_wrapper(simdata$Bhat,npc=npc,adjust = 'lb',md.method = 'teem')
    mash.out2 = mash_wrapper(simdata$Bhat,npc=npc,adjust = 'const',md.method = 'teem')
    # mashonmaskedZ.out = mash_mask(simdata$Bhat,simdata$P,p.thresh = 0.5,npc=npc)
    # mashonmaskedZ.out2 = mash_mask(simdata$Bhat,simdata$P,p.thresh = 0.4,npc=npc)
    # mashonmaskedZ.out3 = mash_mask(simdata$Bhat,simdata$P,p.thresh = 0.3,npc=npc)
    # mashonmaskedZ.out4 = mash_mask(simdata$Bhat,simdata$P,p.thresh = 0.2,npc=npc)
    # mashonmaskedZ.out5 = mash_mask(simdata$Bhat,simdata$P,p.thresh = 0.1,npc=npc)
    # mashonmaskedZ.out6 = mash_mask(simdata$Bhat,simdata$P,p.thresh = 0.05,npc=npc)
    #
    # maskedmash.out = maskedmash_wrapper(simdata$Bhat,simdata$P,p.thresh = 0.5,npc=npc)
    # maskedmash.out2 = maskedmash_wrapper(simdata$Bhat,simdata$P,p.thresh=0.4,npc=npc)
    # maskedmash.out3 = maskedmash_wrapper(simdata$Bhat,simdata$P,p.thresh=0.3,npc=npc)
    # maskedmash.out4 = maskedmash_wrapper(simdata$Bhat,simdata$P,p.thresh=0.2,npc=npc)
    maskedmash.out5 = maskedmash_wrapper(simdata$Bhat,simdata$P,p.thresh=0.1,npc=npc,adjust = 'lb')
    maskedmash.out6 = maskedmash_wrapper(simdata$Bhat,simdata$P,p.thresh=0.1,npc=npc,adjust = 'const')
    # maskedmash.out6 = maskedmash_wrapper(simdata$Bhat,simdata$P,p.thresh=0.05,npc=npc)
    list(data=simdata,
         mash.out = mash.out,
         mash.out2 = mash.out2,
         maskedmash.out5 = maskedmash.out5,
         maskedmash.out6 = maskedmash.out6
         )

    # list(data = simdata,
    #      mash.out=mash.out,
    #      mashonmaskedZ.out=mashonmaskedZ.out,
    #      mashonmaskedZ.out2=mashonmaskedZ.out2,
    #      mashonmaskedZ.out3=mashonmaskedZ.out3,
    #      mashonmaskedZ.out4=mashonmaskedZ.out4,
    #      mashonmaskedZ.out5=mashonmaskedZ.out5,
    #      mashonmaskedZ.out6=mashonmaskedZ.out6,
    #      maskedmash.out=maskedmash.out,
    #      maskedmash.out2=maskedmash.out2,
    #      maskedmash.out3=maskedmash.out3,
    #      maskedmash.out4=maskedmash.out4,
    #      maskedmash.out5=maskedmash.out5,
    #      maskedmash.out6=maskedmash.out6)

  },mc.cores = mc.cores)

  result

}


res = sim_study(Ulist,prior_weight =  rep(1/length(Ulist),length(Ulist)),N=length(Ulist)*500,prior='uniform',df=10,
                seed=12345,nreps = 20,mc.cores = 2,npc=5,mean.range=4)
