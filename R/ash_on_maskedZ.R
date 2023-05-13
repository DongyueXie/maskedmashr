

#'@title run ash condition by condition on masked Z scores
#'@description The procedure masks the z scores, by taking the bigger one(absolute value) in masked set, then fit ash condition by condition
#'@export

ash_on_maskedZ = function(Z,P=NULL,p.thresh=0.5,Shat=NULL,verbose=FALSE,eps=1e-10){

  t0 = Sys.time()
  if(is.null(Z)){
    stop('Z must be provided')
  }
  if(is.null(P)){
    P = 2*(1-pnorm(abs(Z)))
  }

  R = ncol(Z)
  N = nrow(Z)

  P.temp = P
  P.temp = pmax(P.temp,eps)
  P.temp = pmin(P.temp,1-eps)

  Z = sign(Z)*qnorm(1-P.temp/2)
  if(is.null(p.thresh)){
    p.thresh = 0.5
  }
  thresh = qnorm(1-p.thresh/2)
  Z.mask = t(apply(Z,1,function(x){
    x.alt = mask.func(x)
    is.mask = is.mask.z(x,thresh)
    is.mask*sign(x)*pmax(abs(x),abs(x.alt)) + (1-is.mask)*x
  }))

  lfsr = matrix(nrow = N,ncol=R)
  for(r in 1:R){
    lfsr[,r] = ash(Z.mask[,r],1)$result$lfsr
  }

  out = list(result = list(lfsr=lfsr))
  out$P = P
  out$p.thresh = p.thresh
  t1 = Sys.time()
  out$run_time = difftime(t1,t0)
  out
}
