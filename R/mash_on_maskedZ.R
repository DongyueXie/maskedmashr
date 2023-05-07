

#'@title mask on masked Z scores
#'@description The procedure masks the z scores, by taking the bigger one(absolute value) in masked set, then fit mash.
#'@export

mash_on_maskedZ = function(Z,P=NULL,p.thresh=0.5,npc=5,Shat=NULL,verbose=FALSE,eps=1e-10){

  if(is.null(Z)){
    stop('Z must be provided')
  }
  if(is.null(P)){
    P = 2*(1-pnorm(abs(Z)))
  }

  R = ncol(Z)
  N = nrow(Z)

  # if(!is.null(P)){
  #   Z = sign(Bhat)*qnorm(1-P/2)
  # }else{
  #   Z = Bhat
  # }

  P.temp = P
  P.temp = pmax(P.temp,eps)
  P.temp = pmin(P.temp,1-eps)

  Z = sign(Z)*qnorm(1-P.temp/2)
  if(is.null(p.thresh)){
    p.thresh = 0.5
  }
  thresh = qnorm(1-p.thresh/2)
  # mask z scores
  Z.mask = t(apply(Z,1,function(x){
    x.alt = mask.func(x)
    is.mask = is.mask.z(x,thresh)
    is.mask*sign(x)*pmax(abs(x),abs(x.alt)) + (1-is.mask)*x
  }))

  out = mash_wrapper(Bhat=Z.mask,npc=npc,Shat=Shat,verbose=verbose)
  out$P = P
  out$p.thresh = p.thresh
  out
}
