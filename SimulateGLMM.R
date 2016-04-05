# Intent: Simulate GLMM data
simulateGLMM<-function(ngroup=2, nblock=4, nsample=140, alpha=7, beta1=.05, sd=.01){
  
  # generate values of group covariate
  group <-1:ngroup
  
  # generate random effect of block
  r.block <-rnorm(n=nblock,mean=0,sd=sd)
  
  eta<-exp(alpha+beta1*group+r.block)
  
  Y<-rpois(nsample, lambda=eta) 
  
  return(list(ngroup=ngroup, nblock=nblock, alpha=alpha, beta1=beta1, sd=sd, eta=eta, Y=Y))
    }

simulateGLMM()
