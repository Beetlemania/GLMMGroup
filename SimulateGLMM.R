# Intent: Simulate GLMM data
simulateGLMM<-function(ngroup=2, nblock=4, alpha=7, sd.beta=2, sd=.01){
  
   # generate numbers of samples
  nsample<-ngroup*nblock*100 
  
  # generate values of group covariate
  group <-1:ngroup
  
   # generate blocking
  block<-1:nblock
  
  # generate betas
  beta<-rnorm(n=ngroup,mean=0,sd=sd.beta)
  
  data<- data.frame(Group=rep(group, (nsample/ngroup)) , Block= rep(block,(nsample/nblock)), Y=NA)
  
  # generate random effect of block
  r.block <-rnorm(n=nblock,mean=0,sd=sd)
  
  eta<-exp(alpha+beta*group+r.block)
  
  data$Y<-rpois(nsample, lambda=eta) 
  
  return(list(data=data, nsample=nsample, alpha=alpha, beta=beta, r.block=r.block)) 
    }

sim.glm<-simulateGLMM()

library(lme4)
sim.glm$data$Group <-factor(sim.glm$data$Group)
glmm.fit <-glmer(Y~Group+(1|Block),family="poisson", data=sim.glm$data)
summary(glmm.fit)
