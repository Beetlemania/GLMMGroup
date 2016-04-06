# Intent: Simulate GLMM data
simulateGLMM<-function(ngroup=2, nblock=4, reps.per.group = 100, alpha=7, sd.beta=2, sd=.01){

    # 100 is a lot! -Perry
   # generate numbers of samples
  nsample<-ngroup*nblock*reps.per.group
  
  # generate values of group covariate
  group <-1:ngroup
  
   # generate blocking
  block<-1:nblock
  
  # generate betas
  beta<-rnorm(n=ngroup,mean=0,sd=sd.beta)

  # I think we want the groups and the blocks crossed. -Perry
  data<- data.frame(Group=rep(group, (nsample/ngroup)) , Block= rep(block, each=(nsample/nblock)), Y=NA)

  # Ok and simplest for now to simulate data here but also will be possible to do it through nimble later. -Perry
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
