# Intent: Simulate GLMM data
simulateGLMM<-function(ngroup=2, nblock=4, reps.per.group = 20, alpha=7, sd.beta=2, sd=.01){

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
  data<- data.frame(Sample.Num=rep(1:reps.per.group,times=(nsample/nblock)), Group=rep(group, (nsample/ngroup)) , Block = rep(block, each=(nsample/nblock)), Eta=NA, Y=NA)

  # Ok and simplest for now to simulate data here but also will be possible to do it through nimble later. -Perry
  # generate random effect of block
  r.block <-rnorm(n=nblock,mean=0,sd=sd)
  
  data$Eta<- exp(rep(alpha,nsample)+beta[data$Group]+r.block[data$Block])
  
  data$Y<-rpois(nsample, lambda=data$Eta) 
  
  return(list(data=data, nsample=nsample, alpha=alpha, beta=beta, r.block=r.block)) 
    }

sim.glm<-simulateGLMM()
head(sim.glm$data)
table(sim.glm$data$Group, sim.glm$data$Block)

#library(lme4)
sim.glm$data$Group <-factor(sim.glm$data$Group)
glmm.fit <-glmer(Y~Group+(1|Block),family="poisson", data=sim.glm$data)
summary(glmm.fit)

### Attempt to Simulate with Nimble ###
library("Nimble")
