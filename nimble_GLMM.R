library(nimble)

glmmCode <- nimbleCode({
    for(i in 1:ngroup) {
        group.intercept[i] ~ dnorm(0, sd = 1000)
    }

    for(i in 1:nblock) {
        r.block[i] ~ dnorm(0, sd = block.sd)
    }
    block.sd ~ dunif(0, 10)
    
    for(i in 1:nsample) {
        eta[i] <- group.intercept[ group[i] ] + inprod(X[i,1:nslope], slopes[1:nslope]) + r.block[ block[i] ]
        y[i] ~ dpois(exp(eta[i]))
    }
})

## I'm going to rework simulateGLMM to provide information (but not simulated data) for use in the nimble model
##
## Let's add explanatory variables to this
## Instead of alpha and sd.beta, lets call them group.mean and group.sd
## Instead of sd, let's call it block.sd
## Let's call slopes a set of coefficients for X values that we'll generate as uniformly distributed between -1 and 1
## Let's call betas the group.effects
setupGLMMparams <- function(ngroup=2, nblock=4, reps.per.group = 20, group.mean=7, group.sd=2, slopes = c(.2, .4), block.sd=.01){

    ##
    nslope <- length(slopes)
    
   ## generate numbers of samples
    nsample<-ngroup*nblock*reps.per.group
    
    ## generate values of group covariate
    group <-1:ngroup
  
   ## generate blocking
    block<-1:nblock
  
    ## generate group.intercepts (formerly called betas)
    group.intercept <- rnorm(n=ngroup,mean=group.mean,sd=group.sd)

    Xmatrix <- matrix(runif(nslope * nsample, -1, 1), nrow = nsample)

    ## I'm just going to return information needed for the BUGS code, not simulated data itself.
    setupInfo <- list(data = list(X = Xmatrix, y = rep(0, nsample)),
                      constants = list(nblock = nblock, ngroup = ngroup, nsample = nsample, nslope = nslope, group = rep(group, nsample/ngroup), block = rep(block, each = nsample/nblock)),
                      inits = list(slopes = slopes, group.intercept = group.intercept, block.sd = block.sd))

    return(setupInfo)


    
  ##   ## I think we want the groups and the blocks crossed. -Perry
  ##   data <- data.frame(Sample.Num=rep(1:reps.per.group,times=(nsample/nblock)), Group=rep(group, (nsample/ngroup)) , Block = rep(block, each=(nsample/nblock)), Eta=NA, Y=NA)

  ## # Ok and simplest for now to simulate data here but also will be possible to do it through nimble later. -Perry
  ## # generate random effect of block
  ## r.block <-rnorm(n=nblock,mean=0,sd=sd)
  
  ## data$Eta<- exp(rep(alpha,nsample)+beta[data$Group]+r.block[data$Block])
  
  ## data$Y<-rpois(nsample, lambda=data$Eta) 
  
  ## return(list(data=data, nsample=nsample, alpha=alpha, beta=beta, r.block=r.block)) 
}

myGLMMparams <- setupGLMMparams() ## use all defaults

## set up a mode in nimble
glmmModel <- nimbleModel(glmmCode, constants = myGLMMparams$constants, data = myGLMMparams$data, inits = myGLMMparams$inits)

## As an example, we can see our parameter values entered:
glmmModel$slopes

## This function uses the model to simulate data.
## The simulated data will be in the model.  The function also returns it in case we want it for some other purpose.
simulateGLMM <- function(glmmModel) {
    parameters <- glmmModel$getNodeNames(topOnly = TRUE)
    nodesToSimulate <- glmmModel$getDependencies(parameters, self = FALSE, includeData = TRUE)
    glmmModel$simulate( nodesToSimulate, includeData = TRUE )
    data <- values(glmmModel, glmmModel$getNodeNames( dataOnly = TRUE))
    data
}

## This function compiles and runs the MCMC.
## We could keep it simpler by not compiling it each time, but this is quick'n'dirty to make it run
runMCMC <- function(glmmModel, niter = 10000) {
    glmmMCMCconf <- configureMCMC(glmmModel)
    glmmMCMCconf$addMonitors('r.block')
    glmmMCMC <- buildMCMC(glmmMCMCconf)
    compiled <- compileNimble(glmmModel, glmmMCMC)
    compiled$glmmMCMC$run(niter)
    samples <- as.matrix(compiled$glmmMCMC$mvSamples)
    samples
}

## simulate a new data set
simulateGLMM(glmmModel)
## run the MCMC
mcmcSample <- runMCMC(glmmModel, niter = 100000)
