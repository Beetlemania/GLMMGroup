# Intent: Simulate GLMM data

 eta=.5
simulateGLMMM<-function(eta, ngroup, nblock,nsample, ncovs){
  
  Y<-data.frame(rep(NA,nsample))
  Y<-rpois(nsample, lambda=eta)
    
  for(i in 1:nblock){
  REsigma[i]<-dunif(x=nsample,0,1000)
  }
}
 
 Ni <- rpois(50, lambda = 4)
