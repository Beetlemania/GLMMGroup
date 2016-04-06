#### Bayesian GLMM analysis

#### Preliminaries
rm(list = ls())
my_packages<-c('data.table', 'tidyr', 'lattice', 'dplyr',
               'snow', 'dclone', 'rjags', 'R2jags')
lapply(my_packages, require, character.only=T)

## Data restructuring functions
source("functions_glmm_analysis.R")

#### Load count data
countDataList <- readRDS("potato_psyllid_count_data_for_GLMM.rds")
countData <- countDataList[[2]]
# Remove NAs from covariates
countData <- dplyr::filter(countData, !is.na(aet) & !is.na(cwd) & !is.na(tmn) & !is.na(tmx))

### There are too many zeros (zero-inflated data)
### In future iterations, could model zero-inflated Poisson distribution
### For now, just subsample zeros (absences)
### Use SDM rule-of-thumb of absences = 10*presences

# Data with detections
detections <- countData[countData$count > 0,]
# Number of desired nondetections
nabsence <- detections %>% nrow()*10
# Subsample nondetections and recombine with detections
countData <- countData[countData$count == 0,] %>% 
  sample_n(., size = nabsence, replace = FALSE) %>% rbind(.,detections)


#### Transform data set for JAGS model
# Create a year-month variable
countData$ymo <- paste(countData$year, countData$month, sep = "-")
# Make a data matrix for the count (response) data
countMatrix <- makeEcoDataMatrix("count")

## Make separate data matrices for the covariates, use standardized covariates
# make standardized year and month
countData$stdyear <- standardize(countData$year)
countData$stdmonth <- standardize(countData$month)
covars <- names(countData)[c(9:13,16:17)]
# Covariate matrices can't have NAs for jags to run, fill missing cells with 0
jagsData <- lapply(covars, function(x) makeEcoDataMatrix(x, fill = 0))
names(jagsData) <- covars

# Combine all matrices into a list
jagsGLMMdata <- list(countMatrix = countMatrix,
                     year = jagsData$stdyear,
                     month = jagsData$stdmonth,
                     list_length = jagsData$stdlnlist_length,
                     aet = jagsData$stdaet,
                     cwd = jagsData$stdcwd,
                     tmn = jagsData$stdtmn,
                     tmx = jagsData$stdtmx,
                     nlist = nrow(countMatrix),
                     nsite = ncol(countMatrix))


################################################################################
#### Binomial GLMM Model specification
sink("climate_glmm_model.jags")
cat("
    model {    
    ## Priors
    
    # For the site random effect
    for(j in 1:nsite) { 
      alpha[j] ~ dnorm(mu.alpha, tau.alpha) 
    }
    mu.alpha ~ dnorm(0, 0.001)
    tau.alpha <- 1 / (sigma.alpha * sigma.alpha)
    sigma.alpha ~ dunif(0, 5)
    
    # Grand mean
    mu ~ dnorm(0, 0.01)

    # For the fixed effect coefficients
    beta1 ~ dnorm(0, 0.001)
    beta2 ~ dnorm(0, 0.001)
    beta3 ~ dnorm(0, 0.001)
    beta4 ~ dnorm(0, 0.001)
    beta5 ~ dnorm(0, 0.001)
    beta6 ~ dnorm(0, 0.001)
    beta7 ~ dnorm(0, 0.001)
    beta8 ~ dnorm(0, 0.001)

    # Likelihood
    for (i in 1:nlist){ # i = events (year-months)
      for(j in 1:nsite) { # j = sites
        detectionMatrix[i,j] ~ dbern(p[i,j])  # Distribution for random part
        logit(p[i,j]) <- mu + beta1*year[i,j] + beta2*month[i,j] + beta3*pow(month[i,j],2) +
                          beta4*list_length[i,j] + 
                          beta5*aet[i,j] + beta6*cwd[i,j] + beta7*tmn[i,j] + beta8*tmx[i,j] +
			  alpha[j]
      } #j
    } #i
    }",fill = TRUE)
sink()

#########################################################################################

#### Specifications of JAGS run
# Initial values
# Specify initial values for mu.alpha, sigma.alpha, and beta1
inits <- function() list(mu.alpha = runif(1, -3, 3),
                         sigma.alpha = runif(1, 0, 5),
                         mu = runif(1, -3, 3),
			 beta1 = runif(1, -3, 3),
                         beta2 = runif(1, -3, 3),
                         beta3 = runif(1, -3, 3),
                         beta4 = runif(1, -3, 3),
                         beta5 = runif(1, -3, 3),
                         beta6 = runif(1, -3, 3),
                         beta7 = runif(1, -3, 3),
                         beta8 = runif(1, -3, 3))
# Monitored parameters
params <- c('mu', 'beta1', 'beta2', 'beta3', 'beta4', 'beta5', 'beta6', 'beta7', 'beta8', 'alpha')
# MCMC specifications
ni=81000; nt=10; nc=3
# for jags.parfit(), burn-in iterations = n.adapt + n.update
n.adapt <- 500; n.update <- 500
nb <- n.adapt + n.update


#### Non-parallel jags()
# glmmOutput <- jags(data = jagsTestdata, 
#                 inits = inits, 
#                 parameters.to.save = params, 
#                 model.file = "climate_glmm_model.jags", 
#                 n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, 
#                 working.directory = getwd())    


#### Parallel JAGS
# Make a SOCK cluster using snow
cl <- makeCluster(3, type = "SOCK")
date()
# Call to jags.parfit
glmmOutput <- jags.parfit(cl, data = jagsGLMMdata,
                            params = params,
                            model = "climate_glmm_model.jags",
                            inits = inits,
                            n.adapt = n.adapt, n.update = n.update,
                            n.iter = ni, thin = nt, n.chains = nc)
date()
stopCluster(cl) # Close the cluster
#### Compute statistics and save output
saveRDS(glmmOutput, file = "climate_glmm_jags_out_full.rds")
glmmdctab <- dctable(glmmOutput)
glmmResults <- data.frame(rbindlist(glmmdctab))
row.names(glmmResults) <- names(glmmdctab)
saveRDS(glmmResults, file = "climate_glmm_jags_out_params.rds")
print("SUCCESS!")

