#########################################################################
# GLMM Group 2016
#########################################################################
#
# BPA with JAGS
#
# Based on the book "Bayesian population analysis using WinBUGS - a hierarchical perspective" by Marc K?ry & Michael Schaub (2012, Academic Press)
#
#########################################################################

# Michael Schaub, 24.1.2012, revised 9.7.2012, revised 10.12.2013

# This is the BPA code "translated" into JAGS. In many cases there are no differences at all compared to the BUGS code, but sometimes  small changes either in the model code or in the function to provide initial values were necessary. When such changes occur, they are most often commented in the code.

#########################################################################

# Remove all previously loaded objects; start with a clean slate
rm(list = ls())

# Load packages
library(lattice)
library(coda)
library(R2WinBUGS)
library(R2jags)

#########################################################################

# Make sure that you have installed JAGS version 3.2.0

#setwd("C:/Users/Adam/Documents/GitHub/GLMMGroup")     # Optional


############################################################################
#
# 4. Introduction to random effects: Conventional Poisson GLMM for count data
#
##############################################################################
library(lme4)

# 4.1. Introduction
# 4.1.1. An example
# Define and plot data
mass <- c(25, 14, 68, 79, 64, 139, 49, 119, 111)
pop <- factor(c(1, 1, 1, 2, 2, 2, 3, 3, 3))
length <- c(1, 14, 22, 2, 9, 20, 2, 13, 22)
plot(length, mass, col = c(rep("red", 3), rep("blue", 3), rep("green", 3)), xlim = c(-1, 25), ylim = c(0, 140), cex = 1.5, lwd = 2, frame.plot = FALSE, las = 1, pch = 16, xlab = "Length", ylab = "Mass")

# Fit fixed-effects model, print regression parameter estimates and plot regression lines
summary(lm <- lm(mass ~ pop-1 + length))
abline(lm$coef[1], lm$coef[4], col = "red", lwd = 3, lty = 2)
abline(lm$coef[2], lm$coef[4], col = "blue", lwd = 3, lty = 2)
abline(lm$coef[3], lm$coef[4], col = "green", lwd = 3, lty = 2)

# Fit mixed model, print random effects and plot regression lines
summary(lmm <- lmer(mass ~ length + (1|pop)))
ranef(lmm)
abline((lmm@fixef[1]+ranef(lmm)$pop)[1,], lmm@fixef[2], col = "red", lwd = 3)
abline((lmm@fixef[1]+ranef(lmm)$pop)[2,], lmm@fixef[2], col = "blue", lwd = 3)
abline((lmm@fixef[1]+ranef(lmm)$pop)[3,], lmm@fixef[2], col = "green", lwd = 3)


# 4.2. Accounting for overdispersion by random effects-modeling in R and WinBUGS
# 4.2.1. Generation and analysis of simulated data
data.fn <- function(n = 40, alpha = 3.5576, beta1 = -0.0912, beta2 = 0.0091, beta3 = -0.00014, sd = 0.1){
   # n: Number of years
   # alpha, beta1, beta2, beta3: coefficients of a 
   #    cubic polynomial of count on year
   # sd: standard deviation of normal distribution assumed for year effects

   # Generate values of time covariate
   year <- 1:n

   # First level of noise: generate random year effects
   eps <- rnorm(n = n, mean = 0, sd = sd)

   # Signal (plus first level of noise): build up systematic part of the GLM and add the random year effects
   log.expected.count <- alpha + beta1 * year + beta2 * year^2 + beta3 * year^3 + eps
   expected.count <- exp(log.expected.count)

   # Second level of noise: generate random part of the GLM: Poisson noise around expected counts
   C <- rpois(n = n, lambda = expected.count)

   # Plot simulated data
   plot(year, C, type = "b", lwd = 2, main = "", las = 1, ylab = "Population size", xlab = "Year", ylim = c(0, 1.1*max(C)))
   lines(year, expected.count, type = "l", lwd = 3, col = "red")

   return(list(n = n, alpha = alpha, beta1 = beta1, beta2 = beta2, beta3 = beta3, year = year, sd = sd, expected.count = expected.count, C = C))
   }

data <- data.fn()

library(lme4)
yr <- factor(data$year)         # Create a factor year
glmm.fit <- lmer(C ~ (1 | yr) + year + I(year^2) + I(year^3), family = poisson, data = data)

mny <- mean(data$year)
sdy <- sd(data$year)
cov1 <- (data$year - mny) / sdy
cov2 <- cov1 * cov1
cov3 <- cov1 * cov1 * cov1
glmm.fit <- lmer(C ~ (1 | yr) + cov1 + cov2 + cov3, family = poisson, data = data)
glmm.fit

R.predictions <- exp(fixef(glmm.fit)[1] + fixef(glmm.fit)[2]*cov1 + fixef(glmm.fit)[3]*cov2 + fixef(glmm.fit)[4]*cov3 + unlist(ranef(glmm.fit)))
lines(data$year, R.predictions, col = "green", lwd = 2, type = "l")

# Specify model in BUGS language
sink("GLMM_Poisson.jags")
cat("
model {

# Priors
alpha ~ dunif(-20, 20)
beta1 ~ dunif(-10, 10)
beta2 ~ dunif(-10, 10)
beta3 ~ dunif(-10, 10)
tau <- 1 / (sd*sd)
sd ~ dunif(0, 5)

# Likelihood: note key components of a GLM in one line each
for (i in 1:n){
   C[i] ~ dpois(lambda[i])          # 1. Distribution for random part
   log(lambda[i]) <- log.lambda[i]  # 2. Link function
   log.lambda[i] <- alpha + beta1 * year[i] + beta2 * pow(year[i],2) + beta3 * pow(year[i],3) + eps[i]       # 3. Linear predictor incl. random year effect
   eps[i] ~ dnorm(0, tau)     # 4. Definition of random effects dist
   }
}
",fill = TRUE)
sink()

# Bundle data
win.data <- list(C = data$C, n = length(data$C), year = cov1)

# Initial values
inits <- function() list(alpha = runif(1, -2, 2), beta1 = runif(1, -3, 3), sd = runif(1, 0,1))

# Parameters monitored
params <- c("alpha", "beta1", "beta2", "beta3", "lambda", "sd", "eps")

# MCMC settings
ni <- 30000
nt <- 10
nb <- 20000
nc <- 3

# Call JAGS from R (BRT <1 min)
out <- jags(win.data, inits, params, "GLMM_Poisson.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

# Summarize posteriors
print(out, dig = 2)

JAGS.predictions <- out$BUGSoutput$mean$lambda
lines(data$year, JAGS.predictions, col = "blue", lwd = 2, type = "l", lty = 2)

glm.fit <- glm(C ~ cov1 + cov2 + cov3, family = poisson, data = data)
summary(glm.fit)
summary(glmm.fit)


# 4.2.2. Analysis of real data
# Read data again
peregrine <- read.table("falcons.txt", header = TRUE)

yr <- factor(peregrine$Year)
mny <- mean(peregrine$Year)
sdy <- sd(peregrine$Year)
cov1 <- (peregrine$Year - mny) / sdy
cov2 <- cov1 * cov1
cov3 <- cov1 * cov1 * cov1
glmm <- lmer(peregrine$Pairs ~ (1 | yr) + cov1 + cov2 + cov3, family = poisson, data = peregrine)
glmm

# Bundle data
win.data <- list(C = peregrine$Pairs, n = length(peregrine$Pairs), year = cov1)

# MCMC settings (may have to adapt)
ni <- 30000
nt <- 10
nb <- 20000
nc <- 3

# Call JAGS from R (BRT < 1 min)
out <- jags(win.data, inits, params, "GLMM_Poisson.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

# Summarize posteriors
print(out, dig = 3)


# 4.3. Mixed models with random effects for variability among groups (site and year effects)
# 4.3.1. Generation and analysis of simulated data
data.fn <- function(nsite = 5, nyear = 40, alpha = 4.18456, beta1 = 1.90672, beta2 = 0.10852, beta3 = -1.17121, sd.site = 0.5, sd.year = 0.2){
   # nsite: Number of populations
   # nyear: Number of years
   # alpha, beta1, beta2, beta3: cubic polynomial coefficients of year
   # sd.site: standard deviation of the normal distribution assumed for the population intercepts alpha
   # sd.year: standard deviation of the normal distribution assumed for the year effects
   # We standardize the year covariate so that it runs from about 1 to 1

   # Generate data structure to hold counts and log(lambda)
   C <- log.expected.count <- array(NA, dim = c(nyear, nsite))

   # Generate covariate values
   year <- 1:nyear
   yr <- (year-20)/20	# Standardize
   site <- 1:nsite

   # Draw two sets of random effects from their respective distribution
   alpha.site <- rnorm(n = nsite, mean = alpha, sd = sd.site)
   eps.year <- rnorm(n = nyear, mean = 0, sd = sd.year)

   # Loop over populations
   for (j in 1:nsite){
      # Signal (plus first level of noise): build up systematic part of the GLM including random site and year effects
      log.expected.count[,j] <- alpha.site[j] + beta1 * yr + beta2 * yr^2 + beta3 * yr^3 + eps.year
      expected.count <- exp(log.expected.count[,j])

      # Second level of noise: generate random part of the GLM: Poisson noise around expected counts
      C[,j] <- rpois(n = nyear, lambda = expected.count)
      }

   # Plot simulated data
   matplot(year, C, type = "l", lty = 1, lwd = 2, main = "", las = 1, ylab = "Population size", xlab = "Year")

   return(list(nsite = nsite, nyear = nyear, alpha.site = alpha.site, beta1 = beta1, beta2 = beta2, beta3 = beta3, year = year, sd.site = sd.site, sd.year = sd.year, expected.count = expected.count, C = C))
   }

data <- data.fn(nsite = 100, nyear = 40, sd.site = 0.3, sd.year = 0.2)

# Specify model in BUGS language
sink("GLMM_Poisson.jags")
cat("
model {

# Priors
for (j in 1:nsite){
   alpha[j] ~ dnorm(mu, tau.alpha)		# 4. Random site effects
   }
mu ~ dnorm(0, 0.01)				# Hyperparameter 1
tau.alpha <- 1 / (sd.alpha*sd.alpha)	        # Hyperparameter 2
sd.alpha ~ dunif(0, 2)
for (p in 1:3){
   beta[p] ~ dnorm(0, 0.01)
   }

tau.year <- 1 / (sd.year*sd.year)
sd.year ~ dunif(0, 1)				# Hyperparameter 3

# Likelihood
for (i in 1:nyear){
   eps[i] ~ dnorm(0, tau.year)                # 4. Random year effects
   for (j in 1:nsite){
      C[i,j] ~ dpois(lambda[i,j])             # 1. Distribution for random part
      lambda[i,j] <- exp(log.lambda[i,j])     # 2. Link function
      log.lambda[i,j] <- alpha[j] + beta[1] * year[i] + beta[2] * pow(year[i],2) + beta[3] * pow(year[i],3) + eps[i]    # 3. Linear predictor including random site and random year effects
      }  #j
   }  #i
}
",fill = TRUE)
sink()


# Bundle data
win.data <- list(C = data$C, nsite = ncol(data$C), nyear = nrow(data$C), year = (data$year-20) / 20) # Note year standardized

# Initial values
inits <- function() list(mu = runif(1, 0, 2), alpha = runif(data$nsite, -1, 1), beta = runif(3, -1, 1), sd.alpha = runif(1, 0, 0.1), sd.year = runif(1, 0, 0.1))

# Parameters monitored (may want to add "lambda")
params <- c("mu", "alpha", "beta", "sd.alpha", "sd.year")

# MCMC settings (may have to adapt)
ni <- 100000
nt <- 50
nb <- 50000
nc <- 3

# Call JAGS from R (BRT 98 min)
out <- jags(win.data, inits, params, "GLMM_Poisson.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

# Summarize posteriors
print(out, dig = 3)


# 4.3.2. Analysis of real data set
# Read in the tit data and have a look at them
tits <- read.table("tits.txt", header = TRUE)
str(tits)
C <- as.matrix(tits[5:13])
obs <- as.matrix(tits[14:22])
first <- as.matrix(tits[23:31])

matplot(1999:2007, t(C), type = "l", lty = 1, lwd = 2, main = "", las = 1, ylab = "Territory counts", xlab = "Year", ylim = c(0, 80), frame = FALSE)

table(obs)
length(table(obs))

apply(first, 2, sum, na.rm = TRUE)

a <- as.numeric(levels(factor(obs)))     # All the levels, numeric
newobs <- obs                            # Gets ObsID from 1:271
for (j in 1:length(a)){newobs[which(obs==a[j])] <- j }
table(newobs)

newobs[is.na(newobs)] <- 272
table(newobs)
first[is.na(first)] <- 0
table(first)

#  (a) Null or intercept-only model
# Specify model in BUGS language
sink("GLM0.jags")
cat("
model {

# Prior
alpha ~ dnorm(0, 0.01)    # log(mean count)

# Likelihood
for (i in 1:nyear){
   for (j in 1:nsite){
      C[i,j] ~ dpois(lambda[i,j])
      lambda[i,j] <- exp(log.lambda[i,j])
      log.lambda[i,j] <- alpha
      }  #j
   }  #i
}
",fill = TRUE)
sink()

# Bundle data
win.data <- list(C = t(C), nsite = nrow(C), nyear = ncol(C))

# Initial values
inits <- function() list(alpha = runif(1, -10, 10))

# Parameters monitored
params <- c("alpha")

# MCMC settings
ni <- 1200
nt <- 2
nb <- 200
nc <- 3

# Call JAGS from R (BRT < 1 min)
out0 <- jags(win.data, inits, params, "GLM0.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

# Summarize posteriors
print(out0, dig = 3)


#  (b) Fixed site effects
# Specify model in BUGS language
sink("GLM1.jags")
cat("
model {

# Priors
for (j in 1:nsite){
   alpha[j] ~ dnorm(0, 0.01)     # Site effects
   }

# Likelihood
for (i in 1:nyear){
   for (j in 1:nsite){
      C[i,j] ~ dpois(lambda[i,j])
      lambda[i,j] <- exp(log.lambda[i,j])
      log.lambda[i,j] <- alpha[j]
      }  #j
   }  #i
} 
",fill = TRUE)
sink()

# Bundle data
win.data <- list(C = t(C), nsite = nrow(C), nyear = ncol(C))

# Initial values (not required for all)
inits <- function() list(alpha = runif(235, -1, 1))

# Parameters monitored
params <- c("alpha")

# MCMC settings
ni <- 1200
nt <- 2
nb <- 200
nc <- 3

# Call JAGS from R (BRT < 1 min)
out1 <- jags(win.data, inits, params, "GLM1.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

# Summarize posteriors
print(out1, dig = 2)


#  (c) Fixed site and fixed year effects
# Specify model in BUGS language
sink("GLM2.jags")
cat("
model {

# Priors
for (j in 1:nsite){           # site effects
   alpha[j] ~ dnorm(0, 0.01)
   }

for (i in 2:nyear){           # nyear-1 year effects
   eps[i] ~ dnorm(0, 0.01)
   }
eps[1] <- 0                   # Aliased

# Likelihood
for (i in 1:nyear){
   for (j in 1:nsite){
      C[i,j] ~ dpois(lambda[i,j])
      lambda[i,j] <- exp(log.lambda[i,j])
      log.lambda[i,j] <- alpha[j] + eps[i]
      }  #j
   }  #i
}
",fill = TRUE)
sink()

# Bundle data
win.data <- list(C = t(C), nsite = nrow(C), nyear = ncol(C))

# Initial values
inits <- function() list(alpha = runif(235, -1, 1), eps = c(NA, runif(8, -1, 1)))

# Parameters monitored
params <- c("alpha", "eps")

# MCMC settings
ni <- 1200
nt <- 2
nb <- 200
nc <- 3

# Call JAGS from R (BRT < 1 min)
out2 <- jags(win.data, inits, params, "GLM2.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

# Summarize posteriors
print(out2, dig = 2)


#  (d) Random site effects (no year effects)
# Specify model in BUGS language
sink("GLMM1.jags")
cat("
model {

# Priors
for (j in 1:nsite){
   alpha[j] ~ dnorm(mu.alpha, tau.alpha)   # Random site effects
   }
mu.alpha ~ dnorm(0, 0.01)
tau.alpha <- 1/ (sd.alpha * sd.alpha)
sd.alpha ~ dunif(0, 5)

# Likelihood
for (i in 1:nyear){
   for (j in 1:nsite){
      C[i,j] ~ dpois(lambda[i,j])
      lambda[i,j] <- exp(log.lambda[i,j])
      log.lambda[i,j] <- alpha[j]
      }  #j
   }  #i
}
",fill = TRUE)
sink()

# Bundle data
win.data <- list(C = t(C), nsite = nrow(C), nyear = ncol(C))

# Initial values
inits <- function() list(mu.alpha = runif(1, 2, 3))

# Parameters monitored
params <- c("alpha", "mu.alpha", "sd.alpha")

# MCMC settings
ni <- 1200
nt <- 2
nb <- 200
nc <- 3

# Call JAGS from R (BRT < 1 min)
out3 <- jags(win.data, inits, params, "GLMM1.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

# Summarize posteriors
print(out3, dig = 2)


#  (e) Random site and random year effects
# Specify model in BUGS language
sink("GLMM2.jags")
cat("
model {

# Priors
mu ~ dnorm(0, 0.01)                  # Grand mean

for (j in 1:nsite){
   alpha[j] ~ dnorm(0, tau.alpha)    # Random site effects
   }
tau.alpha <- 1/ (sd.alpha * sd.alpha)
sd.alpha ~ dunif(0, 5)

for (i in 1:nyear){
   eps[i] ~ dnorm(0, tau.eps)        # Random year effects
   }
tau.eps <- 1/ (sd.eps * sd.eps)
sd.eps ~ dunif(0, 3)

# Likelihood
for (i in 1:nyear){
   for (j in 1:nsite){
       C[i,j] ~ dpois(lambda[i,j])
      lambda[i,j] <- exp(log.lambda[i,j])
      log.lambda[i,j] <- mu + alpha[j] + eps[i]
      }  #j
   }  #i
}
",fill = TRUE)
sink()

# Bundle data
win.data <- list(C = t(C), nsite = nrow(C), nyear = ncol(C))

# Initial values (not required for all)
inits <- function() list(mu = runif(1, 0, 4), alpha = runif(235, -2, 2), eps = runif(9, -1, 1))

# Parameters monitored
params <- c("mu", "alpha", "eps", "sd.alpha", "sd.eps")

# MCMC settings
ni <- 6000
nt <- 5
nb <- 1000
nc <- 3

# Call JAGS from R (BRT 3 min)
out4 <- jags(win.data, inits, params, "GLMM2.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

# Summarize posteriors
print(out4, dig = 2)


# (f) Random site and random year effects and first-year fixed observer effect
# Specify model in BUGS language
sink("GLMM3.jags")
cat("
model {

# Priors
mu ~ dnorm(0, 0.01)                 # Overall mean
beta2 ~ dnorm(0, 0.01)              # First-year observer effect

for (j in 1:nsite){
   alpha[j] ~ dnorm(0, tau.alpha)   # Random site effects
   }
tau.alpha <- 1/ (sd.alpha * sd.alpha)
sd.alpha ~ dunif(0, 5)

for (i in 1:nyear){
   eps[i] ~ dnorm(0, tau.eps)      # Random year effects
   }
tau.eps <- 1/ (sd.eps * sd.eps)
sd.eps ~ dunif(0, 5)

# Likelihood
for (i in 1:nyear){
   for (j in 1:nsite){
      C[i,j] ~ dpois(lambda[i,j])
      lambda[i,j] <- exp(log.lambda[i,j])
      log.lambda[i,j] <- mu + beta2 * first[i,j] + alpha[j] + eps[i]
      }  #j
   }  #i
}
",fill = TRUE)
sink()

# Bundle data
win.data <- list(C = t(C), nsite = nrow(C), nyear = ncol(C), first = t(first))

# Initial values
inits <- function() list(mu = runif(1, 0, 4), beta2 = runif(1, -1, 1), alpha = runif(235, -2, 2), eps = runif(9, -1, 1))

# Parameters monitored
params <- c("mu", "beta2", "alpha", "eps", "sd.alpha", "sd.eps")

# MCMC settings
ni <- 6000
nt <- 5
nb <- 1000
nc <- 3

# Call JAGS from R (BRT 3 min)
out5 <- jags(win.data, inits, params, "GLMM3.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

# Summarize posteriors
print(out5, dig = 2)


#  (g) Random site and random year effects, first-year fixed observer effect and overall linear time trend
# Specify model in BUGS language
sink("GLMM4.jags")
cat("
model {

# Priors
mu ~ dnorm(0, 0.01)                  # Overall intercept
beta1 ~ dnorm(0, 0.01)               # Overall trend 
beta2 ~ dnorm(0, 0.01)               # First-year observer effect

for (j in 1:nsite){
   alpha[j] ~ dnorm(0, tau.alpha)    # Random site effects
   }
tau.alpha <- 1/ (sd.alpha * sd.alpha)
sd.alpha ~ dunif(0, 5)

for (i in 1:nyear){
   eps[i] ~ dnorm(0, tau.eps)        # Random year effects
   }
tau.eps <- 1/ (sd.eps * sd.eps)
sd.eps ~ dunif(0, 3)

# Likelihood
for (i in 1:nyear){
   for (j in 1:nsite){
      C[i,j] ~ dpois(lambda[i,j])
      lambda[i,j] <- exp(log.lambda[i,j])
      log.lambda[i,j] <- mu + beta1 * year[i] + beta2 * first[i,j] + alpha[j] + eps[i]
      }  #j
   }  #i
}
",fill = TRUE)
sink()

# Bundle data
win.data <- list(C = t(C), nsite = nrow(C), nyear = ncol(C), first = t(first), year = ((1:9)-5) / 4)

# Initial values
inits <- function() list(mu = runif(1, 0, 4), beta1 = runif(1, -1, 1), beta2 = runif(1, -1, 1), alpha = runif(235, -2, 2), eps = runif(9, -1, 1))

# Parameters monitored
params <- c("mu", "beta1", "beta2", "alpha", "eps", "sd.alpha", "sd.eps")

# MCMC settings
ni <- 12000
nt <- 6
nb <- 6000
nc <- 3

# Call JAGS from R (BRT 7 min)
out6 <- jags(win.data, inits, params, "GLMM4.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

# Summarize posteriors
print(out6, dig = 2)


# (h) The full model 
# Specify model in BUGS language
sink("GLMM5.jags")
cat("
model {

# Priors
mu ~ dnorm(0, 0.01)                  # Overall intercept
beta1 ~ dnorm(0, 0.01)               # Overall trend 
beta2 ~ dnorm(0, 0.01)               # First-year observer effect

for (j in 1:nsite){
   alpha[j] ~ dnorm(0, tau.alpha)    # Random site effects
   }
tau.alpha <- 1/ (sd.alpha * sd.alpha)
sd.alpha ~ dunif(0, 3)

for (i in 1:nyear){
   eps[i] ~ dnorm(0, tau.eps)        # Random year effects
   }
tau.eps <- 1/ (sd.eps * sd.eps)
sd.eps ~ dunif(0, 1)

for (k in 1:nobs){
   gamma[k] ~ dnorm(0, tau.gamma)   # Random observer effects
   }
tau.gamma <- 1/ (sd.gamma * sd.gamma)
sd.gamma ~ dunif(0, 1)

# Likelihood
for (i in 1:nyear){
   for (j in 1:nsite){
      C[i,j] ~ dpois(lambda[i,j])
      lambda[i,j] <- exp(log.lambda[i,j])
      log.lambda[i,j] <- mu + beta1 * year[i] + beta2 * first[i,j] + alpha[j] + gamma[newobs[i,j]] + eps[i]
      }  #j
   }  #i
}
",fill = TRUE)
sink()

# Bundle data
win.data <- list(C = t(C), nsite = nrow(C), nyear = ncol(C), nobs = 272, newobs = t(newobs), first = t(first), year = ((1:9)-5) / 4)

# Initial values
inits <- function() list(mu = runif(1, 0, 4), beta1 = runif(1, -1, 1), beta2 = runif(1, -1, 1), alpha = runif(235, -1, 1), gamma = runif(272, -1, 1), eps = runif(9, -1, 1))

# Parameters monitored
params <- c("mu", "beta1", "beta2", "alpha", "gamma", "eps", "sd.alpha", "sd.gamma", "sd.eps")

# MCMC settings
ni <- 12000
nt <- 6
nb <- 6000
nc <- 3

# Call JAGS from R (BRT 11 min)
out7 <- jags(win.data, inits, params, "GLMM5.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

# Summarize posteriors
print(out7, dig = 2)