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


##########################################################################
# 
# 1. Introduction
# 
##########################################################################

# 1.3. The binomial distribution as a canonical description of the observation process
N <- 16                 # Population size of sparrows in the yard
p <- 0.4                # Individual detection probability
rbinom(n = 1, size = N, prob = p)
rbinom(n = 1, size = N, prob = p)
rbinom(n = 1, size = N, prob = p)
C <- rbinom(n = 1000000, size = N, prob = p)
mean(C)
var(C)
sd(C)
hist(C, breaks = 50, col = "gray", main = "", xlab = "Sparrow count", las = 1, freq = FALSE)


# 1.4. Structure and overview of the contents of this book
# Population values for mean and standard deviation of individual lengths
mu <- 65            # Population mean
sigma <- 5          # Population SD

# Draw a single sample of 10 and summarize it
x <- rnorm(n = 10, mean = mu, sd = sigma)
reps <- 10^6
sample.means <- rep(NA, reps)
for (i in 1:reps){
   sample.means[i] <- mean(rnorm(n = 10, mean = mu, sd = sigma))
   }

# Produce figure
par(mfrow = c(1, 2), las = 1)
hist(x, col = "gray", main = "", xlab = "Body length (cm)", las = 1)
abline(v = mu, lwd = 3, col = "red")
abline(v = mean(x), lwd = 3, col = "blue")
hist(sample.means, col = "gray", main = "", xlab = "Body length (cm)", nclass = 50, freq = FALSE, las = 1)
abline(v = mu, lwd = 5, col = "red")
abline(v = mean(sample.means), lwd = 5, col = "blue", lty = 2)
sd(sample.means)




##########################################################################
#
# 2. Very brief introduction to Bayesian statistical modeling
# 
##########################################################################




#########################################################################
# 
# 3. Introduction to the generalized linear model (GLM): The simplest model for count data
#
##########################################################################

# 3.2. Statistical models: Response = Signal + Noise
# 3.2.1. The noise component
plot(density(rbeta(n = 10000, shape1 = 2, shape2 = 4))
hist(rbeta(10^6, 2, 4), nclass = 100, col = "gray")
# 3.2.2. The signal component
# Define and plot data
y <- c(25, 14, 68, 79, 64, 139, 49, 119, 111)
A <- factor(c(1, 1, 1, 2, 2, 2, 3, 3, 3))
X <- c(1, 14, 22, 2, 9, 20, 2, 13, 22)
plot(X, y, col = c(rep("red", 3), rep("blue", 3), rep("green", 3)), xlim = c(-1, 25), ylim = c(0, 140))
summary(fm <- lm(y ~ A-1 + X))
model.matrix(~ A + X)
model.matrix(~ A-1 + X)


# 3.3. Poisson GLM in R and WinBUGS for modeling times series of counts
# 3.3.1. Generation and analysis of simulated data
data.fn <- function(n = 40, alpha = 3.5576, beta1 = -0.0912, beta2 = 0.0091, beta3 = -0.00014){
# n: Number of years
# alpha, beta1, beta2, beta3: coefficients of a 
#    cubic polynomial of count on year

# Generate values of time covariate
year <- 1:n

# Signal: Build up systematic part of the GLM
log.expected.count <- alpha + beta1 * year + beta2 * year^2 + beta3 * year^3
expected.count <- exp(log.expected.count)

# Noise: generate random part of the GLM: Poisson noise around expected counts
C <- rpois(n = n, lambda = expected.count)

# Plot simulated data
plot(year, C, type = "b", lwd = 2, col = "black", main = "", las = 1, ylab = "Population size", xlab = "Year", cex.lab = 1.2, cex.axis = 1.2)
lines(year, expected.count, type = "l", lwd = 3, col = "red")

return(list(n = n, alpha = alpha, beta1 = beta1, beta2 = beta2, beta3 = beta3, year = year, expected.count = expected.count, C = C))
}

data <- data.fn()

fm <- glm(C ~ year + I(year^2) + I(year^3), family = poisson, data = data)
summary(fm)

# Specify model in BUGS language
sink("GLM_Poisson.jags")
cat("
model {

# Priors
alpha ~ dunif(-20, 20)
beta1 ~ dunif(-10, 10)
beta2 ~ dunif(-10, 10)
beta3 ~ dunif(-10, 10)

# Likelihood: Note key components of a GLM on one line each
for (i in 1:n){
   C[i] ~ dpois(lambda[i])          # 1. Distribution for random part
   log(lambda[i]) <- log.lambda[i]  # 2. Link function
   log.lambda[i] <- alpha + beta1 * year[i] + beta2 * pow(year[i],2) + beta3 * pow(year[i],3)                      # 3. Linear predictor
   } #i
}
",fill = TRUE)
sink()

# Bundle data
win.data <- list(C = data$C, n = length(data$C), year = data$year)

# Initial values
inits <- function() list(alpha = runif(1, -2, 2), beta1 = runif(1, -3, 3))

# Parameters monitored
params <- c("alpha", "beta1", "beta2", "beta3", "lambda")

# MCMC settings
ni <- 2000
nt <- 2
nb <- 1000
nc <- 3

# Call JAGS from R
out <- jags(data = win.data, inits = inits, parameters.to.save = params, model.file = "GLM_Poisson.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

# Bundle data
mean.year <- mean(data$year)             # Mean of year covariate
sd.year <- sd(data$year)                 # SD of year covariate
win.data <- list(C = data$C, n = length(data$C), year = (data$year - mean.year) / sd.year)

# Call JAGS from R (BRT < 1 min)
out <- jags(data = win.data, inits = inits, parameters.to.save = params, model.file = "GLM_Poisson.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

print(out, dig = 3)

# New MCMC settings with essentially no burnin
ni <- 100
nt <- 1
nb <- 1

# Call JAGS from R (BRT < 1 min)
tmp <- jags(data = win.data, inits = inits, parameters.to.save = params, model.file = "GLM_Poisson.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

plot(1:40, data$C, type = "b", lwd = 2, col = "black", main = "", las = 1, ylab = "Population size", xlab = "Year")
R.predictions <- predict(glm(C ~ year + I(year^2) + I(year^3), family = poisson, data = data), type = "response")
lines(1:40, R.predictions, type = "l", lwd = 3, col = "green")
JAGS.predictions <- out$BUGSoutput$mean$lambda
lines(1:40, JAGS.predictions, type = "l", lwd = 3, col = "blue", lty = 2)
cbind(R.predictions, JAGS.predictions)

# An option in JAGS to see the traceplots of all parameters is the following:
traceplot(tmp)

# 3.3.2. Analysis of real data set
# Read data
peregrine <- read.table("falcons.txt", header = TRUE)
attach(peregrine)
plot(Year, Pairs, type = "b", lwd = 2, main = "", las = 1, ylab = "Pair count", xlab = "Year", ylim = c(0, 200), pch = 16)

# Bundle data
mean.year <- mean(1:length(Year))        # Mean of year covariate
sd.year <- sd(1:length(Year))            # SD of year covariate
win.data <- list(C = Pairs, n = length(Pairs), year = (1: length(Year) - mean.year) / sd.year)

# Initial values
inits <- function() list(alpha = runif(1, -2, 2), beta1 = runif(1, -3, 3))

# Parameters monitored
params <- c("alpha", "beta1", "beta2", "beta3", "lambda")

# MCMC settings
ni <- 2500
nt <- 2
nb <- 500
nc <- 3

# Call JAGS from R (BRT < 1 min)
out1 <- jags(data = win.data, inits = inits, parameters.to.save = params, model.file = "GLM_Poisson.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

# Summarize posteriors
print(out1, dig = 3) 

JAGS.predictions <- out1$BUGSoutput$mean$lambda
lines(Year, JAGS.predictions, type = "l", lwd = 3, col = "blue", lty = 2)


# 3.4. Poisson GLM for modeling fecundity
plot(Year, Eyasses, type = "b", lwd = 2, main = "", las = 1, ylab = "Nestling count", xlab = "Year", ylim = c(0, 260), pch = 16)

# Bundle data
mean.year <- mean(1:length(Year))   # Mean of year covariate
sd.year <- sd(1:length(Year))       # SD of year covariate
win.data <- list(C = Eyasses, n = length(Eyasses), year = (1: length(Year) - mean.year) / sd.year)

# Call JAGS from R (BRT < 1 min)
out2 <- jags(data = win.data, inits = inits, parameters.to.save = params, model.file = "GLM_Poisson.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

# Plot predictions
JAGS.predictions <- out2$BUGSoutput$mean$lambda
lines(Year, JAGS.predictions, type = "l", lwd = 3, col = "blue")


# 3.5. Binomial GLM for modeling bounded counts or proportions
# 3.5.1. Generation and analysis of simulated data
data.fn <- function(nyears = 40, alpha = 0, beta1 = -0.1, beta2 = -0.9){
# nyears: Number of years
# alpha, beta1, beta2: coefficients

# Generate untransformed and transformed values of time covariate
year <- 1:nyears
YR <- (year-round(nyears/2)) / (nyears / 2)

# Generate values of binomial totals (N)
N <- round(runif(nyears, min = 20, max = 100))

# Signal: build up systematic part of the GLM
exp.p <- plogis(alpha + beta1 * YR + beta2 * (YR^2))

# Noise: generate random part of the GLM: Binomial noise around expected counts (which is N)
C <- rbinom(n = nyears, size = N, prob = exp.p)

# Plot simulated data
plot(year, C/N, type = "b", lwd = 2, col = "black", main = "", las = 1, ylab = "Proportion successful pairs", xlab = "Year", ylim = c(0, 1))
points(year, exp.p, type = "l", lwd = 3, col = "red")

return(list(nyears = nyears, alpha = alpha, beta1 = beta1, beta2 = beta2, year = year, YR = YR, exp.p = exp.p, C = C, N = N))
}

data <- data.fn(nyears = 40, alpha = 1, beta1 = -0.03, beta2 = -0.9)

# Specify model in BUGS language
sink("GLM_Binomial.jags")
cat("
model {

# Priors
alpha ~ dnorm(0, 0.001)
beta1 ~ dnorm(0, 0.001)
beta2 ~ dnorm(0, 0.001)

# Likelihood
for (i in 1:nyears){
   C[i] ~ dbin(p[i], N[i])          # 1. Distribution for random part
   logit(p[i]) <- alpha + beta1 * year[i] + beta2 * pow(year[i],2) # link function and linear predictor
   }
}
",fill = TRUE)
sink()

# Bundle data
win.data <- list(C = data$C, N = data$N, nyears = length(data$C), year = data$YR)

# Initial values
inits <- function() list(alpha = runif(1, -1, 1), beta1 = runif(1, -1, 1), beta2 = runif(1, -1, 1))

# Parameters monitored
params <- c("alpha", "beta1", "beta2", "p")

# MCMC settings
ni <- 2500
nt <- 2
nb <- 500
nc <- 3

# Call JAGS from R (BRT < 1 min)
out <- jags(data = win.data, inits = inits, parameters.to.save = params, model.file = "GLM_Binomial.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

# Plot predictions
JAGS.predictions <- out$BUGSoutput$mean$p
lines(1:length(data$C), JAGS.predictions, type = "l", lwd = 3, col = "blue", lty = 2)


# 3.5.2. Analysis of real data set
# Read data and attach them
peregrine <- read.table("falcons.txt", header = TRUE)
attach(peregrine)

# Bundle data (note yet another standardization for year)
win.data <- list(C = R.Pairs, N = Pairs, nyears = length(Pairs), year = (Year-1985)/ 20)

# Initial values
inits <- function() list(alpha = runif(1, -1, 1), beta1 = runif(1, -1, 1), beta2 = runif(1, -1, 1))

# Parameters monitored
params <- c("alpha", "beta1", "beta2", "p")

# MCMC settings
ni <- 2500
nt <- 2
nb <- 500
nc <- 3

# Call JAGS from R (BRT < 1 min)
out3 <- jags(data = win.data, inits = inits, parameters.to.save = params, model.file = "GLM_Binomial.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

# Summarize posteriors and plot estimates
print(out3, dig = 3)
plot(Year, R.Pairs/Pairs, type = "b", lwd = 2, col = "black", main = "", las = 1, ylab = "Proportion successful pairs", xlab = "Year", ylim = c(0,1))
lines(Year, out3$BUGSoutput$mean$p, type = "l", lwd = 3, col = "blue")




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



########################################################################
#
# 5. State-space models 
#
###########################################################################


# 5.2. A simple model
n.years <- 25           # Number of years
N1 <- 30                # Initial population size
mean.lambda <- 1.02     # Mean annual population growth rate
sigma2.lambda <- 0.02   # Process (temporal) variation of the growth rate
sigma2.y <- 20          # Variance of the observation error

y <- N <- numeric(n.years)
N[1] <- N1
lambda <- rnorm(n.years-1, mean.lambda, sqrt(sigma2.lambda))
for (t in 1:(n.years-1)){
   N[t+1] <- N[t] * lambda[t]
   }

for (t in 1:n.years){
   y[t] <- rnorm(1, N[t], sqrt(sigma2.y))
   }

# Specify model in BUGS language
sink("ssm.jags")
cat("
model { 
# Priors and constraints
N.est[1] ~ dunif(0, 500)            # Prior for initial population size
mean.lambda ~ dunif(0, 10)          # Prior for mean growth rate
sigma.proc ~ dunif(0, 10)           # Prior for sd of state process
sigma2.proc <- pow(sigma.proc, 2)
tau.proc <- pow(sigma.proc, -2)
sigma.obs ~ dunif(0, 100)           # Prior for sd of observation process
sigma2.obs <- pow(sigma.obs, 2)
tau.obs <- pow(sigma.obs, -2)

# Likelihood
# State process
for (t in 1:(T-1)){
   lambda[t] ~ dnorm(mean.lambda, tau.proc) 
   N.est[t+1] <- N.est[t] * lambda[t] 
   }
# Observation process
for (t in 1:T) {
   y[t] ~ dnorm(N.est[t], tau.obs)
   }
}
",fill = TRUE)
sink()

# Bundle data
jags.data <- list(y = y, T = n.years)

# Initial values
inits <- function(){list(sigma.proc = runif(1, 0, 5), mean.lambda = runif(1, 0.1, 2), sigma.obs = runif(1, 0, 10), N.est = c(runif(1, 20, 40), rep(NA, (n.years-1))))} 

# Parameters monitored
parameters <- c("lambda", "mean.lambda", "sigma2.obs", "sigma2.proc", "N.est")

# MCMC settings
ni <- 25000
nt <- 3
nb <- 10000
nc <- 3

# Call JAGS from R (BRT <1 min)
ssm <- jags(jags.data, inits, parameters, "ssm.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())


# Define function to draw a graph to summarize results
graph.ssm <- function(ssm, N, y){
   fitted <- lower <- upper <- numeric()
   n.years <- length(y)
   for (i in 1:n.years){
      fitted[i] <- mean(ssm$BUGSoutput$sims.list$N.est[,i])
      lower[i] <- quantile(ssm$BUGSoutput$sims.list$N.est[,i], 0.025)
      upper[i] <- quantile(ssm$BUGSoutput$sims.list$N.est[,i], 0.975)}
   m1 <- min(c(y, fitted, N, lower))
   m2 <- max(c(y, fitted, N, upper))
   par(mar = c(4.5, 4, 1, 1), cex = 1.2)
   plot(0, 0, ylim = c(m1, m2), xlim = c(0.5, n.years), ylab = "Population size", xlab = "Year", las = 1, col = "black", type = "l", lwd = 2, frame = FALSE, axes = FALSE)
   axis(2, las = 1)
   axis(1, at = seq(0, n.years, 5), labels = seq(0, n.years, 5))
   axis(1, at = 0:n.years, labels = rep("", n.years + 1), tcl = -0.25)
   polygon(x = c(1:n.years, n.years:1), y = c(lower, upper[n.years:1]), col = "gray90", border = "gray90")
   points(N, type = "l", col = "red", lwd = 2)
   points(y, type = "l", col = "black", lwd = 2)
   points(fitted, type = "l", col = "blue", lwd = 2)
   legend(x = 1, y = m2, legend = c("True", "Observed", "Estimated"), lty = c(1, 1, 1), lwd = c(2, 2, 2), col = c("red", "black", "blue"), bty = "n", cex = 1)
}

# Execute function: Produce figure 
graph.ssm(ssm, N, y)


# 5.3. Systematic bias in the observation process
n.years <- 25  # Number of years
N <- rep(50, n.years) 

p <- 0.7
y <- numeric(n.years)
for (t in 1:n.years){
   y[t] <- rbinom(1, N[t], p)
   }
y

# Bundle data
jags.data <- list(y = y, T = n.years)

# Initial values
inits <- function(){list(sigma.proc = runif(1, 0, 5), mean.lambda = runif(1, 0.1, 2), sigma.obs = runif(1, 0, 10), N.est = c(runif(1, 30, 60), rep(NA, (n.years-1))))}

# Parameters monitored
parameters <- c("lambda", "mean.lambda", "sigma2.obs", "sigma2.proc", "N.est")

# MCMC settings
ni <- 25000
nt <- 3
nb <- 10000
nc <- 3

# Call JAGS from R (BRT <1 min)
ssm <- jags(jags.data, inits, parameters, "ssm.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

# Summarize posteriors
print(ssm, digits = 3)

# Produce figure
graph.ssm(ssm, N, y)



n.years <- 25  # Number of years
N <- rep(50, n.years)

lp <- -0.5 + 0.1*(1:n.years)  # Increasing trend of logit p
p <- plogis(lp)
y <- numeric(n.years)
for (t in 1:n.years){
   y[t] <- rbinom(1, N[t], p[t])
   }

# Bundle data
jags.data <- list(y = y, T = n.years)

# Call JAGS from R (BRT <1 min)
ssm <- jags(jags.data, inits, parameters, "ssm.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

# Produce figure
graph.ssm(ssm, N, y)
points(N*p, col = "black", type = "l", lwd = 2, lty = 2)
legend(x = 1, y = 45.5, legend = "Np", lwd = 2, col = "black", lty = 2, bty = "n")


# 5.4. Real example: House martin population counts in the village of Magden
# Specify model in BUGS language
sink("ssm.jags")
cat("
model {
# Priors and constraints
logN.est[1] ~ dnorm(5.6, 0.01)       # Prior for initial population size
mean.r ~ dnorm(1, 0.001)             # Prior for mean growth rate
sigma.proc ~ dunif(0, 1)             # Prior for sd of state process
sigma2.proc <- pow(sigma.proc, 2)
tau.proc <- pow(sigma.proc, -2)
sigma.obs ~ dunif(0, 1)              # Prior for sd of observation process
sigma2.obs <- pow(sigma.obs, 2)
tau.obs <- pow(sigma.obs, -2)

# Likelihood
# State process
for (t in 1:(T-1)){
   r[t] ~ dnorm(mean.r, tau.proc)
   logN.est[t+1] <- logN.est[t] + r[t]
   }
# Observation process
for (t in 1:T) {
   y[t] ~ dnorm(logN.est[t], tau.obs)
   }

# Population sizes on real scale
for (t in 1:T) {
   N.est[t] <- exp(logN.est[t])
   }
}
",fill = TRUE)
sink()

# House martin population data from Magden
pyears <- 6 # Number of future years with predictions
hm <- c(271, 261, 309, 318, 231, 216, 208, 226, 195, 226, 233, 209, 226, 192, 191, 225, 245, 205, 191, 174, rep(NA, pyears))
year <- 1990:(2009 + pyears)

# Bundle data
jags.data <- list(y = log(hm), T = length(year))

# Initial values
inits <- function(){list(sigma.proc = runif(1, 0, 1), mean.r = rnorm(1), sigma.obs = runif(1, 0, 1), logN.est = c(rnorm(1, 5.6, 0.1), rep(NA, (length(year)-1))))}

# Parameters monitored
parameters <- c("r", "mean.r", "sigma2.obs", "sigma2.proc", "N.est")

# MCMC settings
ni <- 200000
nt <- 6
nb <- 100000
nc <- 3

# Call JAGS from R (BRT 3 min)
hm.ssm <- jags(jags.data, inits, parameters, "ssm.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

# Summarize posteriors
print(hm.ssm, digits = 3)


# Draw figure
fitted <- lower <- upper <- numeric()
year <- 1990:2015
n.years <- length(hm)
for (i in 1:n.years){
   fitted[i] <- mean(hm.ssm$BUGSoutput$sims.list$N.est[,i])
   lower[i] <- quantile(hm.ssm$BUGSoutput$sims.list$N.est[,i], 0.025)
   upper[i] <- quantile(hm.ssm$BUGSoutput$sims.list$N.est[,i], 0.975)}
m1 <- min(c(fitted, hm, lower), na.rm = TRUE)
m2 <- max(c(fitted, hm, upper), na.rm = TRUE)
par(mar = c(4.5, 4, 1, 1))
plot(0, 0, ylim = c(m1, m2), xlim = c(1, n.years), ylab = "Population size", xlab = "Year", col = "black", type = "l", lwd = 2, axes = FALSE, frame = FALSE)
axis(2, las = 1)
axis(1, at = 1:n.years, labels = year)
polygon(x = c(1:n.years, n.years:1), y = c(lower, upper[n.years:1]), col = "gray90", border = "gray90")
points(hm, type = "l", col = "black", lwd = 2)
points(fitted, type = "l", col = "blue", lwd = 2)
legend(x = 1, y = 150, legend = c("Counts", "Estimates"), lty = c(1, 1), lwd = c(2, 2), col = c("black", "blue"), bty = "n", cex = 1)

# Probability of N(2015) < N(2009)
mean(hm.ssm$BUGSoutput$sims.list$N.est[,26] < hm.ssm$BUGSoutput$mean$N.est[20])




#######################################################################
#
# 6. Estimation of the size of a closed population
# 
#######################################################################

# 6.2. Generation and analysis of simulated data with data augmentation
# 6.2.1. Introduction to data augmentation for the simplest case: model M0
# Define function to simulate data under M0
data.fn <- function(N = 100, p = 0.5, T = 3){
   yfull <- yobs <- array(NA, dim = c(N, T))
   for (j in 1:T){
      yfull[,j] <- rbinom(n = N, size = 1, prob = p)
      }
   ever.detected <- apply(yfull, 1, max)
   C <- sum(ever.detected)
   yobs <- yfull[ever.detected == 1,]
   cat(C, "out of", N, "animals present were detected.\n")
   return(list(N = N, p = p, C = C, T = T, yfull = yfull, yobs = yobs))
   }

data <- data.fn()

str(data)

# Augment data set by 150 potential individuals
nz <- 150
yaug <- rbind(data$yobs, array(0, dim = c(nz, data$T)))

# Specify model in BUGS language
sink("model.jags")
cat("
model {

# Priors
omega ~ dunif(0, 1)
p ~ dunif(0, 1)

# Likelihood
for (i in 1:M){
   z[i] ~ dbern(omega)			# Inclusion indicators
   for (j in 1:T){
      yaug[i,j] ~ dbern(p.eff[i,j])
      p.eff[i,j] <- z[i] * p		# Can only be detected if z=1
      } #j
   } #i

# Derived quantities
N <- sum(z[])
}
",fill = TRUE)
sink()

# Bundle data
win.data <- list(yaug = yaug, M = nrow(yaug), T = ncol(yaug))

# Initial values
inits <- function() list(z = rep(1, nrow(yaug)), p = runif(1, 0, 1))

# Parameters monitored
params <- c("N", "p", "omega")

# MCMC settings
ni <- 2500
nt <- 2
nb <- 500
nc <- 3

# Call JAGS from R (BRT <1 min)
out <- jags(win.data, inits, params, "model.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

# Summarize posteriors
print(out, dig = 3)
hist(out$BUGSoutput$sims.list$N, nclass = 50, col = "gray", main = "", xlab = "Population size", las = 1, xlim = c(80, 150))
abline(v = data$C, lwd = 3)

nz <- 5
yaug <- rbind(data$yobs, array(0, dim = c(nz, data$T)))
win.data <- list(yaug = yaug, M = dim(yaug)[1], T = dim(yaug)[2])
out5 <- jags(win.data, inits, params, "model.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())
print(out5, dig = 3)
par(mfrow = c(3, 1))
hist(out5$BUGSoutput$sims.list$N, nclass = 30, col = "gray", main = "Augmentation by 5", xlab = "Population size", las = 1, xlim = c(80, 140))
abline(v = data$C, col = "black", lwd = 3)
abline(v = mean(out5$BUGSoutput$sims.list$N), col = "blue", lwd = 3)

nz <- 150
yaug <- rbind(data$yobs, array(0, dim = c(nz, data$T)))
win.data <- list(yaug = yaug, M = dim(yaug)[1], T = dim(yaug)[2])
out150 <- jags(win.data, inits, params, "model.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())
print(out150, dig = 3)
hist(out$BUGSoutput$sims.list$N, nclass = 30, col = "gray", main = "Augmentation by 150", xlab = "Population size", las = 1, xlim = c(80, 140))
abline(v = data$C, col = "black", lwd = 3)
abline(v = mean(out150$BUGSoutput$sims.list$N), col = "blue", lwd = 3)

nz <- 1500
yaug <- rbind(data$yobs, array(0, dim = c(nz, data$T)))
win.data <- list(yaug = yaug, M = dim(yaug)[1], T = dim(yaug)[2])
out1500 <- jags(win.data, inits, params, "model.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())
print(out1500, dig = 3)
hist(out1500$BUGSoutput$sims.list$N, nclass = 30, col = "gray", main = "Augmentation by 1500", xlab = "Population size", las = 1, xlim = c(80, 140))
abline(v = data$C, col = "black", lwd = 3)
abline(v = mean(out1500$BUGSoutput$sims.list$N), col = "blue", lwd = 3)


# 6.2.2. Time effects: model Mt
# Define function to simulate data under Mt
data.fn <- function(N = 100, mean.p = 0.5, T = 3, time.eff = runif(T, -2, 2)){
   yfull <- yobs <- array(NA, dim = c(N, T) )
   p.vec <- array(NA, dim = T)
   for (j in 1:T){
      p <- plogis(log(mean.p / (1-mean.p)) + time.eff[j])
      yfull[,j] <- rbinom(n = N, size = 1, prob = p)
      p.vec[j] <- p
      }
   ever.detected <- apply(yfull, 1, max)
   C <- sum(ever.detected)
   yobs <- yfull[ever.detected == 1,]
   cat(C, "out of", N, "animals present were detected.\n")
   return(list(N = N, p.vec = p.vec, C = C, T = T, yfull = yfull, yobs = yobs))
   }

data <- data.fn()

# Augment data set
nz <- 150
yaug <- rbind(data$yobs, array(0, dim = c(nz, data$T)))

# Specify model in BUGS language
sink("model.jags")
cat("
model {
# Priors
omega ~ dunif(0, 1)
for (i in 1:T){
   p[i] ~ dunif(0, 1)
   }

# Likelihood
for (i in 1:M){
   z[i] ~ dbern(omega)
   for (j in 1:T){
      yaug[i,j] ~ dbern(p.eff[i,j])
      p.eff[i,j] <- z[i] * p[j]
      } #j
   } #i

# Derived quantities
N <- sum(z[])
} # end model
",fill = TRUE)
sink()

# Bundle data
win.data <- list(yaug = yaug, M = nrow(yaug), T = ncol(yaug))

# Initial values
inits <- function() list(z = rep(1, nrow(yaug)), p = runif(data$T, 0, 1))

# Parameters monitored
params <- c("N", "p", "omega")

# MCMC settings
ni <- 2500
nt <- 2
nb <- 500
nc <- 3

# Call JAGS from R (BRT <1 min)
out <- jags(win.data, inits, params, "model.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

# Summarize posteriors
print(out, dig = 3)
hist(out$BUGSoutput$sims.list$N, nclass = 40, col = "gray", main = "", xlab = "Population size", las = 1, xlim = c(70, 150))
abline(v = data$C, col = "black", lwd = 3)


# 6.2.3. Behavioral or memory effects: model Mb
# Define function to simulate data under Mb
data.fn <- function(N = 200, T = 5, p = 0.3, c = 0.4){
   yfull <- yobs <- array(NA, dim = c(N, T) )
   p.eff <- array(NA, dim = N)

   # First capture occasion
   yfull[,1] <- rbinom(n = N, size = 1, prob = p)

   # Later capture occasions
   for (j in 2:T){
      p.eff <- (1 - yfull[,(j-1)]) * p + yfull[,(j-1)] * c
      yfull[,j] <- rbinom(n = N, size = 1, prob = p.eff)
      }

   ever.detected <- apply(yfull, 1, max)
   C <- sum(ever.detected)
   yobs <- yfull[ever.detected == 1,]
   cat(C, "out of", N, "animals present were detected.\n")
   return(list(N = N, p = p, c = c, C = C, T = T, yfull = yfull, yobs = yobs))
   }

data <- data.fn(N = 200)

# Augment data set
nz <- 150
yaug <- rbind(data$yobs, array(0, dim = c(nz, data$T)))

# Specify model in BUGS language
sink("model.jags")
cat("
model {
# Priors
omega ~ dunif(0, 1)
p ~ dunif(0, 1)     # Cap prob when caught at t-1
c ~ dunif(0, 1)     # Cap prob when not caught at t-1

# Likelihood
for (i in 1:M){
   z[i] ~ dbern(omega)

   # First occasion
   yaug[i,1] ~ dbern(p.eff[i,1])
   p.eff[i,1] <- z[i] * p

   # All subsequent occasions
   for (j in 2:T){
      yaug[i,j] ~ dbern(p.eff[i,j])
      p.eff[i,j] <- z[i] * ( (1-yaug[i,(j-1)]) * p + yaug[i,(j-1)] * c )
      } #j
   } #i

# Derived quantities
N <- sum(z[])
trap.response <- c - p
} # end model
",fill = TRUE)
sink()

# Bundle data
win.data <- list(yaug = yaug, M = nrow(yaug), T = ncol(yaug))

# Initial values
inits <- function() list(z = rep(1, nrow(yaug)), p = runif(1, 0, 1))

# Parameters monitored
params <- c("N", "p", "c", "trap.response", "omega")

# MCMC settings
ni <- 2500
nt <- 2
nb <- 500
nc <- 3

# Call JAGS from R (BRT <1 min)
out <- jags(win.data, inits, params, "model.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

# Summarize posteriors
print(out, dig = 3)

hist(out$BUGSoutput$sims.list$N, nclass = 40, col = "gray", main = "", xlab = "Population size", las = 1, xlim = c(150, 300))
abline(v= data$C, col = "black", lwd = 3)


# 6.2.4. Individual (random) effects: the heterogeneity model Mh
# Define function to simulate data under Mh
data.fn <- function(N = 100, mean.p = 0.4, T = 5, sd = 1){
   yfull <- yobs <- array(NA, dim = c(N, T) )
   mean.lp <- log(mean.p / (1-mean.p))
   p.vec <- plogis(mean.lp+ rnorm(N, 0, sd))

   for (i in 1:N){
      yfull[i,] <- rbinom(n = T, size = 1, prob = p.vec[i])
      }

   ever.detected <- apply(yfull, 1, max)
   C <- sum(ever.detected)
   yobs <- yfull[ever.detected == 1,]
   cat(C, "out of", N, "animals present were detected.\n")
   hist(p.vec, xlim = c(0,1), nclass = 20, col = "gray", main = "", xlab = "Detection probability", las = 1)
   return(list(N = N, p.vec = p.vec, mean.lp = mean.lp, C = C, T = T, yfull = yfull, yobs = yobs))
   }

data <- data.fn()

# Aggregate capture-histories and augment data set
y <- sort(apply(data$yobs, 1, sum), decreasing = TRUE)
nz <- 300
yaug <- c(y, rep(0, nz))
yaug

# Specify model in BUGS language
sink("model.jags")
cat("
model {

# Priors
omega ~ dunif(0, 1)
mean.lp <- logit(mean.p)
mean.p ~ dunif(0, 1)
tau <- 1 / (sd * sd)
sd ~ dunif(0, 5)

# Likelihood
for (i in 1:M){
   z[i] ~ dbern(omega)
   logit(p[i]) <- eps[i]
   eps[i] ~ dnorm(mean.lp, tau)T(-16, 16)	# See web appendix A in Royle (2009): but JAGS needs T instead of I
   p.eff[i] <- z[i] * p[i]
   y[i] ~ dbin(p.eff[i], T)
   }

# Derived quantities
N <- sum(z[])
}
",fill = TRUE)
sink()

# Bundle data
win.data <- list(y = yaug, M = length(yaug), T = ncol(data$yobs))

# Initial values
inits <- function() list(z = rep(1, length(yaug)), sd = runif(1, 0.1, 0.9))

# Parameters monitored
params <- c("N", "mean.p", "sd", "omega")

# MCMC settings
ni <- 25000
nt <- 2
nb <- 5000
nc <- 3

# Call JAGS from R (BRT 6 min)
out <- jags(win.data, inits, params, "model.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

# Summarize posteriors
print(out, dig = 3)
hist(out$BUGSoutput$sims.list$N, nclass = 50, col = "gray", main = "", xlab = "Population size", las = 1, xlim = c(80, 250))
abline(v = data$C, col = "black", lwd = 3)


# 6.2.5. Combined effects: model Mth
# Define function to simulate data under Mth
data.fn <- function(N = 100, T = 5, mean.p = 0.4, time.effects = runif(T, -1, 1), sd = 1){
   yfull <- yobs <- p <- array(NA, dim = c(N, T) )
   mean.lp <- log(mean.p / (1-mean.p))         # mean p on logit scale
   eps <- rnorm(N, 0, sd)                      # Individual effects

   for (j in 1:T){
      pp <- p[,j] <- plogis(mean.lp + time.effects[j] + eps)
      yfull[,j] <- rbinom(n = N, size = 1, prob = pp)
      }

   ever.detected <- apply(yfull, 1, max)
   C <- sum(ever.detected)
   yobs <- yfull[ever.detected == 1,]
   cat(C, "out of", N, "animals present were detected.\n")
   cat("Mean p per occasion:", round(apply(p, 2, mean), 2), "\n")
   par(mfrow = c(2,1))
   plot(plogis(mean.lp + time.effects), xlab = "Occasion", type = "b", main = "Approx. mean p at each occasion", ylim = c(0, 1))
   hist(plogis(mean.lp + eps), xlim = c(0, 1), col = "gray", main = "Approx. distribution of p at average occasion")
   return(list(N = N, mean.lp = mean.lp, time.effects = time.effects, sd = sd, eps = eps, C = C, T = T, yfull = yfull, yobs = yobs))
   }

data <- data.fn()

# data<-data.fn(T = 10, mean.p = 0.2, time.effects = runif(10, 0, 0), sd = 0)	# M0
# data<-data.fn(T = 10, mean.p = 0.5, time.effects = runif(10, 0, 0), sd = 1)	# Mh
# data <- data.fn(T = 10, sd = 0)	# Mt

# Augment data set
nz <- 300
yaug <- rbind(data$yobs, array(0, dim=c(nz, data$T)))

# Specify model in BUGS language
sink("model.jags")
cat("
model {

# Priors
omega ~ dunif(0, 1)
for (j in 1:T){
   mean.lp[j] <- log(mean.p[j] / (1 - mean.p[j]) ) # Define logit 
   mean.p[j] ~ dunif(0, 1)
   }
tau <- 1 / (sd * sd)
sd ~ dunif(0,5)

# Likelihood
for (i in 1:M){
   z[i] ~ dbern(omega)
   eps[i] ~ dnorm(0, tau)T(-16, 16)	# See web appendix A in Royle (2009): but JAGS needs T instead of I
   for (j in 1:T){
      lp[i,j] <- mean.lp[j] + eps[i]
      p[i,j] <- 1 / (1 + exp(-lp[i,j])) 			# Define logit
      p.eff[i,j] <- z[i] * p[i,j]
      y[i,j] ~ dbern(p.eff[i,j])
      } #j
   } #i

# Derived quantities
N <- sum(z[])
} 
",fill = TRUE)
sink()

# Bundle data
win.data <- list(y = yaug, M = nrow(yaug), T = ncol(yaug))

# Initial values
inits <- function() list(z = rep(1, nrow(yaug)), sd = runif(1, 0.1, 0.9))

# Parameters monitored
params <- c("N", "mean.p", "mean.lp", "sd", "omega")

# MCMC settings
ni <- 25000
nt <- 2
nb <- 5000
nc <- 3

# Call JAGS from R (BRT 47 min)
out <- jags(win.data, inits, params, "model.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

# Summarize posteriors
print(out, dig = 3)
hist(out$BUGSoutput$sims.list$N, nclass = 50, col = "gray", main = "", xlab = "Population size", las = 1, xlim = c(80, 200))
abline(v = data$C, col = "black", lwd = 3)


# 6.3. Analysis of a real data set: model Mtbh for species richness estimation
# Read in data and look at them
p610 <- read.table("p610.txt", header = TRUE)
y <- p610[,5:9]                           # Grab counts
y[y > 1] <- 1                             # Counts to det-nondetections
C <- sum(apply(y, 1, max)) ; print(C)     # Number of observed species
table(apply(y, 1, sum))                   # Capture-frequencies

# Specify model in BUGS language
sink("M_tbh.jags")
cat("
model {

# Priors
omega ~ dunif(0, 1)
for  (j in 1:T){
   alpha[j] <- log(mean.p[j] / (1-mean.p[j])) # Define logit 
   mean.p[j] ~ dunif(0, 1) 	# Detection intercepts
   }
gamma ~ dnorm(0, 0.01)
tau <- 1 / (sd * sd)
sd ~ dunif(0, 3)

# Likelihood
for (i in 1:M){
   z[i] ~ dbern(omega)
   eps[i] ~ dnorm(0, tau)T(-16, 16)

   # First occasion: no term for recapture (gamma)
   y[i,1] ~ dbern(p.eff[i,1])
   p.eff[i,1] <- z[i] * p[i,1]
   p[i,1] <- 1 / (1 + exp(-lp[i,1]))
   lp[i,1] <- alpha[1] + eps[i]

   # All subsequent occasions: includes recapture term (gamma)
   for (j in 2:T){
      y[i,j] ~ dbern(p.eff[i,j])
      p.eff[i,j] <- z[i] * p[i,j]
      p[i,j] <- 1 / (1 + exp(-lp[i,j]))   
      lp[i,j] <- alpha[j] + eps[i] + gamma * y[i,(j-1)]
      } #j
   } #i

# Derived quantities
N <- sum(z[])
} 
",fill = TRUE)
sink()

# Bundle data
win.data <- list(y = as.matrix(y), M = nrow(y), T = ncol(y))

# Initial values
inits <- function() list(z = rep(1, nrow(y)), sd = runif(1, 0.1, 0.9))

# Parameters monitored
params <- c("N", "mean.p", "gamma", "sd", "omega")

# MCMC settings
ni <- 50000
nt <- 4
nb <- 10000
nc <- 3

# Call JAGS from R (BRT 24 min)
out <- jags(win.data, inits, params, "M_tbh.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

# Summarize posteriors and plot posterior for N
print(out, dig = 3)
par(mfrow = c(1,2))
hist(out$BUGSoutput$sims.list$N, breaks = 35, col = "gray", main = "", xlab = "Community size", las = 1, xlim = c(30, 100), freq = FALSE)
abline(v = C, col = "black", lwd = 3)

# Define model
sink("M0.jags")
cat("
model {

# Priors
omega ~ dunif(0, 1)
p ~ dunif(0, 1)

# Likelihood
for (i in 1:M){
   z[i] ~ dbern(omega)
   for (j in 1:T){
      y[i,j] ~ dbern(p.eff[i,j])
      p.eff[i,j] <- z[i] * p
      } #j
   } #i

# Derived quantities
N <- sum(z[])
} # end model
",fill = TRUE)
sink()

# Initial values
inits <- function() list(z = rep(1, nrow(y)))

# Define parameters to be monitored
params <- c("N", "p", "omega")

# MCMC settings
ni <- 50000
nt <- 4
nb <- 10000
nc <- 3

# Call JAGS from R (BRT 1 min)
out0 <- jags(win.data, inits, params, "M0.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,  working.directory = getwd())

# Inspect output
print(out0, dig = 3)


# 6.4. Capture-recapture models with individual covariates: model Mt+X
# 6.4.1. Individual covariate model for species richness estimation
p610 <- read.table("p610.txt", header = TRUE)
y <- p610[,5:9]                         # Grab counts
y[y > 1] <- 1                           # Convert to det-nondetections
ever.observed <- apply(y, 1, max)
wt <- p610$bm[ever.observed == 1]       # Body mass
yy <- as.matrix(y[ever.observed == 1,]) # Detection histories
dimnames(yy) <- NULL

mlog <- mean(log(p610$bm^(1/3)))
sdlog <- sd(log(p610$bm^(1/3)))
hist(p610$bm^(1/3), xlim = c(0, 30), nclass = 25, freq = FALSE, col = "gray")
lines(density(rlnorm(n = 10^6, meanlog = mlog, sdlog = sdlog)), col = "blue", lwd = 3)

# Augment both data sets
nz = 150
yaug <- rbind(yy, array(0, dim = c(nz, ncol(yy))))
logwt3 <- c(log(wt^(1/3)), rep(NA, nz))

# Specify model in BUGS language
sink("M_t+X.jags")
cat("
model {

# Priors
omega ~ dunif(0, 1)
for (j in 1:T){
   alpha[j] <- log(mean.p[j] / (1-mean.p[j]))
   mean.p[j] ~ dunif(0, 1)
   }
beta ~ dnorm(0, 0.01)
mu.size ~ dnorm(0, 0.01)
tau.size <- 1 / pow(sd.size, 2)
sd.size ~ dunif(0, prior.sd.upper)   # Provide upper bound as data

# Likelihood
for (i in 1:M){  # Loop over individuals
   z[i] ~ dbern(omega)
   size[i] ~ dnorm(mu.size, tau.size)T(-6,6)
   for (j in 1:T){  # Loop over occasions
      y[i,j] ~ dbern(p.eff[i,j])
      p.eff[i,j] <- z[i] * p[i,j]
      p[i,j] <- 1 / (1 + exp(-lp[i,j]))   
      lp[i,j] <- alpha[j] + beta * size[i]
      } #j
   } #i

# Derived quantities
N <- sum(z[])
}
",fill = TRUE)
sink()

# Bundle data
win.data <- list(y = yaug, size = logwt3 - mean(logwt3, na.rm = TRUE), M = nrow(yaug), T = ncol(yaug), prior.sd.upper = 3)

# Initial values
inits <- function() list(z = rep(1, nrow(yaug)), beta = runif(1, 0, 1), mu.size = rnorm(1, 0, 1))

# Parameters monitored
params <- c("N", "mean.p", "beta", "omega", "mu.size", "sd.size")

# MCMC settings
ni <- 50000
nt <- 4
nb <- 10000
nc <- 3

# Call JAGS from R (BRT 19 min)
outX <- jags(win.data, inits, params, "M_t+X.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

# Summarize posteriors and plot posterior for N
print(outX, dig = 3)

hist(outX$BUGSoutput$sims.list$N, breaks = 100, col = "gray", main = "", xlab = "Community size", las = 1, xlim = c(30, 100), freq = FALSE)
abline(v = 31, col = "black", lwd = 3)

pred.wt <- seq(5, 2000, length.out = 100)	# Cov. vals for prediction
pred.wt.st <- log(pred.wt^(1/3))- mean(logwt3, na.rm = TRUE) # Transform them in the same was as in the analysis
pred.p<- plogis(log(mean(outX$BUGSoutput$mean$mean.p)/(1- mean(outX$BUGSoutput$mean$mean.p))) + outX$BUGSoutput$mean$beta * pred.wt.st) # Compute predicted response
plot(pred.wt, pred.p, type = "l", lwd = 3, col = "blue", las = 1, frame.plot = FALSE, ylim = c(0, 0.5))


# 6.4.2. Individual covariate model for population size estimation
# Read in data and look at shell width distribution
pinna <- read.table("pinna.txt", header = TRUE)
y <- cbind(pinna$d1, pinna$d2)
size <- pinna$width
hist(size, col = "gray", nclass = 50, xlim = c(0, 30), freq = FALSE)
lines(density(rnorm(10^6, mean = mean(size), sd = sd(size))), col = "blue", lwd = 3)

# Augment both data sets
nz = 150
yaug <- rbind(y, array(0, dim = c(nz, ncol(y))))
size <- c(size, rep(NA, nz))

# Bundle data
win.data <- list(y = yaug, size = size - mean(size, na.rm = TRUE), M = nrow(yaug), T = ncol(yaug), prior.sd.upper = 5)

# Specify model in BUGS language. Is exactly the same model as before, but without a truncation of the normal distribution for size
sink("M_t+X.jags")
cat("
model {

# Priors
omega ~ dunif(0, 1)
for (j in 1:T){
   alpha[j] <- log(mean.p[j] / (1-mean.p[j]))
   mean.p[j] ~ dunif(0, 1)
   }
beta ~ dnorm(0, 0.01)
mu.size ~ dnorm(0, 0.01)
tau.size <- 1 / pow(sd.size, 2)
sd.size ~ dunif(0, prior.sd.upper)   # Provide upper bound as data

# Likelihood
for (i in 1:M){  # Loop over individuals
   z[i] ~ dbern(omega)
   size[i] ~ dnorm(mu.size, tau.size)
   for (j in 1:T){  # Loop over occasions
      y[i,j] ~ dbern(p.eff[i,j])
      p.eff[i,j] <- z[i] * p[i,j]
      p[i,j] <- 1 / (1 + exp(-lp[i,j]))   
      lp[i,j] <- alpha[j] + beta * size[i]
      } #j
   } #i

# Derived quantities
N <- sum(z[])
}
",fill = TRUE)
sink()


# MCMC settings
ni <- 2500
nt <- 2
nb <- 500
nc <- 3

# Call JAGS from R (BRT 1 min)
outXX <- jags(win.data, inits, params, "M_t+X.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

# Summarize posteriors
print(outXX, dig = 2)

# Plot posterior for N and prediction of p
par(mfrow = c(1,2), mar = c(4.5, 4, 2, 1))
hist(outXX$BUGSoutput$sims.list$N, breaks = 30, col = "gray", main = "", xlab = "Population size", las = 1, xlim = c(143, 220), freq = FALSE)
abline(v = 143, col = "black", lwd = 3)

pred.size <- seq(0, 30, length.out = 1000)	# Cov. vals for prediction
pred.size.st <- pred.size - mean(size, na.rm = TRUE) # Transform them
pred.p<- plogis(log(mean(outXX$BUGSoutput$mean$mean.p)/(1- mean(outXX$BUGSoutput$mean$mean.p))) + outXX$BUGSoutput$mean$beta * pred.size.st) # Compute predicted detection prob.
plot(pred.size, pred.p, type = "l", lwd = 3, col = "blue", las = 1, frame.plot = FALSE, ylim = c(0, 1), xlab = "Shell width (cm)", ylab = "Predicted detection probability")




####################################################################
#
# 7. Estimation of survival probabilities using capture-recapture data
#
#####################################################################

# 7.3. Models with constant parameters
# Define parameter values
n.occasions <- 6                   # Number of capture occasions
marked <- rep(50, n.occasions-1)   # Annual number of newly marked individuals
phi <- rep(0.65, n.occasions-1)
p <- rep(0.4, n.occasions-1)

# Define matrices with survival and recapture probabilities
PHI <- matrix(phi, ncol = n.occasions-1, nrow = sum(marked))
P <- matrix(p, ncol = n.occasions-1, nrow = sum(marked))

# Define function to simulate a capture-history (CH) matrix
simul.cjs <- function(PHI, P, marked){
   n.occasions <- dim(PHI)[2] + 1
   CH <- matrix(0, ncol = n.occasions, nrow = sum(marked))
   # Define a vector with the occasion of marking
   mark.occ <- rep(1:length(marked), marked[1:length(marked)])
   # Fill the CH matrix
   for (i in 1:sum(marked)){
      CH[i, mark.occ[i]] <- 1       # Write an 1 at the release occasion
      if (mark.occ[i]==n.occasions) next
      for (t in (mark.occ[i]+1):n.occasions){
         # Bernoulli trial: does individual survive occasion?
         sur <- rbinom(1, 1, PHI[i,t-1])
         if (sur==0) break		# If dead, move to next individual 
         # Bernoulli trial: is individual recaptured? 
         rp <- rbinom(1, 1, P[i,t-1])
         if (rp==1) CH[i,t] <- 1
         } #t
      } #i
   return(CH)
   }

# Execute function
CH <- simul.cjs(PHI, P, marked)

# Create vector with occasion of marking
get.first <- function(x) min(which(x!=0))
f <- apply(CH, 1, get.first)

# Specify model in BUGS language
sink("cjs-c-c.jags")
cat("
model {

# Priors and constraints
for (i in 1:nind){
   for (t in f[i]:(n.occasions-1)){
      phi[i,t] <- mean.phi
      p[i,t] <- mean.p
      } #t
   } #i

mean.phi ~ dunif(0, 1)         # Prior for mean survival
mean.p ~ dunif(0, 1)           # Prior for mean recapture

# Likelihood 
for (i in 1:nind){
   # Define latent state at first capture
   z[i,f[i]] <- 1
   for (t in (f[i]+1):n.occasions){
      # State process
      z[i,t] ~ dbern(mu1[i,t])
      mu1[i,t] <- phi[i,t-1] * z[i,t-1]
      # Observation process
      y[i,t] ~ dbern(mu2[i,t])
      mu2[i,t] <- p[i,t-1] * z[i,t]
      } #t
   } #i
}
",fill = TRUE)
sink()

# Bundle data
jags.data <- list(y = CH, f = f, nind = dim(CH)[1], n.occasions = dim(CH)[2])


# Initial values
# In JAGS we have to give good initial values for the latent state z. At all occasions when an individual was observed, its state is z = 1 for sure. In addition, if an individual was not observed at an occasion, but was alive for sure, because it was observed before and thereafter (i.e. has a capture history of e.g. {101} or {10001}), then we know that the individual was alive at all of these occasions, and thus z = 1. Therefore, we should provide initial values of z = 1 at these positions as well. The following function provides such initial values from the observed capture histories:
known.state.cjs <- function(ch){
   state <- ch
   for (i in 1:dim(ch)[1]){
      n1 <- min(which(ch[i,]==1))
      n2 <- max(which(ch[i,]==1))
      state[i,n1:n2] <- 1
      state[i,n1] <- NA
      }
   state[state==0] <- NA
   return(state)
   }

# (Note that the function known.state.cjs is used in section 7.3.1 as well for another purpose) 

inits <- function(){list(mean.phi = runif(1, 0, 1), mean.p = runif(1, 0, 1), z = known.state.cjs(CH))}

# Parameters monitored
parameters <- c("mean.phi", "mean.p")

# MCMC settings
ni <- 10000
nt <- 6
nb <- 5000
nc <- 3

# Call JAGS from R (BRT 1 min)
cjs.c.c <- jags(jags.data, inits, parameters, "cjs-c-c.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

# Summarize posteriors
print(cjs.c.c, digits = 3)


# 7.3.1. Inclusion of information about latent state variable
# Function to create a matrix with information about known latent state z
known.state.cjs <- function(ch){
   state <- ch
   for (i in 1:dim(ch)[1]){
      n1 <- min(which(ch[i,]==1))
      n2 <- max(which(ch[i,]==1))
      state[i,n1:n2] <- 1
      state[i,n1] <- NA
      }
   state[state==0] <- NA
   return(state)
   }

# Bundle data
jags.data <- list(y = CH, f = f, nind = dim(CH)[1], n.occasions = dim(CH)[2], z = known.state.cjs(CH))

# Function to create a matrix of initial values for latent state z
cjs.init.z <- function(ch,f){
   for (i in 1:dim(ch)[1]){
      if (sum(ch[i,])==1) next
      n2 <- max(which(ch[i,]==1))
      ch[i,f[i]:n2] <- NA
      }
   for (i in 1:dim(ch)[1]){
   ch[i,1:f[i]] <- NA
   }
   return(ch)
   }

# Initial values
inits <- function(){list(z = cjs.init.z(CH, f), mean.phi = runif(1, 0, 1), mean.p = runif(1, 0, 1))}  

# Parameters monitored
parameters <- c("mean.phi", "mean.p")

# MCMC settings
ni <- 10000
nt <- 6
nb <- 5000
nc <- 3

# Call JAGS from R (BRT <1 min)
cjs.c.c <- jags(jags.data, inits, parameters, "cjs-c-c.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

# Summarize posteriors
print(cjs.c.c, digits = 3)


# 7.4. Models with time-variation
# 7.4.1. Fixed time effects
# Priors and constraints
for (i in 1:nind){
   for (t in f[i]:(n.occasions-1)){
      phi[i,t] <- alpha[t]
      p[i,t] <- beta[t]
      } #t
   } #i
for (t in 1:(n.occasions-1)){
   alpha[t] ~ dunif(0, 1)        # Priors for time-spec. survival
   beta[t] ~ dunif(0, 1)         # Priors for time-spec. recapture
   }


# 7.4.2. Random time effects
# Define parameter values
n.occasions <- 20                  # Number of capture occasions
marked <- rep(30, n.occasions-1)   # Annual number of newly marked individuals
mean.phi <- 0.65
var.phi <- 1                       # Temporal variance of survival
p <- rep(0.4, n.occasions-1)

# Determine annual survival probabilities
logit.phi <- rnorm(n.occasions-1, qlogis(mean.phi), var.phi^0.5)
phi <- plogis(logit.phi)

# Define matrices with survival and recapture probabilities
PHI <- matrix(phi, ncol = n.occasions-1, nrow = sum(marked), byrow = TRUE)
P <- matrix(p, ncol = n.occasions-1, nrow = sum(marked))

# Simulate capture-histories
CH <- simul.cjs(PHI, P, marked)

# Create vector with occasion of marking
get.first <- function(x) min(which(x!=0))
f <- apply(CH, 1, get.first)

# Specify model in BUGS language
sink("cjs-temp-raneff.jags")
cat("
model {

# Priors and constraints
for (i in 1:nind){
   for (t in f[i]:(n.occasions-1)){
      logit(phi[i,t]) <- mu + epsilon[t]
      p[i,t] <- mean.p
      } #t
   } #i
for (t in 1:(n.occasions-1)){
   epsilon[t] ~ dnorm(0, tau)
   }

#mu ~ dnorm(0, 0.001)                    # Prior for logit of mean survival
#mean.phi <- 1 / (1+exp(-mu))            # Logit transformation
mean.phi ~ dunif(0, 1)                   # Prior for mean survival
mu <- log(mean.phi / (1-mean.phi))       # Logit transformation
sigma ~ dunif(0, 10)                     # Prior for standard deviation
tau <- pow(sigma, -2)
sigma2 <- pow(sigma, 2)                  # Temporal variance
mean.p ~ dunif(0, 1)                     # Prior for mean recapture

# Likelihood 
for (i in 1:nind){
   # Define latent state at first capture
   z[i,f[i]] <- 1
   for (t in (f[i]+1):n.occasions){
      # State process
      z[i,t] ~ dbern(mu1[i,t])
      mu1[i,t] <- phi[i,t-1] * z[i,t-1]
      # Observation process
      y[i,t] ~ dbern(mu2[i,t])
      mu2[i,t] <- p[i,t-1] * z[i,t]
      } #t
   } #i
}
",fill = TRUE)
sink()

# Bundle data
jags.data <- list(y = CH, f = f, nind = dim(CH)[1], n.occasions = dim(CH)[2], z = known.state.cjs(CH))

# Initial values
inits <- function(){list(z = cjs.init.z(CH, f), mean.phi = runif(1, 0, 1), sigma = runif(1, 0, 10), mean.p = runif(1, 0, 1))}  

# Parameters monitored
parameters <- c("mean.phi", "mean.p", "sigma2")

# MCMC settings
ni <- 10000
nt <- 6
nb <- 5000
nc <- 3

# Call JAGS from R (BRT 17 min)
cjs.ran <- jags(jags.data, inits, parameters, "cjs-temp-raneff.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

# Summarize posteriors
print(cjs.ran, digits = 3)

# Produce histogram
hist(cjs.ran$BUGSoutput$sims.list$sigma2, col = "gray", nclass = 35, las = 1, xlab = expression(sigma^2), main = "")
abline(v = var.phi, col = "red", lwd = 2)


# 7.4.3. Temporal covariates
# Define parameter values
n.occasions <- 20                  # Number of capture occasions
marked <- rep(15, n.occasions-1)   # Annual number of newly marked individuals
mean.phi <- 0.65
p <- rep(0.4, n.occasions-1)
beta <- -0.3                       # Slope of survival-winter relationship	
r.var <- 0.2                       # Residual temporal variance

# Draw annual survival probabilities
winter <- rnorm(n.occasions-1, 0, 1^0.5)
logit.phi <- qlogis(mean.phi) + beta*winter + rnorm(n.occasions-1, 0, r.var^0.5)
phi <- plogis(logit.phi)

# Define matrices with survival and recapture probabilities
PHI <- matrix(phi, ncol = n.occasions-1, nrow = sum(marked), byrow = TRUE)
P <- matrix(p, ncol = n.occasions-1, nrow = sum(marked))

# Simulate capture-histories
CH <- simul.cjs(PHI, P, marked)

# Create vector with occasion of marking
get.first <- function(x) min(which(x!=0))
f <- apply(CH, 1, get.first)


# Specify model in BUGS language
sink("cjs-cov-raneff.jags")
cat("
model {

# Priors and constraints
for (i in 1:nind){
   for (t in f[i]:(n.occasions-1)){
      logit(phi[i,t]) <- mu + beta*x[t] + epsilon[t]
      p[i,t] <- mean.p
      } #t
   } #i
for (t in 1:(n.occasions-1)){
   epsilon[t] ~ dnorm(0, tau)
   phi.est[t] <- 1 / (1+exp(-mu-beta*x[t]-epsilon[t])) # Yearly survival
   }
mu ~ dnorm(0, 0.001)                     # Prior for logit of mean survival
mean.phi <- 1 / (1+exp(-mu))             # Logit transformation
beta ~ dnorm(0, 0.001)I(-10, 10)         # Prior for slope parameter
sigma ~ dunif(0, 10)                     # Prior on standard deviation
tau <- pow(sigma, -2)
sigma2 <- pow(sigma, 2)                  # Residual temporal variance
mean.p ~ dunif(0, 1)                     # Prior for mean recapture

# Likelihood 
for (i in 1:nind){
   # Define latent state at first capture
   z[i,f[i]] <- 1
   for (t in (f[i]+1):n.occasions){
      # State process
      z[i,t] ~ dbern(mu1[i,t])
      mu1[i,t] <- phi[i,t-1] * z[i,t-1]
      # Observation process
      y[i,t] ~ dbern(mu2[i,t])
      mu2[i,t] <- p[i,t-1] * z[i,t]
      } #t
   } #i
}
",fill = TRUE)
sink()

# Bundle data
jags.data <- list(y = CH, f = f, nind = dim(CH)[1], n.occasions = dim(CH)[2], z = known.state.cjs(CH), x = winter)

# Initial values
inits <- function(){list(z = cjs.init.z(CH, f), mu = rnorm(1), sigma = runif(1, 0, 5), beta = runif(1, -5, 5), mean.p = runif(1, 0, 1))}  

# Parameters monitored
parameters <- c("mean.phi", "mean.p", "phi.est", "sigma2", "beta")

# MCMC settings
ni <- 20000
nt <- 6
nb <- 10000
nc <- 3

# Call JAGS from R (BRT 12 min)
cjs.cov <- jags(jags.data, inits, parameters, "cjs-cov-raneff.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

# Summarize posteriors
print(cjs.cov, digits = 3)

# Produce graph
par(mfrow = c(1, 2), las = 1)
hist(cjs.cov$BUGSoutput$sims.list$beta, nclass = 25, col = "gray", main = "", xlab = expression(beta), ylab = "Frequency")
abline(v = -0.3, col = "red", lwd = 2)
hist(cjs.cov$BUGSoutput$sims.list$sigma2, nclass = 50, col = "gray", main = "", xlab = expression(sigma^2), ylab = "Frequency", xlim=c(0, 3))
abline(v = 0.2, col = "red", lwd = 2)


# 7.5. Models with individual variation
# 7.5.1. Fixed group effects
# Define parameter values
n.occasions <- 12                  # Number of capture occasions
marked <- rep(30, n.occasions-1)   # Annual number of newly marked individuals
phi.f <- rep(0.65, n.occasions-1)  # Survival of females
p.f <- rep(0.6, n.occasions-1)     # Recapture of females
phi.m <- rep(0.8, n.occasions-1)   # Survival of males
p.m <- rep(0.3, n.occasions-1)     # Reacpture of males

# Define matrices with survival and recapture probabilities
PHI.F <- matrix(phi.f, ncol = n.occasions-1, nrow = sum(marked))
P.F <- matrix(p.f, ncol = n.occasions-1, nrow = sum(marked))
PHI.M <- matrix(phi.m, ncol = n.occasions-1, nrow = sum(marked))
P.M <- matrix(p.m, ncol = n.occasions-1, nrow = sum(marked))

# Simulate capture-histories
CH.F <- simul.cjs(PHI.F, P.F, marked)
CH.M <- simul.cjs(PHI.M, P.M, marked)

# Merge capture-histories by row
CH <- rbind(CH.F, CH.M)

# Create group variable
group <- c(rep(1, dim(CH.F)[1]), rep(2, dim(CH.M)[1]))

# Create vector with occasion of marking
get.first <- function(x) min(which(x!=0))
f <- apply(CH, 1, get.first)

# Specify model in JAGS language
sink("cjs-group.jags")
cat("
model {

# Priors and constraints
for (i in 1:nind){
   for (t in f[i]:(n.occasions-1)){
      phi[i,t] <- phi.g[group[i]]
      p[i,t] <- p.g[group[i]]
      } #t
   } #i
for (u in 1:g){
   phi.g[u] ~ dunif(0, 1)              # Priors for group-specific survival
   p.g[u] ~ dunif(0, 1)                # Priors for group-specific recapture
   }

# Likelihood 
for (i in 1:nind){
   # Define latent state at first capture
   z[i,f[i]] <- 1
   for (t in (f[i]+1):n.occasions){
      # State process
      z[i,t] ~ dbern(mu1[i,t])
      mu1[i,t] <- phi[i,t-1] * z[i,t-1]
      # Observation process
      y[i,t] ~ dbern(mu2[i,t])
      mu2[i,t] <- p[i,t-1] * z[i,t]
      } #t
   } #i
}
",fill = TRUE)
sink()

# Bundle data
jags.data <- list(y = CH, f = f, nind = dim(CH)[1], n.occasions = dim(CH)[2], z = known.state.cjs(CH), g = length(unique(group)), group = group)

# Initial values
inits <- function(){list(z = cjs.init.z(CH, f), phi.g = runif(length(unique(group)), 0, 1), p.g = runif(length(unique(group)), 0, 1))}  

# Parameters monitored
parameters <- c("phi.g", "p.g")

# MCMC settings
ni <- 5000
nt <- 3
nb <- 2000
nc <- 3

# Call JAGS from R (BRT 2 min)
cjs.group <- jags(jags.data, inits, parameters, "cjs-group.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

# Summarize posteriors
print(cjs.group, digits = 3)


# 7.5.2. Random group effects
# Priors and constraints
for (i in 1:nind){
   for (t in f[i]:(n.occasions-1)){
      logit(phi[i,t]) <- beta[group[i]]
      p[i,t] <- mean.p
      } #t
   } #i
for (u in 1:g){
   beta[u] ~ dnorm(mean.beta, tau)
   phi.g[u] <- 1 / (1+exp(-beta[u])) # Back-transformed group-specific survival
   }
mean.beta ~ dnorm(0, 0.001)         # Prior for logit of mean survival
mean.phi <- 1 / (1+exp(-mean.beta)) # Back-transformed mean survival
sigma ~ dunif(0, 10)                # Prior for sd of logit of survival variability
tau <- pow(sigma, -2)
mean.p ~ dunif(0, 1)                # Prior for mean recapture


# 7.5.3. Individual random effects
# Define parameter values
n.occasions <- 20                 # Number of capture occasions
marked <- rep(30, n.occasions-1)  # Annual number of newly marked individuals
mean.phi <- 0.65
p <- rep(0.4, n.occasions-1)
v.ind <- 0.5

# Draw annual survival probabilities
logit.phi <- rnorm(sum(marked), qlogis(mean.phi), v.ind^0.5)
phi <- plogis(logit.phi)

# Define matrices with survival and recapture probabilities
PHI <- matrix(phi, ncol = n.occasions-1, nrow = sum(marked), byrow = FALSE)
P <- matrix(p, ncol = n.occasions-1, nrow = sum(marked))

# Simulate capture-histories
CH <- simul.cjs(PHI, P, marked)

# Create vector with occasion of marking
get.first <- function(x) min(which(x!=0))
f <- apply(CH, 1, get.first)

# Specify model in BUGS language
sink("cjs-ind-raneff.jags")
cat("
model {

# Priors and constraints
for (i in 1:nind){
   for (t in f[i]:(n.occasions-1)){
      logit(phi[i,t]) <- mu + epsilon[i]
      p[i,t] <- mean.p
      } #t
   } #i
for (i in 1:nind){
   epsilon[i] ~ dnorm(0, tau)
   }
mean.phi ~ dunif(0, 1)                   # Prior for mean survival
mu <- log(mean.phi / (1-mean.phi))       # Logit transformation
sigma ~ dunif(0, 5)                      # Prior for standard deviation
tau <- pow(sigma, -2)
sigma2 <- pow(sigma, 2)
mean.p ~ dunif(0, 1)                     # Prior for mean recapture 

# Likelihood 
for (i in 1:nind){
   # Define latent state at first capture
   z[i,f[i]] <- 1
   for (t in (f[i]+1):n.occasions){
      # State process
      z[i,t] ~ dbern(mu1[i,t])
      mu1[i,t] <- phi[i,t-1] * z[i,t-1]
      # Observation process
      y[i,t] ~ dbern(mu2[i,t])
      mu2[i,t] <- p[i,t-1] * z[i,t]
      } #t
   } #i
}
",fill = TRUE)
sink()

# Bundle data
jags.data <- list(y = CH, f = f, nind = dim(CH)[1], n.occasions = dim(CH)[2], z = known.state.cjs(CH))

# Initial values 
inits <- function(){list(z = cjs.init.z(CH, f), mean.phi = runif(1, 0, 1), mean.p = runif(1, 0, 1), sigma = runif(1, 0, 2))}  

# Parameters monitored
parameters <- c("mean.phi", "mean.p", "sigma2")

# MCMC settings
ni <- 50000
nt <- 6
nb <- 20000
nc <- 3

# Call JAGS from R (BRT 73 min)
cjs.ind <- jags(jags.data, inits, parameters, "cjs-ind-raneff.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

# Summarize posteriors
print(cjs.ind, digits = 3)

# Produce graph
par(mfrow = c(1, 2), las = 1)
hist(cjs.ind$BUGSoutput$sims.list$mean.phi, nclass = 25, col = "gray", main = "", xlab = expression(bar(phi)), ylab = "Frequency")
abline(v = mean.phi, col = "red", lwd = 2)
hist(cjs.ind$BUGSoutput$sims.list$sigma2, nclass = 15, col = "gray", main = "", xlab = expression(sigma^2), ylab = "Frequency", xlim = c(0, 3))
abline(v = v.ind, col = "red", lwd = 2)

# Priors and constraints
for (i in 1:nind){
   for (t in f[i]:(n.occasions-1)){
      logit(phi[i,t]) <- mu + beta*x[i] + epsilon[i]
      p[i,t] <- mean.p
      } #t
   } #i
for (i in 1:nind){
   epsilon[i] ~ dnorm(0, tau)
   }
mean.phi ~ dunif(0, 1)                   # Prior for mean survival 
mu <- log(mean.phi / (1-mean.phi))       # Logit transformation
beta ~ dnorm(0, 0.001)                   # Prior for covariate slope
sigma ~ dunif(0, 5)                      # Prior for standard deviation
tau <- pow(sigma, -2)
sigma2 <- pow(sigma, 2)
mean.p ~ dunif(0, 1)                     # Prior for mean recapture 


# 7.6. Models with time and group effects
# 7.6.1. Fixed group and time effects
# Define parameter values
n.occasions <- 12                  # Number of capture occasions
marked <- rep(50, n.occasions-1)   # Annual number of newly marked individuals
phi.f <- c(0.6, 0.5, 0.55, 0.6, 0.5, 0.4, 0.6, 0.5, 0.55, 0.6, 0.7)
p.f <- rep(0.6, n.occasions-1)
diff <- 0.5     # Difference between male and female survival on logit scale
phi.m <- plogis(qlogis(phi.f) + diff)
p.m <- rep(0.3, n.occasions-1)

# Define matrices with survival and recapture probabilities
PHI.F <- matrix(rep(phi.f, sum(marked)), ncol = n.occasions-1, nrow = sum(marked), byrow = TRUE)
P.F <- matrix(rep(p.f, sum(marked)), ncol = n.occasions-1, nrow = sum(marked), byrow = TRUE)
PHI.M <- matrix(rep(phi.m, sum(marked)), ncol = n.occasions-1, nrow = sum(marked), byrow = TRUE)
P.M <- matrix(rep(p.m, sum(marked)), ncol = n.occasions-1, nrow = sum(marked), byrow = TRUE)

# Simulate capture-histories
CH.F <- simul.cjs(PHI.F, P.F, marked)
CH.M <- simul.cjs(PHI.M, P.M, marked)

# Merge capture-histories
CH <- rbind(CH.F, CH.M)

# Create group variable
group <- c(rep(1, dim(CH.F)[1]), rep(2, dim(CH.M)[1]))

# Create vector with occasion of marking
get.first <- function(x) min(which(x!=0))
f <- apply(CH, 1, get.first)

# Specify model in BUGS language
sink("cjs-add.jags")
cat("
model {

# Priors and constraints
for (i in 1:nind){
   for (t in f[i]:(n.occasions-1)){
      logit(phi[i,t]) <- beta[group[i]] + gamma[t]
      p[i,t] <- p.g[group[i]]
      } #t
   } #i
# for survival parameters
for (t in 1:(n.occasions-1)){
   gamma[t] ~ dnorm(0, 0.01)I(-10, 10)          # Priors for time effects
   phi.g1[t] <- 1 / (1+exp(-gamma[t]))          # Back-transformed survival of males
   phi.g2[t] <- 1 / (1+exp(-gamma[t]-beta[2]))  # Back-transformed survival of females 
   }
beta[1] <- 0                            # Corner constraint
beta[2] ~ dnorm(0, 0.01)I(-10, 10)      # Prior for difference in male and female survival
# for recapture parameters
for (u in 1:g){
   p.g[u] ~ dunif(0, 1)                 # Priors for group-spec. recapture
   }

# Likelihood 
for (i in 1:nind){
   # Define latent state at first capture
   z[i,f[i]] <- 1
   for (t in (f[i]+1):n.occasions){
      # State process
      z[i,t] ~ dbern(mu1[i,t])
      mu1[i,t] <- phi[i,t-1] * z[i,t-1]
      # Observation process
      y[i,t] ~ dbern(mu2[i,t])
      mu2[i,t] <- p[i,t-1] * z[i,t]
      } #t
   } #i
}
",fill = TRUE)
sink()

# Bundle data
jags.data <- list(y = CH, f = f, nind = dim(CH)[1], n.occasions = dim(CH)[2], z = known.state.cjs(CH), g = length(unique(group)), group = group)

# Initial values
inits <- function(){list(z = cjs.init.z(CH, f), gamma = rnorm(n.occasions-1), beta = c(NA, rnorm(1)), p.g = runif(length(unique(group)), 0, 1))}  

# Parameters monitored
parameters <- c("phi.g1", "phi.g2", "p.g", "beta")

# MCMC settings
ni <- 5000
nt <- 3
nb <- 2000
nc <- 3

# Call JAGS from R (BRT 7 min)
cjs.add <- jags(jags.data, inits, parameters, "cjs-add.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

# Summarize posteriors
print(cjs.add, digits = 3)   

# Figure of male and female survival
lower.f <- upper.f <- lower.m <- upper.m <- numeric()
for (t in 1:(n.occasions-1)){
   lower.f[t] <- quantile(cjs.add$BUGSoutput$sims.list$phi.g1[,t], 0.025)
   upper.f[t] <- quantile(cjs.add$BUGSoutput$sims.list$phi.g1[,t], 0.975)
   lower.m[t] <- quantile(cjs.add$BUGSoutput$sims.list$phi.g2[,t], 0.025)
   upper.m[t] <- quantile(cjs.add$BUGSoutput$sims.list$phi.g2[,t], 0.975)
   }
plot(x=(1:(n.occasions-1))-0.1, y = cjs.add$BUGSoutput$mean$phi.g1, type = "b", pch = 16, ylim = c(0.2, 1), ylab = "Survival probability", xlab = "Year", bty = "n", cex = 1.5, axes = FALSE)
axis(1, at = 1:11, labels = rep(NA,11), tcl = -0.25)
axis(1, at = seq(2,10,2), labels = c("2","4","6","8","10"))
axis(2, at = seq(0.2, 1, 0.1), labels = c("0.2", NA, "0.4", NA, "0.6", NA, "0.8", NA, "1.0"), las = 1)
segments((1:(n.occasions-1))-0.1, lower.f, (1:(n.occasions-1))-0.1, upper.f)
points(x = (1:(n.occasions-1))+0.1, y = cjs.add$BUGSoutput$mean$phi.g2, type = "b", pch = 1, lty = 2, cex = 1.5)
segments((1:(n.occasions-1))+0.1, lower.m, (1:(n.occasions-1))+0.1, upper.m)


# Priors and constraints
for (i in 1:nind){
   for (t in f[i]:(n.occasions-1)){
      phi[i,t] <- eta.phi[group[i],t]
      p[i,t] <- p.g[group[i]]
      } #t
   } #i
# for survival parameters
for (u in 1:g){
   for (t in 1:(n.occasions-1)){
      eta.phi[u,t] ~ dunif(0, 1)     # Prior for time and group-spec. survival
      } #t
   } #g
# for recapture parameters
for (u in 1:g){
   p.g[u] ~ dunif(0, 1)              # Priors for group-spec. recapture
   }


# 7.6.2. Fixed group and random time effects
# Priors and constraints
for (i in 1:nind){
   for (t in f[i]:(n.occasions-1)){
      logit(phi[i,t]) <- eta.phi[group[i],t]
      p[i,t] <- p.g[group[i]]
      } #t
   } #i
# for survival parameters
for (u in 1:g){
   for (t in 1:(n.occasions-1)){
      eta.phi[u,t] <- mu.phi[u] + epsilon[u,t]
      epsilon[u,t] ~ dnorm(0, tau[u])
      } #t
      mean.phi[u] ~ dunif(0, 1)      # Priors on mean group-spec. survival     
      mu.phi[u] <- log(mean.phi[u] / (1-mean.phi[u]))
      sigma[u] ~ dunif(0, 10)        # Priors for group-spec. sd 
      tau[u] <- pow(sigma[u], -2)
      sigma2[u] <- pow(sigma[u], 2)
   } #g
# for recapture parameters
for (u in 1:g){
   p.g[u] ~ dunif(0,1)               # Priors for group-spec. recapture
   }

# Specify model in BUGS language
sink("cjs-temp-corr.jags")
cat("
model {

# Priors and constraints
for (i in 1:nind){
   for (t in f[i]:(n.occasions-1)){
      logit(phi[i,t]) <- eta.phi[t,group[i]]
      p[i,t] <- p.g[group[i]]
      } #t
   } #i
# for survival parameters
for (t in 1:(n.occasions-1)){
   eta.phi[t,1:g] ~ dmnorm(mu.phi[], Omega[,])
   } #t
for (u in 1:g){      
   mean.phi[u] ~ dunif(0, 1)    # Priors on mean group-spec. survival
   mu.phi[u] <- log(mean.phi[u] / (1-mean.phi[u]))
   } #g
Omega[1:g, 1:g] ~ dwish(R[,], df)  # Priors for variance-covariance matrix
Sigma[1:g, 1:g] <- inverse(Omega[,])

# for recapture parameters
for (u in 1:g){
   p.g[u] ~ dunif(0, 1)            # Priors for group-spec. recapture
   }

# Likelihood 
for (i in 1:nind){
   # Define latent state at first capture
   z[i,f[i]] <- 1
   for (t in (f[i]+1):n.occasions){
      # State process
      z[i,t] ~ dbern(mu1[i,t])
      mu1[i,t] <- phi[i,t-1] * z[i,t-1]
      # Observation process
      y[i,t] ~ dbern(mu2[i,t])
      mu2[i,t] <- p[i,t-1] * z[i,t]
      } #t
   } #i
}
",fill = TRUE)
sink()

# Bundle data
jags.data <- list(y = CH, f = f, nind = dim(CH)[1], n.occasions = dim(CH)[2], z = known.state.cjs(CH), g = length(unique(group)), group = group, R = matrix(c(5, 0, 0, 1), ncol = 2), df = 3)

# Initial values
inits <- function(){list(z = cjs.init.z(CH, f), p.g = runif(length(unique(group)), 0, 1), Omega = matrix(c(1, 0, 0, 1), ncol = 2))}  

# Parameters monitored
parameters <- c("eta.phi", "p.g", "Sigma", "mean.phi")

# MCMC settings
ni <- 5000
nt <- 3
nb <- 2000
nc <- 3

# Call JAGS from R (BRT 5 min)
cjs.corr <- jags(jags.data, inits, parameters, "cjs-temp-corr.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

# Summarize posteriors
print(cjs.corr, digits = 3)   


# 7.7. Models with age effects
# Define parameter values
n.occasions <- 10                   # Number of capture occasions
marked.j <- rep(200, n.occasions-1) # Annual number of newly marked juveniles
marked.a <- rep(30, n.occasions-1)  # Annual number of newly marked adults
phi.juv <- 0.3                      # Juvenile annual survival
phi.ad <- 0.65                      # Adult annual survival
p <- rep(0.5, n.occasions-1)        # Recapture
phi.j <- c(phi.juv, rep(phi.ad, n.occasions-2))
phi.a <- rep(phi.ad, n.occasions-1)

# Define matrices with survival and recapture probabilities
PHI.J <- matrix(0, ncol = n.occasions-1, nrow = sum(marked.j))
for (i in 1:length(marked.j)){
   PHI.J[(sum(marked.j[1:i])-marked.j[i]+1):sum(marked.j[1:i]),i:(n.occasions-1)] <- matrix(rep(phi.j[1:(n.occasions-i)],marked.j[i]), ncol = n.occasions-i, byrow = TRUE)
   }
P.J <- matrix(rep(p, sum(marked.j)), ncol = n.occasions-1, nrow = sum(marked.j), byrow = TRUE)
PHI.A <- matrix(rep(phi.a, sum(marked.a)), ncol = n.occasions-1, nrow = sum(marked.a), byrow = TRUE)
P.A <- matrix(rep(p, sum(marked.a)), ncol = n.occasions-1, nrow = sum(marked.a), byrow = TRUE)

# Apply simulation function
CH.J <- simul.cjs(PHI.J, P.J, marked.j)
CH.A <- simul.cjs(PHI.A, P.A, marked.a) 

# Create vector with occasion of marking
get.first <- function(x) min(which(x!=0))
f.j <- apply(CH.J, 1, get.first)
f.a <- apply(CH.A, 1, get.first)

# Create matrices X indicating age classes
x.j <- matrix(NA, ncol = dim(CH.J)[2]-1, nrow = dim(CH.J)[1])
x.a <- matrix(NA, ncol = dim(CH.A)[2]-1, nrow = dim(CH.A)[1])
for (i in 1:dim(CH.J)[1]){
   for (t in f.j[i]:(dim(CH.J)[2]-1)){
      x.j[i,t] <- 2
      x.j[i,f.j[i]] <- 1   
      } #t
   } #i
for (i in 1:dim(CH.A)[1]){
   for (t in f.a[i]:(dim(CH.A)[2]-1)){
      x.a[i,t] <- 2
      } #t
   } #i

CH <- rbind(CH.J, CH.A)
f <- c(f.j, f.a)
x <- rbind(x.j, x.a)

# Specify model in BUGS language
sink("cjs-age.jags")
cat("
model {
# Priors and constraints
for (i in 1:nind){
   for (t in f[i]:(n.occasions-1)){
      phi[i,t] <- beta[x[i,t]]
      p[i,t] <- mean.p
      } #t
   } #i
for (u in 1:2){
   beta[u] ~ dunif(0, 1)              # Priors for age-specific survival
   }
mean.p ~ dunif(0, 1)                  # Prior for mean recapture
# Likelihood 
for (i in 1:nind){
   # Define latent state at first capture
   z[i,f[i]] <- 1
   for (t in (f[i]+1):n.occasions){
      # State process
      z[i,t] ~ dbern(mu1[i,t])
      mu1[i,t] <- phi[i,t-1] * z[i,t-1]
      # Observation process
      y[i,t] ~ dbern(mu2[i,t])
      mu2[i,t] <- p[i,t-1] * z[i,t]
      } #t
   } #i
}
",fill = TRUE)
sink()

# Bundle data
jags.data <- list(y = CH, f = f, nind = dim(CH)[1], n.occasions = dim(CH)[2], z = known.state.cjs(CH), x = x)

# Initial values
inits <- function(){list(z = cjs.init.z(CH, f), beta = runif(2, 0, 1), mean.p = runif(1, 0, 1))}  

# Parameters monitored
parameters <- c("beta", "mean.p")

# MCMC settings
ni <- 2000
nt <- 3
nb <- 1000
nc <- 3

# Call JAGS from R (BRT 3 min)
cjs.age <- jags(jags.data, inits, parameters, "cjs-age.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

print(cjs.age, digits = 3)

# Create matrix X indicating age classes
x <- matrix(NA, ncol = dim(CH)[2]-1, nrow = dim(CH)[1])
for (i in 1:dim(CH)[1]){
   for (t in f[i]:(dim(CH)[2]-1)){
      x[i,t] <- t-f[i]+1
      } #t 
   } #i

# Priors and constraints
for (i in 1:nind){
   for (t in f[i]:(n.occasions-1)){
      logit(phi[i,t]) <- mu + beta*x[i,t]
      p[i,t] <- mean.p
      } #t
   } #i
mu ~ dnorm(0, 0.01)             # Prior for mean of logit survival
beta ~ dnorm(0, 0.01)           # Prior for slope parameter
for (i in 1:(n.occasions-1)){
   phi.age[i] <- 1 / (1+exp(-mu -beta*i))   # Logit back-transformation 
   }
mean.p ~ dunif(0, 1)                # Prior for mean recapture


# 7.8. Immediate trap response in recapture probability
# Import data
CH <- as.matrix(read.table(file = "trap.txt", sep = " "))

# Compute vector with occasion of first capture
get.first <- function(x) min(which(x!=0))
f <- apply(CH, 1, get.first)

# Create matrix m indicating when an individual was captured
m <- CH[,1:(dim(CH)[2]-1)]
u <- which(m==0)
m[u] <- 2

# Specify model in BUGS language
sink("cjs-trap.jags")
cat("
model {

# Priors and constraints
for (i in 1:nind){
   for (t in f[i]:(n.occasions-1)){
      phi[i,t] <- mean.phi
      p[i,t] <- beta[m[i,t]]
      } #t
   } #i
mean.phi ~ dunif(0, 1)                # Prior for mean survival
for (u in 1:2){
   beta[u] ~ dunif(0, 1)              # Priors for recapture
   }

# Likelihood components
for (i in 1:nind){
   # Define latent state at first capture
   z[i,f[i]] <- 1
   for (t in (f[i]+1):n.occasions){
      # State process
      z[i,t] ~ dbern(mu1[i,t])
      mu1[i,t] <- phi[i,t-1] * z[i,t-1]
      # Observation process
      y[i,t] ~ dbern(mu2[i,t])
      mu2[i,t] <- p[i,t-1] * z[i,t]
      } #t
   } #i
}
",fill = TRUE)
sink()

# Bundle data
jags.data <- list(y = CH, f = f, nind = dim(CH)[1], n.occasions = dim(CH)[2], z = known.state.cjs(CH), m = m)

# Initial values
inits <- function(){list(z = cjs.init.z(CH, f), mean.phi = runif(1, 0, 1), beta = runif(2, 0, 1))}

# Parameters monitored
parameters <- c("mean.phi", "beta")

# MCMC settings
ni <- 20000
nt <- 3
nb <- 10000
nc <- 3

# Call JAGS from R (BRT 1 min)
cjs.trap <- jags(jags.data, inits, parameters, "cjs-trap.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

print(cjs.trap, digits = 3)


# 7.9. Parameter identifiability
# Define parameter values
n.occasions <- 12                  # Number of capture occasions
marked <- rep(30, n.occasions-1)   # Annual number of newly marked individuals
phi <- c(0.6, 0.5, 0.55, 0.6, 0.5, 0.4, 0.6, 0.5, 0.55, 0.6, 0.7)
p <- c(0.4, 0.65, 0.4, 0.45, 0.55, 0.68, 0.66, 0.28, 0.55, 0.45, 0.35)

# Define matrices with survival and recapture probabilities 
PHI <- matrix(phi, ncol = n.occasions-1, nrow = sum(marked), byrow = TRUE)
P <- matrix(p, ncol = n.occasions-1, nrow = sum(marked), byrow = TRUE)

# Simulate capture-histories
CH <- simul.cjs(PHI, P, marked)

# Create vector with occasion of marking
get.first <- function(x) min(which(x!=0))
f <- apply(CH, 1, get.first)

# Specify model in BUGS language
sink("cjs-t-t.jags")
cat("
model {

# Priors and constraints
for (i in 1:nind){
   for (t in f[i]:(n.occasions-1)){
      phi[i,t] <- phi.t[t]
      p[i,t] <- p.t[t]
      } #t
   } #i
for (t in 1:(n.occasions-1)){
   phi.t[t] ~ dunif(0, 1)          # Priors for time-spec. survival
   p.t[t] ~ dunif(0, 1)            # Priors for time-spec. recapture
   }

# Likelihood 
for (i in 1:nind){
   # Define latent state at first capture
   z[i,f[i]] <- 1
   for (t in (f[i]+1):n.occasions){
      # State process
      z[i,t] ~ dbern(mu1[i,t])
      mu1[i,t] <- phi[i,t-1] * z[i,t-1]
      # Observation process
      y[i,t] ~ dbern(mu2[i,t])
      mu2[i,t] <- p[i,t-1] * z[i,t]
      } #t
   } #i
}
",fill = TRUE)
sink()

# Bundle data
jags.data <- list(y = CH, f = f, nind = dim(CH)[1], n.occasions = dim(CH)[2], z = known.state.cjs(CH))

# Initial values
inits <- function(){list(z = cjs.init.z(CH, f), phi.t = runif((dim(CH)[2]-1), 0, 1), p.t = runif((dim(CH)[2]-1), 0, 1))}

# Parameters monitored
parameters <- c("phi.t", "p.t")

# MCMC settings
ni <- 25000
nt <- 3
nb <- 20000
nc <- 3

# Call JAGS from R (BRT 7 min)
cjs.t.t <- jags(jags.data, inits, parameters, "cjs-t-t.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

# Plot posterior distributions of some phi and p
par(mfrow = c(2, 2), cex = 1.2, las = 1, mar=c(5, 4, 2, 1))
plot(density(cjs.t.t$BUGSoutput$sims.list$phi.t[,6]), xlim = c(0, 1), ylim = c(0, 5), main = "", xlab = expression(phi[6]), ylab = "Density", frame = FALSE, lwd = 2)
abline(h = 1, lty = 2, lwd = 2)
par(mar=c(5, 3, 2, 2))
plot(density(cjs.t.t$BUGSoutput$sims.list$phi.t[,11]), xlim = c(0, 1), ylim = c(0, 5), main = "", xlab = expression(phi[11]), ylab ="", frame = FALSE, lwd = 2)
abline(h = 1, lty = 2, lwd = 2)
par(mar=c(5, 4, 2, 1))
plot(density(cjs.t.t$BUGSoutput$sims.list$p.t[,6]), xlim = c(0, 1), ylim = c(0, 5), main = "", xlab = expression(p[6]), ylab = "Density", frame = FALSE, lwd = 2)
abline(h = 1, lty = 2, lwd = 2)
par(mar=c(5, 3, 2, 2))
plot(density(cjs.t.t$BUGSoutput$sims.list$p.t[,11]), xlim = c(0, 1), ylim = c(0, 5), main = "", xlab = expression(p[11]), ylab ="", frame = FALSE, lwd = 2)
abline(h = 1, lty = 2, lwd = 2) 


# 7.10. Fitting the CJS to data in the m-array format: the multinomial likelihood
# 7.10.1. Introduction
# Function to create a m-array based on capture-histories (CH)
marray <- function(CH){
   nind <- dim(CH)[1]
   n.occasions <- dim(CH)[2]
   m.array <- matrix(data = 0, ncol = n.occasions+1, nrow = n.occasions)
   # Calculate the number of released individuals at each time period
   for (t in 1:n.occasions){
      m.array[t,1] <- sum(CH[,t])
      }
   for (i in 1:nind){
      pos <- which(CH[i,]!=0)
      g <- length(pos)
      for (z in 1:(g-1)){
         m.array[pos[z],pos[z+1]] <- m.array[pos[z],pos[z+1]] + 1
         } #z
      } #i
   # Calculate the number of individuals that is never recaptured
   for (t in 1:n.occasions){
      m.array[t,n.occasions+1] <- m.array[t,1] - sum(m.array[t,2:n.occasions])
      }
   out <- m.array[1:(n.occasions-1),2:(n.occasions+1)]
   return(out)
   }


# 7.10.2. Time-dependent models
# Specify model in BUGS language
sink("cjs-mnl.jags")
cat("
model {
# Priors and constraints
for (t in 1:(n.occasions-1)){
   phi[t] ~ dunif(0, 1)         # Priors for survival
   p[t] ~ dunif(0, 1)           # Priors for recapture
   }
# Define the multinomial likelihood
for (t in 1:(n.occasions-1)){
   marr[t,1:n.occasions] ~ dmulti(pr[t, ], rel[t])
   }
# Define the cell probabilities of the m-array
# Main diagonal
for (t in 1:(n.occasions-1)){
   q[t] <- 1-p[t]                # Probability of non-recapture
   pr[t,t] <- phi[t]*p[t]
   # Above main diagonal
   for (j in (t+1):(n.occasions-1)){
      pr[t,j] <- prod(phi[t:j])*prod(q[t:(j-1)])*p[j]
      } #j
   # Below main diagonal
   for (j in 1:(t-1)){
      pr[t,j] <- 0
      } #j
   } #t
# Last column: probability of non-recapture
for (t in 1:(n.occasions-1)){
   pr[t,n.occasions] <- 1-sum(pr[t,1:(n.occasions-1)])
   } #t

# Assess model fit using Freeman-Tukey statistic
# Compute fit statistics for observed data
for (t in 1:(n.occasions-1)){
   for (j in 1:n.occasions){
      expmarr[t,j] <- rel[t]*pr[t,j]
      E.org[t,j] <- pow((pow(marr[t,j], 0.5)-pow(expmarr[t,j], 0.5)), 2)
      } #j
   } #t
# Generate replicate data and compute fit stats from them
for (t in 1:(n.occasions-1)){
   marr.new[t,1:n.occasions] ~ dmulti(pr[t, ], rel[t])
   for (j in 1:n.occasions){
      E.new[t,j] <- pow((pow(marr.new[t,j], 0.5)-pow(expmarr[t,j], 0.5)), 2)
      } #j
   } #t
fit <- sum(E.org[,])
fit.new <- sum(E.new[,])
}
",fill = TRUE)
sink()

# Create the m-array from the capture-histories
marr <- marray(CH)

# Bundle data
jags.data <- list(marr = marr, n.occasions = dim(marr)[2], rel = rowSums(marr))

# Initial values
inits <- function(){list(phi = runif(dim(marr)[2]-1, 0, 1), p = runif(dim(marr)[2]-1, 0, 1))}  

# Parameters monitored
parameters <- c("phi", "p", "fit", "fit.new")

# MCMC settings
ni <- 10000
nt <- 3
nb <- 5000
nc <- 3

# Call JAGS from R (BRT 1 min)
cjs <- jags(jags.data, inits, parameters, "cjs-mnl.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

print(cjs, digits = 3) 

# Evaluation of fit
plot(cjs$BUGSoutput$sims.list$fit, cjs$BUGSoutput$sims.list$fit.new, xlab = "Discrepancy actual data", ylab = "Discrepancy replicate data", las = 1,  
ylim = c(5, 25), xlim = c(5, 25), bty ="n") 
abline(0, 1, col = "black", lwd = 2)
mean(cjs$BUGSoutput$sims.list$fit.new > cjs$BUGSoutput$sims.list$fit)


# 7.10.3. Age-dependent models
# Define parameter values
n.occasions <- 12                    # Number of capture occasions
marked.j <- rep(200, n.occasions-1)  # Annual number of newly marked juveniles
marked.a <- rep(30, n.occasions-1)   # Annual number of newly marked adults
phi.juv <- 0.3                       # Juvenile annual survival
phi.ad <- 0.65                       # Adult annual survival
p <- rep(0.5, n.occasions-1)         # Recapture
phi.j <- c(phi.juv, rep(phi.ad,n.occasions-2))
phi.a <- rep(phi.ad, n.occasions-1)

# Define matrices with survival and recapture probabilities
PHI.J <- matrix(0, ncol = n.occasions-1, nrow = sum(marked.j))
for (i in 1:(length(marked.j)-1)){
   PHI.J[(sum(marked.j[1:i])-marked.j[i]+1):sum(marked.j[1:i]),i:(n.occasions-1)] <- matrix(rep(phi.j[1:(n.occasions-i)],marked.j[i]), ncol = n.occasions-i, byrow = TRUE)
   }
P.J <- matrix(rep(p, n.occasions*sum(marked.j)), ncol = n.occasions-1, nrow = sum(marked.j), byrow = TRUE)
PHI.A <- matrix(rep(phi.a, sum(marked.a)), ncol = n.occasions-1, nrow = sum(marked.a), byrow = TRUE)
P.A <- matrix(rep(p, sum(marked.a)), ncol = n.occasions-1, nrow = sum(marked.a), byrow = TRUE)

# Apply simulation function
CH.J <- simul.cjs(PHI.J, P.J, marked.j)
CH.A <- simul.cjs(PHI.A, P.A, marked.a) 

cap <- apply(CH.J, 1, sum)
ind <- which(cap >= 2)
CH.J.R <- CH.J[ind,]    # Juvenile CH recaptured at least once
CH.J.N <- CH.J[-ind,]   # Juvenile CH never recaptured
# Remove first capture
first <- numeric()
for (i in 1:dim(CH.J.R)[1]){
   first[i] <- min(which(CH.J.R[i,]==1))
   }
CH.J.R1 <- CH.J.R
for (i in 1:dim(CH.J.R)[1]){
   CH.J.R1[i,first[i]] <- 0
   }
# Add grown-up juveniles to adults and create m-array
CH.A.m <- rbind(CH.A, CH.J.R1)
CH.A.marray <- marray(CH.A.m)
# Create CH matrix for juveniles, ignoring subsequent recaptures
second <- numeric()
for (i in 1:dim(CH.J.R1)[1]){
   second[i] <- min(which(CH.J.R1[i,]==1))
   }
CH.J.R2 <- matrix(0, nrow = dim(CH.J.R)[1], ncol = dim(CH.J.R)[2])
for (i in 1:dim(CH.J.R)[1]){
   CH.J.R2[i,first[i]] <- 1
   CH.J.R2[i,second[i]] <- 1
   }
# Create m-array for these
CH.J.R.marray <- marray(CH.J.R2)
# The last column ought to show the number of juveniles not recaptured again and should all be zeros, since all of them are released as adults
CH.J.R.marray[,dim(CH.J)[2]] <- 0
# Create the m-array for juveniles never recaptured and add it to the previous m-array
CH.J.N.marray <- marray(CH.J.N)
CH.J.marray <- CH.J.R.marray + CH.J.N.marray 

# Specify model in BUGS language
sink("cjs-mnl-age.jags")
cat("
model {
# Priors and constraints
for (t in 1:(n.occasions-1)){
   phi.juv[t] <- mean.phijuv
   phi.ad[t] <- mean.phiad
   p[t] <- mean.p
   }
mean.phijuv ~ dunif(0, 1)          # Prior for mean juv. survival
mean.phiad ~ dunif(0, 1)           # Prior for mean ad. survival
mean.p ~ dunif(0, 1)               # Prior for mean recapture
# Define the multinomial likelihood
for (t in 1:(n.occasions-1)){
   marr.j[t,1:n.occasions] ~ dmulti(pr.j[t,], rel.j[t])
   marr.a[t,1:n.occasions] ~ dmulti(pr.a[t,], rel.a[t])
   }
# Define the cell probabilities of the m-arrays
# Main diagonal
for (t in 1:(n.occasions-1)){
   q[t] <- 1-p[t]            # Probability of non-recapture
   pr.j[t,t] <- phi.juv[t]*p[t]
   pr.a[t,t] <- phi.ad[t]*p[t]
   # Above main diagonal
   for (j in (t+1):(n.occasions-1)){
      pr.j[t,j] <- phi.juv[t]*prod(phi.ad[(t+1):j])*prod(q[t:(j-1)])*p[j]
      pr.a[t,j] <- prod(phi.ad[t:j])*prod(q[t:(j-1)])*p[j]
      } #j
   # Below main diagonal
   for (j in 1:(t-1)){
      pr.j[t,j] <- 0
      pr.a[t,j] <- 0
      } #j
   } #t
# Last column: probability of non-recapture
for (t in 1:(n.occasions-1)){
   pr.j[t,n.occasions] <- 1-sum(pr.j[t,1:(n.occasions-1)])
   pr.a[t,n.occasions] <- 1-sum(pr.a[t,1:(n.occasions-1)])
   } #t
}
",fill = TRUE)
sink()

# Bundle data
jags.data <- list(marr.j = CH.J.marray, marr.a = CH.A.marray, n.occasions = dim(CH.J.marray)[2], rel.j = rowSums(CH.J.marray), rel.a = rowSums(CH.A.marray))

# Initial values
inits <- function(){list(mean.phijuv = runif(1, 0, 1), mean.phiad = runif(1, 0, 1), mean.p = runif(1, 0, 1))}  

# Parameters monitored
parameters <- c("mean.phijuv", "mean.phiad", "mean.p")

# MCMC settings
ni <- 3000
nt <- 3
nb <- 1000
nc <- 3

# Call JAGS from R (BRT <1 min)
cjs.2 <- jags(jags.data, inits, parameters, "cjs-mnl-age.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

par(mfrow = c(1, 3), las = 1)
hist(cjs.2$BUGSoutput$sims.list$mean.phijuv, nclass = 30, col = "gray", main = "", xlab = "Juvenile survival", ylab = "Frequency")
abline(v = phi.juv, col = "red", lwd = 2)
hist(cjs.2$BUGSoutput$sims.list$mean.phiad, nclass = 30, col = "gray", main = "", xlab = "Adult survival", ylab = "")
abline(v = phi.ad, col = "red", lwd = 2)
hist(cjs.2$BUGSoutput$sims.list$mean.p, nclass = 30, col = "gray", main = "", xlab = "Recapture", ylab = "")
abline(v = p[1], col = "red", lwd = 2)
 

# 7.11. Analysis of a real data set: survival of female Leisler?s bats
m.leisleri <- matrix(c(4,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3,
0,5,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3,
0,0,9,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3,
0,0,0,10,2,0,0,0,0,0,0,0,0,0,0,0,0,0,5,
0,0,0,0,10,2,1,0,0,0,0,0,0,0,0,0,0,0,6,
0,0,0,0,0,15,0,0,0,0,0,0,0,0,0,0,0,0,6,
0,0,0,0,0,0,11,2,0,1,0,0,0,0,0,0,0,0,19,
0,0,0,0,0,0,0,12,1,1,0,0,0,0,0,0,0,0,6,
0,0,0,0,0,0,0,0,13,2,0,0,0,0,0,0,0,0,4,
0,0,0,0,0,0,0,0,0,14,0,0,0,0,0,0,0,0,6,
0,0,0,0,0,0,0,0,0,0,13,1,0,0,0,1,0,0,8,
0,0,0,0,0,0,0,0,0,0,0,15,3,1,0,0,0,0,12,
0,0,0,0,0,0,0,0,0,0,0,0,12,4,0,1,0,0,7,
0,0,0,0,0,0,0,0,0,0,0,0,0,19,2,0,0,0,3,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,28,1,0,0,4,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,22,7,2,21,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,12,2,21,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,14,18), ncol = 19, nrow = 18, byrow = TRUE)

# Specify model in BUGS language
sink("cjs-mnl-ran.jags")
cat("
model {
# Priors and constraints
for (t in 1:(n.occasions-1)){
   logit(phi[t]) <- mu + epsilon[t]
   epsilon[t] ~ dnorm(0, tau)
   p[t] <- mean.p
   }
mean.phi ~ dunif(0, 1)             # Prior for mean survival
mu <- log(mean.phi / (1-mean.phi)) # Logit transformation
sigma ~ dunif(0, 5)                # Prior for standard deviation
tau <- pow(sigma, -2)
sigma2 <- pow(sigma, 2)
# Temporal variance on real scale
sigma2.real <- sigma2 * pow(mean.phi, 2) * pow((1-mean.phi), 2) 
mean.p ~ dunif(0, 1)	           # Prior for mean recapture
# Define the multinomial likelihood
for (t in 1:(n.occasions-1)){
   marr[t,1:n.occasions] ~ dmulti(pr[t,], rel[t])
   }
# Define the cell probabilities of the m-array:
# Main diagonal
for (t in 1:(n.occasions-1)){
   q[t] <- 1-p[t]
   pr[t,t] <- phi[t]*p[t]	
   # Above main diagonal
   for (j in (t+1):(n.occasions-1)){
      pr[t,j] <- prod(phi[t:j])*prod(q[t:(j-1)])*p[j]
      } #j	
   # Below main diagonal
   for (j in 1:(t-1)){
      pr[t,j]<-0
      } #j
   } #t
# Last column: probability of non-recapture
for (t in 1:(n.occasions-1)){
   pr[t,n.occasions] <- 1-sum(pr[t,1:(n.occasions-1)])
   } # t
# Assess model fit using Freeman-Tukey statistic
# Compute fit statistics for observed data
for (t in 1:(n.occasions-1)){
   for (j in 1:n.occasions){
      expmarr[t,j] <- rel[t]*pr[t,j]
      E.org[t,j] <- pow((pow(marr[t,j], 0.5)-pow(expmarr[t,j], 0.5)), 2)
      }
   }
# Generate replicate data and compute fit stats from them
for (t in 1:(n.occasions-1)){
   marr.new[t,1:n.occasions] ~ dmulti(pr[t,], rel[t])
   for (j in 1:n.occasions){
      E.new[t,j] <- pow((pow(marr.new[t,j], 0.5)-pow(expmarr[t,j], 0.5)), 2)
      }
   }
fit <- sum(E.org[,])
fit.new <- sum(E.new[,])
}
",fill = TRUE)
sink()

# Bundle data
jags.data <- list(marr = m.leisleri, n.occasions = dim(m.leisleri)[2], rel = rowSums(m.leisleri))

# Initial values
inits <- function(){list(mean.phi = runif(1, 0, 1), sigma = runif(1, 0, 5), mean.p = runif(1, 0, 1))}  

# Parameters monitored
parameters <- c("phi", "mean.p", "mean.phi", "sigma2", "sigma2.real", "fit", "fit.new")

# MCMC settings
ni <- 5000
nt <- 3
nb <- 1000
nc <- 3

# Call JAGS from R (BRT 3 min)
leis.result <- jags(jags.data, inits, parameters, "cjs-mnl-ran.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

# Summarize posteriors
print(leis.result, digits = 3)

# Produce figure of female survival probabilities
par(mfrow = c(1, 2), las = 1, mar=c(4, 4, 2, 2), mgp = c(3, 1, 0))
lower <- upper <- numeric()
T <- dim(m.leisleri)[2]-1
for (t in 1:T){
   lower[t] <- quantile(leis.result$BUGSoutput$sims.list$phi[,t], 0.025)
   upper[t] <- quantile(leis.result$BUGSoutput$sims.list$phi[,t], 0.975)
   }
plot(y = leis.result$BUGSoutput$mean$phi, x = (1:T)+0.5, type = "b", pch = 16, ylim = c(0.3, 1), ylab = "Annual survival probability", xlab = "", axes = F)
axis(1, at = seq(1,(T+1),2), labels = seq(1990,2008,2))
axis(1, at = 1:(T+1), labels = rep("", T+1), tcl = -0.25)
axis(2, las = 1)
mtext("Year", 1, line = 2.25)
text(18, 0.975, "(a)")
segments((1:T)+0.5, lower, (1:T)+0.5, upper)
segments(1, leis.result$BUGSoutput$mean$mean.phi, T+1, leis.result$BUGSoutput$mean$mean.phi, lty = 2, col = "red", lwd = 2)
segments(1, quantile(leis.result$BUGSoutput$sims.list$mean.phi,0.025), T+1, quantile(leis.result$BUGSoutput$sims.list$mean.phi, 0.025), lty = 2, col = "red") 
segments(1, quantile(leis.result$BUGSoutput$sims.list$mean.phi, 0.975), T+1, quantile(leis.result$BUGSoutput$sims.list$mean.phi, 0.975), lty = 2, col = "red") 
hist(leis.result$BUGSoutput$sims.list$sigma2.real, nclass = 45, col = "gray", main = "", las = 1, xlab = "")
mtext("Temporal variance of survival", 1, line = 2.25)
text(0.18, 900, "(b)")

# Evaluation of fit
plot(leis.result$BUGSoutput$sims.list$fit, leis.result$BUGSoutput$sims.list$fit.new, main = "", xlab = "Discrepancy actual data", ylab = "Discrepancy replicate data", las = 1, ylim = c(10, 35), xlim = c(10, 35), frame = FALSE)  
abline(0, 1, col = "black")




##################################################################
#
# 8. Estimation of survival probabilities using mark-recovery data
#
#################################################################

# 8.2. The mark-recovery model as a state-space model
# 8.2.1. Simulation of mark-recovery data
# Define parameter values
n.occasions <- 14                 # Number of release occasions
marked <- rep(50, n.occasions)    # Annual number of newly marked individuals
s <- rep(0.8, n.occasions)
r <- rep(0.2, n.occasions)

# Define matrices with survival and recovery probabilities
S <- matrix(s, ncol = n.occasions, nrow = sum(marked))
R <- matrix(r, ncol = n.occasions, nrow = sum(marked))

# Define function to simulate mark-recovery data
simul.mr <- function(S, R, marked){
   n.occasions <- dim(S)[2]
   MR <- matrix(NA, ncol = n.occasions+1, nrow = sum(marked))
   # Define a vector with the occasion of marking
   mark.occ <- rep(1:n.occasions, marked)
   # Fill the CH matrix
   for (i in 1:sum(marked)){
      MR[i, mark.occ[i]] <- 1    # Write an 1 at the release occasion
      for (t in mark.occ[i]:n.occasions){
         # Bernoulli trial: has individual survived occasion? 
         sur <- rbinom(1, 1, S[i,t])
         if (sur==1) next    # If still alive, move to next occasion 
         # Bernoulli trial: has dead individual been recovered? 
         rp <- rbinom(1, 1, R[i,t])
         if (rp==0){
            MR[i,t+1] <- 0
            break
            }
         if (rp==1){
            MR[i,t+1] <- 1
            break
            }
      } #t
   } #i
      # Replace the NA in the file by 0
      MR[which(is.na(MR))] <- 0
      return(MR)
   }

# Execute function
MR <- simul.mr(S, R, marked)


# 8.2.2. Analysis of a model with constant parameters
# Create vector with occasion of marking 
get.first <- function(x) min(which(x!=0))
f <- apply(MR, 1, get.first)

# Specify model in BUGS language
sink("mr.ss.jags")
cat("
model {

# Priors and constraints
for (i in 1:nind){
   for (t in f[i]:(n.occasions-1)){
      s[i,t] <- mean.s
      r[i,t] <- mean.r
      } #t
   } #i
mean.s ~ dunif(0, 1)          # Prior for mean survival
mean.r ~ dunif(0, 1)          # Prior for mean recapture

# Likelihood 
for (i in 1:nind){
   # Define latent state at first capture
   z[i,f[i]] <- 1
   for (t in (f[i]+1):n.occasions){
      # State process
      z[i,t] ~ dbern(mu1[i,t])
      mu1[i,t] <- s[i,t-1] * z[i,t-1]
      # Observation process
      y[i,t] ~ dbern(mu2[i,t])
      mu2[i,t] <- r[i,t-1] * (z[i,t-1] - z[i,t])
      } #t
   } #i
}
",fill = TRUE)
sink()

# Define function to create a matrix with information about known latent state z
known.state.mr <- function(mr){
   state <- matrix(NA, nrow = dim(mr)[1], ncol = dim(mr)[2])
   rec <- which(rowSums(mr)==2)
   for (i in 1:length(rec)){
      n1 <- min(which(mr[rec[i],]==1))
      n2 <- max(which(mr[rec[i],]==1))
      state[rec[i],n1:n2] <- 1
      state[rec[i],n1] <- NA
      state[rec[i],n2:dim(mr)[2]] <- 0
      }
   return(state)
   }

# Bundle data
jags.data <- list(y = MR, f = f, nind = dim(MR)[1], n.occasions = dim(MR)[2], z = known.state.mr(MR))

# Define function to create a matrix of initial values for latent state z
mr.init.z <- function(mr){
   ch <- matrix(NA, nrow = dim(mr)[1], ncol = dim(mr)[2])
   rec <- which(rowSums(mr)==1)
   for (i in 1:length(rec)){
      n1 <- which(mr[rec[i],]==1)
      ch[rec[i],n1:dim(mr)[2]] <- 0
      ch[rec[i],n1] <- NA
      }
   return(ch)
   }

# Initial values
inits <- function(){list(z = mr.init.z(MR), mean.s = runif(1, 0, 1), mean.r = runif(1, 0, 1))}  

# Parameters monitored
parameters <- c("mean.s", "mean.r")

# MCMC settings
ni <- 5000
nt <- 6
nb <- 2000
nc <- 3

# Call JAGS from R (BRT 4 min)
mr.ss <- jags(jags.data, inits, parameters, "mr.ss.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

print(mr.ss, digits = 3)


# 8.3. The mark-recovery model fitted with the multinomial likelihood 
# 8.3.1. Constant parameters
# Define function to create an m-array based for mark-recovery (MR) data
marray.dead <- function(MR){
   nind <- dim(MR)[1]
   n.occasions <- dim(MR)[2]
   m.array <- matrix(data = 0, ncol = n.occasions+1, nrow = n.occasions)
   # Create vector with occasion of marking 
   get.first <- function(x) min(which(x!=0))
   f <- apply(MR, 1, get.first)
   # Calculate the number of released individuals at each time period
   first <- as.numeric(table(f))
   for (t in 1:n.occasions){
      m.array[t,1] <- first[t]
      }
   # Fill m-array with recovered individuals
   rec.ind <- which(apply(MR, 1, sum)==2)
   rec <- numeric()
   for (i in 1:length(rec.ind)){
      d <- which(MR[rec.ind[i],(f[rec.ind[i]]+1):n.occasions]==1)
      rec[i] <- d + f[rec.ind[i]]
      m.array[f[rec.ind[i]],rec[i]] <- m.array[f[rec.ind[i]],rec[i]] + 1
      }
   # Calculate the number of individuals that are never recovered
   for (t in 1:n.occasions){
      m.array[t,n.occasions+1] <- m.array[t,1]-sum(m.array[t,2:n.occasions])
      }
   out <- m.array[1:(n.occasions-1),2:(n.occasions+1)]
   return(out)
   }

marr <- marray.dead(MR)

# Specify model in BUGS language
sink("mr-mnl.jags")
cat("
model {

# Priors and constraints
for (t in 1:n.occasions){
   s[t] <- mean.s
   r[t] <- mean.r
   }
mean.s ~ dunif(0, 1)              # Prior for mean survival
mean.r ~ dunif(0, 1)              # Prior for mean recovery

# Define the multinomial likelihood
for (t in 1:n.occasions){
   marr[t,1:(n.occasions+1)] ~ dmulti(pr[t,], rel[t])
   }
# Define the cell probabilities of the m-array
# Main diagonal
for (t in 1:n.occasions){
   pr[t,t] <- (1-s[t])*r[t]
   # Above main diagonal
   for (j in (t+1):n.occasions){
      pr[t,j] <- prod(s[t:(j-1)])*(1-s[j])*r[j]
      } #j
   # Below main diagonal
   for (j in 1:(t-1)){
      pr[t,j] <- 0
      } #j
   } #t
# Last column: probability of non-recovery
for (t in 1:n.occasions){
   pr[t,n.occasions+1] <- 1-sum(pr[t,1:n.occasions])
   } #t
}
",fill = TRUE)
sink()

# Bundle data
jags.data <- list(marr = marr, n.occasions = dim(marr)[2]-1, rel = rowSums(marr))

# Initial values
inits <- function(){list(mean.s = runif(1, 0, 1), mean.r = runif(1, 0, 1))}  

# Parameters monitored
parameters <- c("mean.s", "mean.r")

# MCMC settings
ni <- 5000
nt <- 6
nb <- 2000
nc <- 3

# Call JAGS from R (BRT <1 min)
mr <- jags(jags.data, inits, parameters, "mr-mnl.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

# Summarize posteriors
print(mr, digits = 3)

par(mfrow = c(1, 2), las = 1)
hist(mr$BUGSoutput$sims.list$mean.s, nclass = 25, col = "gray", main = "", ylab = "Frequency", xlab = "Survival probability")
abline(v = 0.8, col = "red", lwd = 2)
hist(mr$BUGSoutput$sims.list$mean.r, nclass = 25, col = "gray", main = "", ylab = "", xlab = "Recovery probability")
abline(v = 0.2, col = "red", lwd = 2)

 
# 8.3.2. Age-dependent parameters
n.occasions <- 15                   # Number of occasions
marked.j <- rep(200, n.occasions)   # Annual number of newly marked young
marked.a <- rep(20, n.occasions)    # Annual number of newly marked adults
sjuv <- 0.3                         # Juvenile survival probability
sad <- 0.8                          # Adult survival probability
rjuv <- 0.25                        # Juvenile recovery probability
rad <- 0.15                         # Adult recovery probability
sj <- c(sjuv, rep(sad, n.occasions-1))
rj <- c(rjuv, rep(rad, n.occasions-1))

# Define matrices with survival and recovery probabilities
SJ <- matrix(0, ncol = n.occasions, nrow = sum(marked.j))
for (i in 1:length(marked.j)){
   SJ[(sum(marked.j[1:i])-marked.j[i]+1):sum(marked.j[1:i]),i:n.occasions] <- matrix(rep(sj[1:(n.occasions-i+1)],marked.j[i]), ncol = n.occasions-i+1, byrow = TRUE)
   }
SA <- matrix(sad, ncol = n.occasions, nrow = sum(marked.a))
RJ <- matrix(0, ncol = n.occasions, nrow = sum(marked.j))
for (i in 1:length(marked.j)){
   RJ[(sum(marked.j[1:i])-marked.j[i]+1):sum(marked.j[1:i]),i:n.occasions] <- matrix(rep(rj[1:(n.occasions-i+1)],marked.j[i]), ncol = n.occasions-i+1, byrow = TRUE)
   }
RA <- matrix(rad, ncol = n.occasions, nrow = sum(marked.a))

# Execute simulation function
MRj <- simul.mr(SJ, RJ, marked.j)
MRa <- simul.mr(SA, RA, marked.a)

# Summarize data in m-arrays
marr.j <- marray.dead(MRj)
marr.a <- marray.dead(MRa)

# Specify model in BUGS language
sink("mr-mnl-age.jags")
cat("
model {

# Priors and constraints
for (t in 1:n.occasions){
   sj[t] <- mean.sj
   sa[t] <- mean.sa
   rj[t] <- mean.rj
   ra[t] <- mean.ra
   }
mean.sj ~ dunif(0, 1)              # Prior for mean juv. survival
mean.sa ~ dunif(0, 1)              # Prior for mean ad. survival
mean.rj ~ dunif(0, 1)              # Prior for mean juv. recovery
mean.ra ~ dunif(0, 1)              # Prior for mean ad. recovery

# Define the multinomial likelihoods
# Calculate the number of birds released each year
for (t in 1:n.occasions){
   marr.j[t,1:(n.occasions+1)] ~ dmulti(pr.j[t,], rel.j[t])
   marr.a[t,1:(n.occasions+1)] ~ dmulti(pr.a[t,], rel.a[t])
   }
# Define the cell probabilities of the juvenile m-array
# Main diagonal
for (t in 1:n.occasions){
   pr.j[t,t] <- (1-sj[t])*rj[t]
   # Further above main diagonal
   for (j in (t+2):n.occasions){
      pr.j[t,j] <- sj[t]*prod(sa[(t+1):(j-1)])*(1-sa[j])*ra[j]
      } #j
   # Below main diagonal
   for (j in 1:(t-1)){
      pr.j[t,j] <- 0
      } #j
   } #t
for (t in 1:(n.occasions-1)){
   # One above main diagonal
   pr.j[t,t+1] <- sj[t]*(1-sa[t+1])*ra[t+1] 
   } #t
# Last column: probability of non-recovery
for (t in 1:n.occasions){
   pr.j[t,n.occasions+1] <- 1-sum(pr.j[t,1:n.occasions])
   } #t
# Define the cell probabilities of the adult m-array
# Main diagonal
for (t in 1:n.occasions){
   pr.a[t,t] <- (1-sa[t])*ra[t]
   # Above main diagonal
   for (j in (t+1):n.occasions){
      pr.a[t,j] <- prod(sa[t:(j-1)])*(1-sa[j])*ra[j]
      } #j
   # Below main diagonal
   for (j in 1:(t-1)){
      pr.a[t,j] <- 0
      } #j
   } #t
# Last column: probability of non-recovery
for (t in 1:n.occasions){
   pr.a[t,n.occasions+1] <- 1-sum(pr.a[t,1:n.occasions])
   } #t
}
",fill = TRUE)
sink()

# Bundle data
jags.data <- list(marr.j = marr.j, marr.a = marr.a, n.occasions = dim(marr.j)[2]-1, rel.j = rowSums(marr.j), rel.a = rowSums(marr.a))

# Initial values
inits <- function(){list(mean.sj = runif(1, 0, 1), mean.sa = runif(1, 0, 1), mean.rj = runif(1, 0, 1), mean.ra = runif(1, 0, 1))}  

# Parameters monitored
parameters <- c("mean.sj", "mean.rj", "mean.sa", "mean.ra")

# MCMC settings
ni <- 5000
nt <- 6
nb <- 2000
nc <- 3

# Call JAGS from R (BRT <1 min)
mr.age <- jags(jags.data, inits, parameters, "mr-mnl-age.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

print(mr.age, digits = 3)


# 8.4. Real data example: age-dependent survival in Swiss red kites
marray.juv <- c(42, 18, 5, 7, 4, 3, 2, 1, 2, 2, 1, 0, 1, 3, 0, 0, 1, 1388)
marray.ad <- c(3, 1, 1, 3, 0, 2, 1, 0, 1, 1, 0, 1, 1, 0, 0, 0, 0, 137)

# Specify model in BUGS language
sink("mr-mnl-age3.jags")
cat("
model {

# Priors and constraints
#sjuv ~ dbeta(4.2, 2.8)   # Informative prior for juv. survival: Analysis A
sjuv ~ dunif(0, 1)      # Non-informative for juv. survival prior: Analysis B
ssub ~ dunif(0, 1)       # Prior for subad. survival
sad ~ dunif(0, 1)        # Prior for ad. survival
rjuv ~ dunif(0, 1)       # Prior for juv. recovery
rad ~ dunif(0, 1)        # Prior for ad. recovery

# Define the multinomial likelihoods
marr.j[1:(n.age+1)] ~ dmulti(pr.j[], rel.j)
marr.a[1:(n.age+1)] ~ dmulti(pr.a[], rel.a)
# Define the cell probabilities of the juvenile m-array
# First element
pr.j[1] <- (1-sjuv)*rjuv
# Second element
pr.j[2] <- sjuv*(1-ssub)*rad
# Third and further elements
for (t in 3:n.age){
   pr.j[t] <- sjuv*ssub*pow(sad,(t-3))*(1-sad)*rad
   }
# Probability of non-recovery
pr.j[n.age+1] <- 1-sum(pr.j[1:n.age])
# Define the cell probabilities of the adult m-array
# All elements
for (t in 1:n.age){
   pr.a[t] <- pow(sad,(t-1))*(1-sad)*rad
   }
# Probability of non-recovery
pr.a[n.age+1] <- 1-sum(pr.a[1:n.age])
}
",fill = TRUE)
sink()

# Bundle data
jags.data <- list(marr.j = marray.juv, marr.a = marray.ad, n.age = length(marray.juv)-1, rel.j = sum(marray.juv), rel.a = sum(marray.ad))

# Initial values
inits <- function(){list(sjuv = runif(1, 0, 1), ssub = runif(1, 0, 1), sad = runif(1, 0, 1), rjuv = runif(1, 0, 1), rad = runif(1, 0, 1))}  

# Parameters monitored
parameters <- c("sjuv", "ssub", "sad", "rjuv", "rad")

# MCMC settings
ni <- 30000
nt <- 10
nb <- 10000
nc <- 3

# Call JAGS from R (BRT <1 min)
rk.ageA <- jags(jags.data, inits, parameters, "mr-mnl-age3.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())
rk.ageB <- jags(jags.data, inits, parameters, "mr-mnl-age3.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

print(rk.ageA, digits = 3)
print(rk.ageB, digits = 3)

par(mfrow = c(2, 3), las = 1)
hist(rk.ageB$BUGSoutput$sims.list$sjuv, breaks = 20, col = "gray", main = "", xlab = "Juvenile survival")
hist(rk.ageB$BUGSoutput$sims.list$ssub, breaks = 20, col = "gray", main = "", xlab = "Subadult survival")
hist(rk.ageB$BUGSoutput$sims.list$sad, breaks = 20, col = "gray", main = "", xlab = "Adult survival")
hist(rk.ageB$BUGSoutput$sims.list$rjuv, breaks = 20, col = "gray", main = "", xlab = "Juvenile recovery", xlim = c(0, 0.2))
hist(rk.ageB$BUGSoutput$sims.list$rad, breaks = 20, col = "gray", main = "", xlab = "Adult recovery")

plot(density(rk.ageA$BUGSoutput$sims.list$sjuv), ylim = c(0, 5), lwd = 2, main = "", xlab = "Juvenile survival", las = 1)
points(density(rk.ageB$BUGSoutput$sims.list$sjuv), col = "red", type = "l", lwd = 2)
text(x = 0.5, y = 4.8, "Prior distributions", pos = 4, font = 3)
legend(x = 0.6, y = 4.7, legend = c("U(0,1)", "beta(4.2,2.8)"), lwd = c(2, 2), col = c("black", "red"), bty = "n")

quantile(rk.ageA$BUGSoutput$sims.list$ssub-rk.ageA$BUGSoutput$sims.list$sjuv, prob = c(0.025, 0.975))
quantile(rk.ageA$BUGSoutput$sims.list$sad-rk.ageA$BUGSoutput$sims.list$ssub, prob = c(0.025, 0.975))




#############################################################
#
# 9. Multistate capture-recapture models
#
##############################################################

# 9.2. Estimation of movement between two sites
# 9.2.1. Model description
# 9.2.2. Generation of simulated data
# Define mean survival, transitions, recapture, as well as number of occasions, states, observations and released individuals
phiA <- 0.8
phiB <- 0.7
psiAB <- 0.3
psiBA <- 0.5
pA <- 0.7
pB <- 0.4
n.occasions <- 6
n.states <- 3
n.obs <- 3
marked <- matrix(NA, ncol = n.states, nrow = n.occasions)
marked[,1] <- rep(100, n.occasions)  
marked[,2] <- rep(60, n.occasions)
marked[,3] <- rep(0, n.occasions)

# Define matrices with survival, transition and recapture probabilities
# These are 4-dimensional matrices, with 
   # Dimension 1: state of departure
   # Dimension 2: state of arrival
   # Dimension 3: individual
   # Dimension 4: time

# 1. State process matrix
totrel <- sum(marked)*(n.occasions-1)
PSI.STATE <- array(NA, dim=c(n.states, n.states, totrel, n.occasions-1))
for (i in 1:totrel){
   for (t in 1:(n.occasions-1)){
      PSI.STATE[,,i,t] <- matrix(c(
      phiA*(1-psiAB), phiA*psiAB,     1-phiA,
      phiB*psiBA,     phiB*(1-psiBA), 1-phiB,
      0,              0,              1       ), nrow = n.states, byrow = TRUE)
      } #t
   } #i

# 2.Observation process matrix
PSI.OBS <- array(NA, dim=c(n.states, n.obs, totrel, n.occasions-1))
for (i in 1:totrel){
   for (t in 1:(n.occasions-1)){
      PSI.OBS[,,i,t] <- matrix(c(
      pA, 0,  1-pA,
      0,  pB, 1-pB,
      0,  0,  1       ), nrow = n.states, byrow = TRUE)
      } #t
   } #i

# Define function to simulate multistate capture-recapture data
simul.ms <- function(PSI.STATE, PSI.OBS, marked, unobservable = NA){
   # Unobservable: number of state that is unobservable
   n.occasions <- dim(PSI.STATE)[4] + 1
   CH <- CH.TRUE <- matrix(NA, ncol = n.occasions, nrow = sum(marked))
   # Define a vector with the occasion of marking
   mark.occ <- matrix(0, ncol = dim(PSI.STATE)[1], nrow = sum(marked))
   g <- colSums(marked)
   for (s in 1:dim(PSI.STATE)[1]){
      if (g[s]==0) next  # To avoid error message if nothing to replace
      mark.occ[(cumsum(g[1:s])-g[s]+1)[s]:cumsum(g[1:s])[s],s] <-
      rep(1:n.occasions, marked[1:n.occasions,s])
      } #s
   for (i in 1:sum(marked)){
      for (s in 1:dim(PSI.STATE)[1]){
         if (mark.occ[i,s]==0) next
         first <- mark.occ[i,s]
         CH[i,first] <- s
         CH.TRUE[i,first] <- s
         } #s
      for (t in (first+1):n.occasions){
         # Multinomial trials for state transitions
         if (first==n.occasions) next
         state <- which(rmultinom(1, 1, PSI.STATE[CH.TRUE[i,t-1],,i,t-1])==1)
         CH.TRUE[i,t] <- state
         # Multinomial trials for observation process
         event <- which(rmultinom(1, 1, PSI.OBS[CH.TRUE[i,t],,i,t-1])==1)
         CH[i,t] <- event
         } #t
      } #i
   # Replace the NA and the highest state number (dead) in the file by 0
   CH[is.na(CH)] <- 0
   CH[CH==dim(PSI.STATE)[1]] <- 0
   CH[CH==unobservable] <- 0
   id <- numeric(0)
   for (i in 1:dim(CH)[1]){
      z <- min(which(CH[i,]!=0))
      ifelse(z==dim(CH)[2], id <- c(id,i), id <- c(id))
      }
   return(list(CH=CH[-id,], CH.TRUE=CH.TRUE[-id,]))
   # CH: capture histories to be used
   # CH.TRUE: capture histories with perfect observation
   }

# Execute function
sim <- simul.ms(PSI.STATE, PSI.OBS, marked)
CH <- sim$CH

# Compute vector with occasion of first capture
get.first <- function(x) min(which(x!=0))
f <- apply(CH, 1, get.first)

# Recode CH matrix: note, a 0 is not allowed in WinBUGS!
# 1 = seen alive in A, 2 = seen alive in B, 3 = not seen
rCH <- CH          # Recoded CH
rCH[rCH==0] <- 3


# 9.2.3. Analysis of the model
# Specify model in BUGS language
sink("ms.jags")
cat("
model {

# -------------------------------------------------
# Parameters:
# phiA: survival probability at site A
# phiB: survival probability at site B
# psiAB: movement probability from site A to site B
# psiBA: movement probability from site B to site A
# pA: recapture probability at site A
# pB: recapture probability at site B
# -------------------------------------------------
# States (S):
# 1 alive at A
# 2 alive at B
# 3 dead
# Observations (O):  
# 1 seen at A 
# 2 seen at B
# 3 not seen
# -------------------------------------------------

# Priors and constraints
for (t in 1:(n.occasions-1)){
   phiA[t] <- mean.phi[1]
   phiB[t] <- mean.phi[2]
   psiAB[t] <- mean.psi[1]
   psiBA[t] <- mean.psi[2]
   pA[t] <- mean.p[1]
   pB[t] <- mean.p[2]
   }
for (u in 1:2){
   mean.phi[u] ~ dunif(0, 1)    # Priors for mean state-spec. survival
   mean.psi[u] ~ dunif(0, 1)    # Priors for mean transitions
   mean.p[u] ~ dunif(0, 1)      # Priors for mean state-spec. recapture
   }

# Define state-transition and observation matrices
for (i in 1:nind){  
   # Define probabilities of state S(t+1) given S(t)
   for (t in f[i]:(n.occasions-1)){
      ps[1,i,t,1] <- phiA[t] * (1-psiAB[t])
      ps[1,i,t,2] <- phiA[t] * psiAB[t]
      ps[1,i,t,3] <- 1-phiA[t]
      ps[2,i,t,1] <- phiB[t] * psiBA[t]
      ps[2,i,t,2] <- phiB[t] * (1-psiBA[t])
      ps[2,i,t,3] <- 1-phiB[t]
      ps[3,i,t,1] <- 0
      ps[3,i,t,2] <- 0
      ps[3,i,t,3] <- 1
      
      # Define probabilities of O(t) given S(t)
      po[1,i,t,1] <- pA[t]
      po[1,i,t,2] <- 0
      po[1,i,t,3] <- 1-pA[t]
      po[2,i,t,1] <- 0
      po[2,i,t,2] <- pB[t]
      po[2,i,t,3] <- 1-pB[t]
      po[3,i,t,1] <- 0
      po[3,i,t,2] <- 0
      po[3,i,t,3] <- 1
      } #t
   } #i

# Likelihood 
for (i in 1:nind){
   # Define latent state at first capture
   z[i,f[i]] <- y[i,f[i]]
   for (t in (f[i]+1):n.occasions){
      # State process: draw S(t) given S(t-1)
      z[i,t] ~ dcat(ps[z[i,t-1], i, t-1,])
      # Observation process: draw O(t) given S(t)
      y[i,t] ~ dcat(po[z[i,t], i, t-1,])
      } #t
   } #i
}
",fill = TRUE)
sink()

# Function to create known latent states z
known.state.ms <- function(ms, notseen){
   # notseen: label for ?not seen?
   state <- ms
   state[state==notseen] <- NA
   for (i in 1:dim(ms)[1]){
      m <- min(which(!is.na(state[i,])))
      state[i,m] <- NA
      }
   return(state)
   }

# Function to create initial values for unknown z
ms.init.z <- function(ch, f){
   for (i in 1:dim(ch)[1]){ch[i,1:f[i]] <- NA}
   states <- max(ch, na.rm = TRUE)
   known.states <- 1:(states-1)
   v <- which(ch==states)
   ch[-v] <- NA
   ch[v] <- sample(known.states, length(v), replace = TRUE)
   return(ch)
   }

# Bundle data
jags.data <- list(y = rCH, f = f, n.occasions = dim(rCH)[2], nind = dim(rCH)[1], z = known.state.ms(rCH, 3))

# Initial values
inits <- function(){list(mean.phi = runif(2, 0, 1), mean.psi = runif(2, 0, 1), mean.p = runif(2, 0, 1), z = ms.init.z(rCH, f))}  

# Parameters monitored
parameters <- c("mean.phi", "mean.psi", "mean.p")

# MCMC settings
ni <- 10000
nt <- 6
nb <- 2000
nc <- 3

# Call JAGS from R (BRT 8 min)
ms <- jags(jags.data, inits, parameters, "ms.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())
print(ms, digits = 3)

par(mfrow = c(3, 2), las = 1)
hist(ms$BUGSoutput$sims.list$mean.phi[,1], col = "gray", main = "", xlab = expression(phi[A]), ylim=c(0,1300))
abline(v = phiA, col = "red")	
hist(ms$BUGSoutput$sims.list$mean.phi[,2], col = "gray", main = "", xlab = expression(phi[B]), ylim=c(0,1300), ylab="")
abline(v = phiB, , col="red")
hist(ms$BUGSoutput$sims.list$mean.psi[,1], col = "gray", main = "", xlab = expression(psi[AB]), ylim=c(0,1300))
abline(v = psiAB, col="red")
hist(ms$BUGSoutput$sims.list$mean.psi[,2], col = "gray", main = "", xlab = expression(psi[BA]), ylab="", ylim=c(0,1300))
abline(v = psiBA, col="red")
hist(ms$BUGSoutput$sims.list$mean.p[,1], col = "gray", main = "", xlab = expression(p[A]), ylim=c(0,1300))
abline(v = pA, col = "red")
hist(ms$BUGSoutput$sims.list$mean.p[,2], col = "gray", main = "", xlab = expression(p[B]), ylab="", ylim=c(0,1300))
abline(v = pB, col = "red")

# Specify model in BUGS language
sink("ms.alternative1.jags")
cat("
model {

# Priors and constraints
for (t in 1:(n.occasions-1)){
   phiA[t] <- mean.phi[1]
   phiB[t] <- mean.phi[2]
   psiAB[t] <- mean.psi[1]
   psiBA[t] <- mean.psi[2]
   pA[t] <- mean.p[1]
   pB[t] <- mean.p[2]
   }
for (u in 1:2){
   mean.phi[u] ~ dunif(0, 1)      # Priors for mean state-spec. survival
   mean.psi[u] ~ dunif(0, 1)      # Priors for mean transitions
   mean.p[u] ~ dunif(0, 1)        # Priors for mean state-spec. recapture
   }

# Define state-transition and observation matrices
   # Define probabilities of state S(t+1) given S(t)
   for (t in 1:(n.occasions-1)){
      ps[1,t,1] <- phiA[t] * (1-psiAB[t])
      ps[1,t,2] <- phiA[t] * psiAB[t]
      ps[1,t,3] <- 1-phiA[t]
      ps[2,t,1] <- phiB[t] * psiBA[t]
      ps[2,t,2] <- phiB[t] * (1-psiBA[t])
      ps[2,t,3] <- 1-phiB[t]
      ps[3,t,1] <- 0
      ps[3,t,2] <- 0
      ps[3,t,3] <- 1
      
   # Define probabilities of O(t) given S(t)
      po[1,t,1] <- pA[t]
      po[1,t,2] <- 0
      po[1,t,3] <- 1-pA[t]
      po[2,t,1] <- 0
      po[2,t,2] <- pB[t]
      po[2,t,3] <- 1-pB[t]
      po[3,t,1] <- 0
      po[3,t,2] <- 0
      po[3,t,3] <- 1
      } #t

# Likelihood 
for (i in 1:nind){
   # Define latent state at first capture
   z[i,f[i]] <- y[i,f[i]]
   for (t in (f[i]+1):n.occasions){
      # State process: draw S(t) given S(t-1)
      z[i,t] ~ dcat(ps[z[i,t-1], t-1,])
      # Observation process: draw O(t) given S(t)
      y[i,t] ~ dcat(po[z[i,t], t-1,])
      } #t
   } #i
}
",fill = TRUE)
sink()

# Specify model in BUGS language
sink("ms.alternative2.jags")
cat("
model {

# Priors and constraints	
   phiA ~ dunif(0, 1)    # Prior for mean survival in A
   phiB ~ dunif(0, 1)    # Prior for mean survival in B
   psiAB ~ dunif(0, 1)   # Prior for mean movement from A to B
   psiBA ~ dunif(0, 1)   # Prior for mean movement from B to A
   pA ~ dunif(0, 1)      # Prior for mean recapture in A
   pB ~ dunif(0, 1)      # Prior for mean recapture in B

# Define state-transition and observation matrices
   # Define probabilities of state S(t+1) given S(t)
      ps[1,1] <- phiA * (1-psiAB)
      ps[1,2] <- phiA * psiAB
      ps[1,3] <- 1-phiA
      ps[2,1] <- phiB * psiBA
      ps[2,2] <- phiB * (1-psiBA)
      ps[2,3] <- 1-phiB
      ps[3,1] <- 0
      ps[3,2] <- 0
      ps[3,3] <- 1
      
   # Define probabilities of O(t) given S(t)
      po[1,1] <- pA
      po[1,2] <- 0
      po[1,3] <- 1-pA
      po[2,1] <- 0
      po[2,2] <- pB
      po[2,3] <- 1-pB
      po[3,1] <- 0
      po[3,2] <- 0
      po[3,3] <- 1
      
# Likelihood 
for (i in 1:nind){
   # Define latent state at first capture
   z[i,f[i]] <- y[i,f[i]]
   for (t in (f[i]+1):n.occasions){
      # State process: draw S(t) given S(t-1)
      z[i,t] ~ dcat(ps[z[i,t-1],])
      # Observation process: draw O(t) given S(t)
      y[i,t] ~ dcat(po[z[i,t],])
      } #t
   } #i
}
",fill = TRUE)
sink()


# 9.3. Accounting for temporary emigration
# 9.3.1. Model description
# 9.3.2. Generation of simulated data
# Define mean survival, transitions, recapture, as well as number of occasions, states, observations and released individuals
phi <- 0.85
psiIO <- 0.2
psiOI <- 0.3
p <- 0.7
n.occasions <- 8  
n.states <- 3
n.obs <- 2
marked <- matrix(NA, ncol = n.states, nrow = n.occasions)
marked[,1] <- rep(70, n.occasions)	# Present
marked[,2] <- rep(0, n.occasions)	# Absent
marked[,3] <- rep(0, n.occasions)	# Dead 

# Define matrices with survival, transition and recapture probabilities
# These are 4-dimensional matrices, with 
   # Dimension 1: state of departure
   # Dimension 2: state of arrival
   # Dimension 3: individual
   # Dimension 4: time
# 1. State process matrix
totrel <- sum(marked)*(n.occasions-1)
PSI.STATE <- array(NA, dim=c(n.states, n.states, totrel, n.occasions-1))
for (i in 1:totrel){
   for (t in 1:(n.occasions-1)){
      PSI.STATE[,,i,t] <- matrix(c(
      phi*(1-psiIO), phi*psiIO,     1-phi,
      phi*psiOI,     phi*(1-psiOI), 1-phi,
      0,             0,             1       ), nrow = n.states, byrow = TRUE)
      } #t
   } #i

# 2.Observation process matrix
PSI.OBS <- array(NA, dim=c(n.states, n.obs, totrel, n.occasions-1))
for (i in 1:totrel){
   for (t in 1:(n.occasions-1)){
      PSI.OBS[,,i,t] <- matrix(c(
      p, 1-p,
      0, 1,
      0, 1    ), nrow = n.states, byrow = TRUE)
      } #t
   } #i

# Execute simulation function
sim <- simul.ms(PSI.STATE, PSI.OBS, marked)
CH <- sim$CH

# Compute vector with occasion of first capture
get.first <- function(x) min(which(x!=0))
f <- apply(CH, 1, get.first)

# Recode CH matrix: note, a 0 is not allowed!
# 1 = seen alive, 2 = not seen
rCH <- CH  # Recoded CH
rCH[rCH==0] <- 2

 
# 9.3.3. Analysis of the model
# Specify model in BUGS language
sink("tempemi.jags")
cat("
model {

# ---------------------------------
# Parameters:
# phi: survival probability
# psiIO: probability to emigrate
# psiOI: probability to immigrate
# p: recapture probability
# ---------------------------------
# States (S):
# 1 alive and present
# 2 alive and absent
# 3 dead
# Observations (O):
# 1 seen 
# 2 not seen
# ---------------------------------

# Priors and constraints
for (t in 1:(n.occasions-1)){
   phi[t] <- mean.phi
   psiIO[t] <- mean.psiIO
   psiOI[t] <- mean.psiOI
   p[t] <- mean.p
   }
mean.phi ~ dunif(0, 1)       # Prior for mean survival
mean.psiIO ~ dunif(0, 1)     # Prior for mean temp. emigration
mean.psiOI ~ dunif(0, 1)     # Prior for mean temp. immigration
mean.p ~ dunif(0, 1)         # Prior for mean recapture

# Define state-transition and observation matrices 	
for (i in 1:nind){
   # Define probabilities of state S(t+1) given S(t)
   for (t in f[i]:(n.occasions-1)){
      ps[1,i,t,1] <- phi[t] * (1-psiIO[t])
      ps[1,i,t,2] <- phi[t] * psiIO[t]
      ps[1,i,t,3] <- 1-phi[t]
      ps[2,i,t,1] <- phi[t] * psiOI[t]
      ps[2,i,t,2] <- phi[t] * (1-psiOI[t])
      ps[2,i,t,3] <- 1-phi[t]
      ps[3,i,t,1] <- 0
      ps[3,i,t,2] <- 0
      ps[3,i,t,3] <- 1

      # Define probabilities of O(t) given S(t)
      po[1,i,t,1] <- p[t]
      po[1,i,t,2] <- 1-p[t]
      po[2,i,t,1] <- 0
      po[2,i,t,2] <- 1
      po[3,i,t,1] <- 0
      po[3,i,t,2] <- 1
      } #t
   } #i

# Likelihood
for (i in 1:nind){
   # Define latent state at first capture
   z[i,f[i]] <- y[i,f[i]]
   for (t in (f[i]+1):n.occasions){
      # State process: draw S(t) given S(t-1)
      z[i,t] ~ dcat(ps[z[i,t-1], i, t-1,])
      # Observation process: draw O(t) given S(t)
      y[i,t] ~ dcat(po[z[i,t], i, t-1,])
      } #t
   } #i
}
",fill = TRUE)
sink()

# Bundle data
jags.data <- list(y = rCH, f = f, n.occasions = dim(rCH)[2], nind = dim(rCH)[1], z = known.state.ms(rCH, 2))

# Initial values
inits <- function(){list(mean.phi = runif(1, 0, 1), mean.psiIO = runif(1, 0, 1), mean.psiOI = runif(1, 0, 1), mean.p = runif(1, 0, 1), z = ms.init.z(rCH, f))}  

# Parameters monitored
parameters <- c("mean.phi", "mean.psiIO", "mean.psiOI", "mean.p")

# MCMC settings
ni <- 50000
nt <- 10
nb <- 10000
nc <- 3

# Call JAGS from R (BRT 35 min)
tempemi <- jags(jags.data, inits, parameters, "tempemi.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

print(tempemi, digits = 3)

par(mfrow = c(2, 2), las = 1)
hist(tempemi$BUGSoutput$sims.list$mean.phi, col = "gray", main = "",  xlab = expression(phi))
abline(v = phi, col = "red", lwd = 2)
hist(tempemi$BUGSoutput$sims.list$mean.psiIO, col = "gray", main = "", xlab =  expression(psi[IO]), ylab = "")
abline(v = psiIO, col = "red", lwd = 2)
hist(tempemi$BUGSoutput$sims.list$mean.psiOI, col = "gray", main = "", xlab =  expression(psi[OI]))
abline(v = psiOI, col = "red", lwd = 2)
hist(tempemi$BUGSoutput$sims.list$mean.p, col = "gray", main = "", xlab = expression(p), ylab = "")
abline(v = p, col = "red", lwd = 2)

 
# 9.4. Estimation of age-specific probability of first breeding
# 9.4.1. Model description
# 9.4.2. Generation of simulated data
# Define mean survival, transitions, recapture, as well as number of occasions, states, observations and released individuals
phi.1 <- 0.4
phi.2 <- 0.7
phi.ad <- 0.8
alpha.1 <- 0.2
alpha.2 <- 0.6
p.NB <- 0.5
p.B <- 0.7
n.occasions <- 7  
n.states <- 5
n.obs <- 4
marked <- matrix(0, ncol = n.states, nrow = n.occasions)
marked[,1] <- rep(100, n.occasions)	# Releases only as juveniles

# Define matrices with survival, transition and recapture probabilities
# These are 4-dimensional matrices, with 
   # Dimension 1: state of departure
   # Dimension 2: state of arrival
   # Dimension 3: individual
   # Dimension 4: time
# 1. State process matrix
totrel <- sum(marked)*(n.occasions-1)
PSI.STATE <- array(NA, dim=c(n.states, n.states, totrel, n.occasions-1))
for (i in 1:totrel){
   for (t in 1:(n.occasions-1)){
      PSI.STATE[,,i,t] <- matrix(c(
      0, phi.1*(1-alpha.1), 0,                 phi.1*alpha.1, 1-phi.1,
      0, 0,                 phi.2*(1-alpha.2), phi.2*alpha.2, 1-phi.2,
      0, 0,                 0,                 phi.ad,        1-phi.ad,
      0, 0,                 0,                 phi.ad,        1-phi.ad,
      0, 0,                 0,                 0,             1), nrow = n.states, byrow = TRUE)
      } #t
   } #i

# 2.Observation process matrix
PSI.OBS <- array(NA, dim=c(n.states, n.obs, totrel, n.occasions-1))
for (i in 1:totrel){
   for (t in 1:(n.occasions-1)){
      PSI.OBS[,,i,t] <- matrix(c(
      0, 0,    0,   1,
      0, p.NB, 0,   1-p.NB,
      0, p.NB, 0,   1-p.NB,
      0, 0,    p.B, 1-p.B,
      0, 0,    0,   1), nrow = n.states, byrow = TRUE)
      } #t
   } #i

# Execute simulation function
sim <- simul.ms(PSI.STATE, PSI.OBS, marked)
CH <- sim$CH

# Compute vector with occasion of first capture
get.first <- function(x) min(which(x!=0))
f <- apply(CH, 1, get.first)

# Recode CH matrix: note, a 0 is not allowed!
# 1 = seen as juv, 2 = seen no rep, 3 = seen rep, 4 = not seen
rCH <- CH
rCH[rCH==0] <- 4

# 9.4.3. Analysis of the model
# Specify model in BUGS language
sink("agerecruitment.jags")
cat("
model {

# -------------------------------------------------
# Parameters:
# phi.1: first year survival probability
# phi.2: second year survival probability
# phi.ad: adult survival probability
# alpha.1: probability to start breeding when 1 year old
# alpha.2: probability to start breeding when 2 years old
# p.NB: recapture probability of non-breeders
# p.B: recapture probability of breeders
# -------------------------------------------------
# States (S):
# 1 juvenile
# 2 not yet breeding at age 1 year
# 3 not yet breeding at age 2 years
# 4 breeder
# 5 dead
# Observations (O):
# 1 seen as juvenile
# 2 seen as not yet breeding
# 3 seen breeding
# 4 not seen
# -------------------------------------------------

# Priors and constraints	
for (t in 1:(n.occasions-1)){
   phi.1[t] <- mean.phi1
   phi.2[t] <- mean.phi2
   phi.ad[t] <- mean.phiad
   alpha.1[t] <- mean.alpha1
   alpha.2[t] <- mean.alpha2
   p.NB[t] <- mean.pNB
   p.B[t] <- mean.pB
   }
mean.phi1 ~ dunif(0, 1)     # Prior for mean 1y survival
mean.phi2 ~ dunif(0, 1)     # Prior for mean 2y survival
mean.phiad ~ dunif(0, 1)    # Prior for mean ad survival
mean.alpha1 ~ dunif(0, 1)   # Prior for mean 1y breeding prob.
mean.alpha2 ~ dunif(0, 1)   # Prior for mean 2y breeding prob.
mean.pNB ~ dunif(0, 1)      # Prior for mean recapture non-breeders
mean.pB ~ dunif(0, 1)       # Prior for mean recapture breeders

# Define state-transition and observation matrices 	
for (i in 1:nind){
   # Define probabilities of state S(t+1) given S(t)
   for (t in f[i]:(n.occasions-1)){
      ps[1,i,t,1] <- 0
      ps[1,i,t,2] <- phi.1[t] * (1-alpha.1[t])
      ps[1,i,t,3] <- 0
      ps[1,i,t,4] <- phi.1[t] * alpha.1[t]
      ps[1,i,t,5] <- 1-phi.1[t]
      ps[2,i,t,1] <- 0
      ps[2,i,t,2] <- 0
      ps[2,i,t,3] <- phi.2[t] * (1-alpha.2[t])
      ps[2,i,t,4] <- phi.2[t] * alpha.2[t]
      ps[2,i,t,5] <- 1-phi.2[t]
      ps[3,i,t,1] <- 0
      ps[3,i,t,2] <- 0
      ps[3,i,t,3] <- 0
      ps[3,i,t,4] <- phi.ad[t]
      ps[3,i,t,5] <- 1-phi.ad[t]
      ps[4,i,t,1] <- 0
      ps[4,i,t,2] <- 0
      ps[4,i,t,3] <- 0
      ps[4,i,t,4] <- phi.ad[t]
      ps[4,i,t,5] <- 1-phi.ad[t]
      ps[5,i,t,1] <- 0
      ps[5,i,t,2] <- 0
      ps[5,i,t,3] <- 0
      ps[5,i,t,4] <- 0
      ps[5,i,t,5] <- 1

      # Define probabilities of O(t) given S(t)
      po[1,i,t,1] <- 0
      po[1,i,t,2] <- 0
      po[1,i,t,3] <- 0
      po[1,i,t,4] <- 1
      po[2,i,t,1] <- 0
      po[2,i,t,2] <- p.NB[t]
      po[2,i,t,3] <- 0
      po[2,i,t,4] <- 1-p.NB[t]
      po[3,i,t,1] <- 0
      po[3,i,t,2] <- p.NB[t]
      po[3,i,t,3] <- 0
      po[3,i,t,4] <- 1-p.NB[t]
      po[4,i,t,1] <- 0
      po[4,i,t,2] <- 0
      po[4,i,t,3] <- p.B[t]
      po[4,i,t,4] <- 1-p.B[t]
      po[5,i,t,1] <- 0
      po[5,i,t,2] <- 0
      po[5,i,t,3] <- 0
      po[5,i,t,4] <- 1
      } #t
   } #i

# Likelihood 
for (i in 1:nind){
   # Define latent state at first capture
   z[i,f[i]] <- y[i,f[i]]
   for (t in (f[i]+1):n.occasions){
      # State process: draw S(t) given S(t-1)
      z[i,t] ~ dcat(ps[z[i,t-1], i, t-1,])
      # Observation process: draw O(t) given S(t)
      y[i,t] ~ dcat(po[z[i,t], i, t-1,])
      } #t
   } #i
}
",fill = TRUE)
sink()

# Bundle data
jags.data <- list(y = rCH, f = f, n.occasions = dim(rCH)[2], nind = dim(rCH)[1])


# Initial values
# We have to provide initial values for the true latent state z. Since the observed state numbers do not always match the true state numbers, we cannot use simple the observed multistate capture history as initial value for z. The following function creates initial values for z from our observation under the model we want to analyse.

agefirst.init <- function(ch, f){
   age <- array(NA, dim=dim(ch))
   for (i in 1:nrow(ch)){
      for (t in f[i]:ncol(ch)){
         age[i,t] <- min(c(t-f[i]+1, 4))
         }
      }
   ini <- array(NA, dim=dim(ch))
   for (i in 1:nrow(ch)){
      for (t in f[i]:ncol(ch)){
         if(ch[i,t]==1) ini[i,t] <- 1
         if(ch[i,t]==2&age[i,t]==2) ini[i,t] <- 2
         if(ch[i,t]==3&age[i,t]==2) ini[i,t] <- 4
         if(ch[i,t]==3&age[i,t]==3) ini[i,t] <- 4
         if(ch[i,t]==4&age[i,t]==4) ini[i,t] <- 4
         if(ch[i,t]==2&age[i,t]==3) ini[i,t] <- 3
         }
      }
   ini[which(is.na(ini))] <- age[which(is.na(ini))]
   for (i in 1:nrow(ch)){
      ini[i,f[i]] <- NA
      for (t in f[i]:(ncol(ch)-1)){
         if(ini[i,t]==4 & ini[i,t+1]==3) ini[i,t+1] <- 4
         }
      }
   return(ini)
   }

inits <- function(){list(mean.phi1 = runif(1, 0, 1), mean.phi2 = runif(1, 0, 1), mean.phiad = runif(1, 0, 1), mean.alpha1 = runif(1, 0, 1), mean.alpha2 = runif(1, 0, 1), mean.pNB = runif(1, 0, 1), mean.pB = runif(1, 0, 1), z = agefirst.init(rCH, f))}  

# Parameters monitored
parameters <- c("mean.phi1", "mean.phi2", "mean.phiad", "mean.alpha1", "mean.alpha2", "mean.pNB", "mean.pB")

# MCMC settings
ni <- 2000
nt <- 3
nb <- 1000
nc <- 3

# Call JAGS from R (BRT 2 min)
agefirst <- jags(jags.data, inits, parameters, "agerecruitment.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())


par(mfrow = c(3, 3), las = 1)
hist(agefirst$BUGSoutput$sims.list$mean.phi1, col = "gray", main = "",  xlab = expression(phi[1]))
abline(v = phi.1, col = "red", lwd = 2)
hist(agefirst$BUGSoutput$sims.list$mean.phi2, col = "gray", main = "",  xlab = expression(phi[2]), ylab = "")
abline(v = phi.2, col = "red", lwd = 2)
hist(agefirst$BUGSoutput$sims.list$mean.phiad, col = "gray", main = "",  xlab = expression(phi[ad]) , ylab = "")
abline(v = phi.ad, col = "red", lwd = 2)
hist(agefirst$BUGSoutput$sims.list$mean.alpha1, col = "gray", main = "",  xlab = expression(alpha[1]))
abline(v = alpha.1, col = "red", lwd = 2)
hist(agefirst$BUGSoutput$sims.list$mean.alpha2, col = "gray", main = "",  xlab = expression(alpha[2]) , ylab = "")
abline(v = alpha.2, col = "red", lwd = 2)
plot(0, type = "n", axes = F, ylab = "", xlab = "")
hist(agefirst$BUGSoutput$sims.list$mean.pNB, col = "gray", main = "",  xlab = expression(p[NB]))
abline(v = p.NB, col = "red", lwd = 2)
hist(agefirst$BUGSoutput$sims.list$mean.pB, col = "gray", main = "",  xlab = expression(p[B]) , ylab = "")
abline(v = p.B, col = "red", lwd = 2)


# Add-in for JAGS
# Since we have created a matrix with initial values for the true state z, we can use part of this information as data (see chapter 7.3.1) which can help with convergence and computing time). Here we give those initial values that are based on an actual observation. Since the first observation is deterministic, it must be excluded. The following code constructs the data matrix:

ch <- rCH
ch[ch==4] <- NA
z.known <- agefirst.init(rCH, f)
z.known[is.na(ch)] <- NA
for (i in 1:nrow(ch)){
   z.known[i,f[i]] <- NA
   }
z <- agefirst.init(rCH, f)
z[!is.na(z.known)] <- NA

# Bundle data
jags.data <- list(y = rCH, f = f, n.occasions = dim(rCH)[2], nind = dim(rCH)[1], z = z.known)

# Initial values
inits <- function(){list(mean.phi1 = runif(1, 0, 1), mean.phi2 = runif(1, 0, 1), mean.phiad = runif(1, 0, 1),
 mean.alpha1 = runif(1, 0, 1), mean.alpha2 = runif(1, 0, 1), mean.pNB = runif(1, 0, 1), 
mean.pB = runif(1, 0, 1), z = z)}  

# Parameters monitored
parameters <- c("mean.phi1", "mean.phi2", "mean.phiad", "mean.alpha1", "mean.alpha2", "mean.pNB", "mean.pB")

# MCMC settings
ni <- 2000
nt <- 3
nb <- 1000
nc <- 3

# Call JAGS from R (BRT 2 min)
agefirst <- jags(jags.data, inits, parameters, "agerecruitment.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())




################################################



# 9.5. Joint analysis of capture-recapture and mark-recovery data
# 9.5.1. Model description
# 9.5.2. Generation of simulated data
# Define mean survival, transitions, recapture, as well as number of occasions, states, observations and released individuals 
s <- 0.8
F <- 0.6
r <- 0.1
p <- 0.5
n.occasions <- 10  
n.states <- 4
n.obs <- 3
marked <- matrix(0, ncol = n.states, nrow = n.occasions)
marked[,1] <- rep(100, n.occasions)	# Releases in study area

# Define matrices with survival, transition and recapture probabilities
# These are 4-dimensional matrices, with 
   # Dimension 1: state of departure
   # Dimension 2: state of arrival
   # Dimension 3: individual
   # Dimension 4: time
# 1. State process matrix
totrel <- sum(marked)*(n.occasions-1)
PSI.STATE <- array(NA, dim=c(n.states, n.states, totrel, n.occasions-1))
for (i in 1:totrel){
   for (t in 1:(n.occasions-1)){
      PSI.STATE[,,i,t] <- matrix(c(
      s*F, s*(1-F), 1-s, 0,
      0,   s,       1-s, 0,
      0,   0,       0,   1,
      0,   0,       0,   1), nrow = n.states, byrow = TRUE)
      } #t
   } #i

# 2.Observation process matrix
PSI.OBS <- array(NA, dim=c(n.states, n.obs, totrel, n.occasions-1))
for (i in 1:totrel){
   for (t in 1:(n.occasions-1)){
      PSI.OBS[,,i,t] <- matrix(c(
      p, 0, 1-p,
      0, 0, 1,
      0, r, 1-r,
      0, 0, 1), nrow = n.states, byrow = TRUE)
      } #t
   } #i

# Execute simulation function
sim <- simul.ms(PSI.STATE, PSI.OBS, marked)
CH <- sim$CH

# Compute date of first capture
get.first <- function(x) min(which(x!=0))
f <- apply(CH, 1, get.first)

# Recode CH matrix: note, a 0 is not allowed!
# 1 = alive and in study are, 2 = recovered dead, 3 = not seen or recovered
rCH <- CH  # Recoded CH
rCH[rCH==0] <- 3


# 9.5.3. Analysis of the model
# Specify model in BUGS language
sink("lifedead.jags")
cat("
model {

# -------------------------------------------------
# Parameters:
# s: true survival probability
# F: fidelity probability
# r: recovery probability
# p: recapture/resighting probability
# -------------------------------------------------
# States (S):
# 1 alive in study area
# 2 alive outside study area
# 3 recently dead and recovered
# 4 recently dead, but not recovered, or dead (absorbing)
# Observations (O):
# 1 seen alive
# 2 recovered dead
# 3 neither seen nor recovered
# -------------------------------------------------

# Priors and constraints
for (t in 1:(n.occasions-1)){
   s[t] <- mean.s
   F[t] <- mean.f
   r[t] <- mean.r
   p[t] <- mean.p
   }
mean.s ~ dunif(0, 1)     # Prior for mean survival
mean.f ~ dunif(0, 1)     # Prior for mean fidelity
mean.r ~ dunif(0, 1)     # Prior for mean recovery
mean.p ~ dunif(0, 1)     # Prior for mean recapture

# Define state-transition and observation matrices 	
for (i in 1:nind){
   # Define probabilities of state S(t+1) given S(t)
   for (t in f[i]:(n.occasions-1)){
      ps[1,i,t,1] <- s[t]*F[t]
      ps[1,i,t,2] <- s[t]*(1-F[t])
      ps[1,i,t,3] <- (1-s[t])*r[t]
      ps[1,i,t,4] <- (1-s[t])*(1-r[t])
      ps[2,i,t,1] <- 0
      ps[2,i,t,2] <- s[t]
      ps[2,i,t,3] <- (1-s[t])*r[t]
      ps[2,i,t,4] <- (1-s[t])*(1-r[t])
      ps[3,i,t,1] <- 0
      ps[3,i,t,2] <- 0
      ps[3,i,t,3] <- 0
      ps[3,i,t,4] <- 1
      ps[4,i,t,1] <- 0
      ps[4,i,t,2] <- 0
      ps[4,i,t,3] <- 0
      ps[4,i,t,4] <- 1

      # Define probabilities of O(t) given S(t)
      po[1,i,t,1] <- p[t]
      po[1,i,t,2] <- 0
      po[1,i,t,3] <- 1-p[t]
      po[2,i,t,1] <- 0
      po[2,i,t,2] <- 0
      po[2,i,t,3] <- 1
      po[3,i,t,1] <- 0
      po[3,i,t,2] <- 1
      po[3,i,t,3] <- 0
      po[4,i,t,1] <- 0
      po[4,i,t,2] <- 0
      po[4,i,t,3] <- 1
      } #t
   } #i

# Likelihood 
for (i in 1:nind){
   # Define latent state at first capture
   z[i,f[i]] <- y[i,f[i]]
   for (t in (f[i]+1):n.occasions){
      # State process: draw S(t) given S(t-1)
      z[i,t] ~ dcat(ps[z[i,t-1], i, t-1,])
      # Observation process: draw O(t) given S(t)
      y[i,t] ~ dcat(po[z[i,t], i, t-1,])
      } #t
   } #i
}
",fill = TRUE)
sink()


# Bundle data
jags.data <- list(y = rCH, f = f, n.occasions = dim(rCH)[2], nind = dim(rCH)[1])

# Initial values
# In order to run the model, we must provide sound initial values for the true latent state z. The difficulty is that observed states do not alsways correspond to the true latent state. For example, in our observation the state 2 refers to an individuals whose ring has been reported, while the true state 2 refers to an individuals that is avlive, but outside the study area. Consequently, we cannot use the observed states as the initial values for the true state in JAGS (in BUGS this works well). The function known.ld provides initial values for the state z for our model. The important things are i) that the observations correspond to the true state (i.e. all "2" are converted into "3"), ii) that states after a "3" are all "4", iii) that all non-observations between "1" become "1", and iv) that all remaining originally "3" after the first observation become "1".

ld.init <- function(ch, f){
   ch[ch==3] <- NA
   v2 <- which(ch==2, arr.ind = T)
   ch[v2] <- 3
   for (i in 1:nrow(v2)){
      ifelse(v2[i,2]!=ncol(ch), ch[v2[i,1], (v2[i,2]+1):ncol(ch)] <- 4, next)}
   for (i in 1:nrow(ch)){
      m <- max(which(ch[i,]==1))
      ch[i,f[i]:m] <- 1
      }
   for (i in 1:nrow(v2)){
      u1 <- min(which(ch[v2[i,1],]==1))
      ch[v2[i],u1:(v2[i,2]-1)] <- 1
      }
   for (i in 1:nrow(ch)){
      for (j in f[i]:ncol(ch)){
         if(is.na(ch[i,j])==1) ch[i,j] <- 1
         }
      ch[i,f[i]] <- NA
      }
   return(ch)
   }

inits <- function(){list(mean.s = runif(1, 0, 1), mean.f = runif(1, 0, 1), mean.p = runif(1, 0, 1), mean.r = runif(1, 0, 1), z = ld.init(rCH, f))}  

# Parameters monitored
parameters <- c("mean.s", "mean.f", "mean.r", "mean.p")

# MCMC settings
ni <- 4000
nt <- 3
nb <- 1000
nc <- 3

# Call JAGS from R (BRT 80 min)
lifedead <- jags(jags.data, inits, parameters, "lifedead.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

print(lifedead, digit = 3)

# Note that convergence is hard to get, much longer chains or more informative priors would be necessary to get convergence quicker.


# Add-in for JAGS
# Since we have created a matrix with initial values for the true state z, we can use part of this information as data (see chapter 7.3.1) which can help with convergence and computing time). Here we give those initial values that are based on an actual observation. Since the first observation is deterministic, it must be excluded. The following code constructs the data matrix:

ch <- rCH
ch[ch==3] <- NA
z.known <- ld.init(rCH, f)
z.known[is.na(ch)] <- NA
for (i in 1:nrow(ch)){
   z.known[i,f[i]] <- NA
   }
z <- ld.init(rCH, f)
z[!is.na(ch)] <- NA

# Bundle data
jags.data <- list(y = rCH, f = f, n.occasions = dim(rCH)[2], nind = dim(rCH)[1], z = z.known)

inits <- function(){list(mean.s = runif(1, 0, 1), mean.f = runif(1, 0, 1), mean.p = runif(1, 0, 1), mean.r = runif(1, 0, 1), z = z)}  

# Parameters monitored
parameters <- c("mean.s", "mean.f", "mean.r", "mean.p")

# MCMC settings
ni <- 4000
nt <- 3
nb <- 1000
nc <- 3

# Call JAGS from R (BRT 80 min)
lifedead <- jags(jags.data, inits, parameters, "lifedead.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

print(lifedead, digit = 3)



##############################################################



# 9.6. Estimation of movement among three sites
# 9.6.1. Model description
# 9.6.2. Generation of simulated data

# Define mean survival, transitions, recapture, as well as number of occasions, states, observations and released individuals
phiA <- 0.85
phiB <- 0.75
phiC <- 0.65
psiAB <- 0.3
psiAC <- 0.2
psiBA <- 0.5
psiBC <- 0.1
psiCA <- 0.6
psiCB <- 0.1
pA <- 0.7
pB <- 0.4
pC <- 0.5
n.occasions <- 6
n.states <- 4
n.obs <- 4
marked <- matrix(NA, ncol = n.states, nrow = n.occasions)
marked[,1] <- rep(50, n.occasions)  
marked[,2] <- rep(50, n.occasions)
marked[,3] <- rep(50, n.occasions)
marked[,4] <- rep(0, n.occasions)

# Define matrices with survival, transition and recapture probabilities
# These are 4-dimensional matrices, with 
   # Dimension 1: state of departure
   # Dimension 2: state of arrival
   # Dimension 3: individual
   # Dimension 4: time
# 1. State process matrix
totrel <- sum(marked)*(n.occasions-1)
PSI.STATE <- array(NA, dim=c(n.states, n.states, totrel, n.occasions-1))
for (i in 1:totrel){
   for (t in 1:(n.occasions-1)){
      PSI.STATE[,,i,t] <- matrix(c(
      phiA*(1-psiAB-psiAC), phiA*psiAB,           phiA*psiAC,           1-phiA,
      phiB*psiBA,           phiB*(1-psiBA-psiBC), phiB*psiBC,           1-phiB,
      phiC*psiCA,           phiC*psiCB,           phiC*(1-psiCA-psiCB), 1-phiC,
      0,                    0,                    0,                    1), nrow = n.states, byrow = TRUE)
      } #t
   } #i

# 2.Observation process matrix
PSI.OBS <- array(NA, dim=c(n.states, n.obs, totrel, n.occasions-1))
for (i in 1:totrel){
   for (t in 1:(n.occasions-1)){
      PSI.OBS[,,i,t] <- matrix(c(
      pA, 0,  0,  1-pA,
      0,  pB, 0,  1-pB,
      0,  0,  pC, 1-pC,
      0,  0,  0,  1), nrow = n.states, byrow = TRUE)
      } #t
   } #i

# Execute simulation function
sim <- simul.ms(PSI.STATE, PSI.OBS, marked)
CH <- sim$CH

# Compute vector with occasions of first capture
get.first <- function(x) min(which(x!=0))
f <- apply(CH, 1, get.first)

# Recode CH matrix: note, a 0 is not allowed in WinBUGS!
# 1 = seen alive in A, 2 = seen alive in B, 3, seen alive in C, 4 = not seen
rCH <- CH  # Recoded CH
rCH[rCH==0] <- 4


# 9.6.3. Analysis of the model
# Specify model in BUGS language
sink("ms3-multinomlogit.jags")
cat("
model {

# -------------------------------------------------
# Parameters:
# phiA: survival probability at site A
# phiB: survival probability at site B
# phiC: survival probability at site C
# psiAB: movement probability from site A to site B
# psiAC: movement probability from site A to site C
# psiBA: movement probability from site B to site A
# psiBC: movement probability from site B to site C
# psiCA: movement probability from site C to site A
# psiCB: movement probability from site C to site B
# pA: recapture probability at site A
# pB: recapture probability at site B
# pC: recapture probability at site C 
# -------------------------------------------------
# States (S):
# 1 alive at A
# 2 alive at B
# 3 alive at C
# 4 dead
# Observations (O):
# 1 seen at A 
# 2 seen at B
# 3 seen at C
# 4 not seen
# -------------------------------------------------

# Priors and constraints
   # Survival and recapture: uniform
   phiA ~ dunif(0, 1)
   phiB ~ dunif(0, 1)
   phiC ~ dunif(0, 1)
   pA ~ dunif(0, 1)
   pB ~ dunif(0, 1)
   pC ~ dunif(0, 1)

   # Transitions: multinomial logit
      # Normal priors on logit of all but one transition probas
      for (i in 1:2){
         lpsiA[i] ~ dnorm(0, 0.001)
         lpsiB[i] ~ dnorm(0, 0.001)
         lpsiC[i] ~ dnorm(0, 0.001)
         }
      # Constrain the transitions such that their sum is < 1
      for (i in 1:2){
         psiA[i] <- exp(lpsiA[i]) / (1 + exp(lpsiA[1]) + exp(lpsiA[2]))
         psiB[i] <- exp(lpsiB[i]) / (1 + exp(lpsiB[1]) + exp(lpsiB[2]))
         psiC[i] <- exp(lpsiC[i]) / (1 + exp(lpsiC[1]) + exp(lpsiC[2]))
         }
      # Calculate the last transition probability
      psiA[3] <- 1-psiA[1]-psiA[2]
      psiB[3] <- 1-psiB[1]-psiB[2]
      psiC[3] <- 1-psiC[1]-psiC[2]

# Define state-transition and observation matrices 	
for (i in 1:nind){
   # Define probabilities of state S(t+1) given S(t)
   for (t in f[i]:(n.occasions-1)){
      ps[1,i,t,1] <- phiA * psiA[1]
      ps[1,i,t,2] <- phiA * psiA[2]
      ps[1,i,t,3] <- phiA * psiA[3]
      ps[1,i,t,4] <- 1-phiA
      ps[2,i,t,1] <- phiB * psiB[1]
      ps[2,i,t,2] <- phiB * psiB[2]
      ps[2,i,t,3] <- phiB * psiB[3]
      ps[2,i,t,4] <- 1-phiB
      ps[3,i,t,1] <- phiC * psiC[1]
      ps[3,i,t,2] <- phiC * psiC[2]
      ps[3,i,t,3] <- phiC * psiC[3]
      ps[3,i,t,4] <- 1-phiC
      ps[4,i,t,1] <- 0
      ps[4,i,t,2] <- 0
      ps[4,i,t,3] <- 0
      ps[4,i,t,4] <- 1

      # Define probabilities of O(t) given S(t)
      po[1,i,t,1] <- pA
      po[1,i,t,2] <- 0
      po[1,i,t,3] <- 0
      po[1,i,t,4] <- 1-pA
      po[2,i,t,1] <- 0
      po[2,i,t,2] <- pB
      po[2,i,t,3] <- 0
      po[2,i,t,4] <- 1-pB
      po[3,i,t,1] <- 0
      po[3,i,t,2] <- 0
      po[3,i,t,3] <- pC
      po[3,i,t,4] <- 1-pC
      po[4,i,t,1] <- 0
      po[4,i,t,2] <- 0
      po[4,i,t,3] <- 0
      po[4,i,t,4] <- 1
      } #t
   } #i

# Likelihood 
for (i in 1:nind){
   # Define latent state at first capture   
   z[i,f[i]] <- y[i,f[i]]
   for (t in (f[i]+1):n.occasions){
      # State process: draw S(t) given S(t-1)
      z[i,t] ~ dcat(ps[z[i,t-1], i, t-1,])
      # Observation process: draw O(t) given S(t)
      y[i,t] ~ dcat(po[z[i,t], i, t-1,])
      } #t
   } #i
}
",fill = TRUE)
sink()

# Bundle data
jags.data <- list(y = rCH, f = f, n.occasions = dim(rCH)[2], nind = dim(rCH)[1], z = known.state.ms(rCH, 4))

# Initial values 
inits <- function(){list(phiA = runif(1, 0, 1), phiB = runif(1, 0, 1), phiC = runif(1, 0, 1), lpsiA = rnorm(2), lpsiB = rnorm(2), lpsiC = rnorm(2), pA = runif(1, 0, 1) , pB = runif(1, 0, 1) , pC = runif(1, 0, 1), z = ms.init.z(rCH, f))}  


# Parameters monitored
parameters <- c("phiA", "phiB", "phiC", "psiA", "psiB", "psiC", "pA", "pB", "pC")

# MCMC settings
ni <- 50000
nt <- 6
nb <- 20000
nc <- 3

# Call JAGS from R (BRT 56 min)
ms3 <- jags(jags.data, inits, parameters, "ms3-multinomlogit.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

print(ms3, digits = 3)




# 9.7. Real data example: the showy lady?s slipper
CH <- as.matrix(read.table("orchids.txt", sep=" ", header = F))
n.occasions <- dim(CH)[2]

# Compute vector with occasion of first capture
f <- numeric()
for (i in 1:dim(CH)[1]){f[i] <- min(which(CH[i,]!=0))}

# Recode CH matrix: note, a 0 is not allowed by WinBUGS!
# 1 = seen vegetative, 2 = seen flowering, 3 = not seen
rCH <- CH  # Recoded CH
rCH[rCH==0] <- 3

# Specify model in BUGS language
sink("ladyslipper.jags")
cat("
model {

# -------------------------------------------------
# Parameters:
# s: survival probability
# psiV: transitions from vegetative
# psiF: transitions from flowering
# psiD: transitions from dormant
# -------------------------------------------------
# States (S):
# 1 vegetative
# 2 flowering
# 3 dormant
# 4 dead
# Observations (O):  
# 1 seen vegetative 
# 2 seen flowering
# 3 not seen
# -------------------------------------------------

# Priors and constraints
   # Survival: uniform
   for (t in 1:(n.occasions-1)){  
      s[t] ~ dunif(0, 1)
      }
   # Transitions: gamma priors
   for (i in 1:3){
      a[i] ~ dgamma(1, 1)
      psiD[i] <- a[i]/sum(a[])
      b[i] ~ dgamma(1, 1)
      psiV[i] <- b[i]/sum(b[])
      c[i] ~ dgamma(1, 1)
      psiF[i] <- c[i]/sum(c[])
      }

# Define state-transition and observation matrices 	
for (i in 1:nind){
   # Define probabilities of state S(t+1) given S(t)
   for (t in 1:(n.occasions-1)){
      ps[1,i,t,1] <- s[t] * psiV[1]
      ps[1,i,t,2] <- s[t] * psiV[2]
      ps[1,i,t,3] <- s[t] * psiV[3]
      ps[1,i,t,4] <- 1-s[t]
      ps[2,i,t,1] <- s[t] * psiF[1]
      ps[2,i,t,2] <- s[t] * psiF[2]
      ps[2,i,t,3] <- s[t] * psiF[3]
      ps[2,i,t,4] <- 1-s[t]
      ps[3,i,t,1] <- s[t] * psiD[1]
      ps[3,i,t,2] <- s[t] * psiD[2]
      ps[3,i,t,3] <- s[t] * psiD[3]
      ps[3,i,t,4] <- 1-s[t]
      ps[4,i,t,1] <- 0
      ps[4,i,t,2] <- 0
      ps[4,i,t,3] <- 0
      ps[4,i,t,4] <- 1

      # Define probabilities of O(t) given S(t)
      po[1,i,t,1] <- 1
      po[1,i,t,2] <- 0
      po[1,i,t,3] <- 0
      po[2,i,t,1] <- 0
      po[2,i,t,2] <- 1
      po[2,i,t,3] <- 0
      po[3,i,t,1] <- 0
      po[3,i,t,2] <- 0
      po[3,i,t,3] <- 1
      po[4,i,t,1] <- 0
      po[4,i,t,2] <- 0
      po[4,i,t,3] <- 1
      } #t
   } #i

# Likelihood 
for (i in 1:nind){
   # Define latent state at first capture
   z[i,f[i]] <- y[i,f[i]]
   for (t in (f[i]+1):n.occasions){
      # State process: draw S(t) given S(t-1)
      z[i,t] ~ dcat(ps[z[i,t-1], i, t-1,])
      # Observation process: draw O(t) given S(t)
      y[i,t] ~ dcat(po[z[i,t], i, t-1,])
      } #t
   } #i
}
",fill=TRUE)
sink()

# Bundle data
jags.data <- list(y = rCH, f = f, n.occasions = dim(rCH)[2], nind = dim(rCH)[1], z = known.state.ms(rCH, 3))


# Initial values 
# Note that all initial values of the unknown z must all be 3 (not 1 or 2). Therefore the function ms.init.z need some slight changes.
ms.init.z <- function(ch, f){
   for (i in 1:dim(ch)[1]){ch[i,1:f[i]] <- NA}
   states <- max(ch, na.rm = TRUE)
   v <- which(ch==states)
   ch[-v] <- NA
   ch[v] <- states
   return(ch)
   }

inits <- function(){list(s = runif((dim(rCH)[2]-1), 0, 1), z = ms.init.z(rCH, f))}

# Parameters monitored
parameters <- c("s", "psiV", "psiF", "psiD")

# MCMC settings
ni <- 5000
nt <- 3
nb <- 2000
nc <- 3

# Call JAGS from R (BRT 3 min)
ls <- jags(jags.data, inits, parameters, "ladyslipper.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

print(ls, digits = 3)




#########################################################################################
#
# 10. Estimation of survival, recruitment and population size using the Jolly-Seber model
# 
##########################################################################################

# 10.3. Fitting the JS model with data augmentation
# 10.3.1. The JS model as a restricted dynamic occupancy model
# Specify model in BUGS language
sink("js-rest.occ.jags")
cat("
model {
# Priors and constraints
for (i in 1:M){
   for (t in 1:(n.occasions-1)){
      phi[i,t] <- mean.phi
      } #t
   for (t in 1:n.occasions){
      p[i,t] <- mean.p
      } #t
   } #i
mean.phi ~ dunif(0, 1)
mean.p ~ dunif(0, 1)

for (t in 1:n.occasions){
   gamma[t] ~ dunif(0, 1)
   } #t

# Likelihood
for (i in 1:M){
   # First occasion
   # State process
   z[i,1] ~ dbern(gamma[1])
   # Observation process
   mu1[i] <- z[i,1] * p[i,1]
   y[i,1] ~ dbern(mu1[i])
   # Subsequent occasions
   for (t in 2:n.occasions){
      # State process
      q[i,t-1] <- 1-z[i,t-1]		# Availability for recruitment
      mu2[i,t] <- phi[i,t-1] * z[i,t-1] + gamma[t] * prod(q[i,1:(t-1)]) 
      z[i,t] ~ dbern(mu2[i,t])
      # Observation process
      mu3[i,t] <- z[i,t] * p[i,t]
      y[i,t] ~ dbern(mu3[i,t])
      } #t
   } #i 

# Calculate derived population parameters
for (t in 1:n.occasions){
   qgamma[t] <- 1-gamma[t]
   }
cprob[1] <- gamma[1]
for (t in 2:n.occasions){
   cprob[t] <- gamma[t] * prod(qgamma[1:(t-1)])
   } #t
psi <- sum(cprob[])            # Inclusion probability
for (t in 1:n.occasions){
   b[t] <- cprob[t] / psi      # Entry probability
   } #t
for (i in 1:M){
   recruit[i,1] <- z[i,1]
   for (t in 2:n.occasions){
      recruit[i,t] <- (1-z[i,t-1]) * z[i,t]
      } #t
   } #i
for (t in 1:n.occasions){
   N[t] <- sum(z[1:M,t])        # Actual population size
   B[t] <- sum(recruit[1:M,t])  # Number of entries
   } #t
for (i in 1:M){
   Nind[i] <- sum(z[i,1:n.occasions])
   Nalive[i] <- 1-equals(Nind[i], 0)
   } #i
Nsuper <- sum(Nalive[])         # Superpopulation size
}
",fill=TRUE)
sink()


# 10.3.2. The JS model as a multistate model
# Specify model in BUGS language
sink("js-ms.jags")
cat("
model {

#--------------------------------------
# Parameters:
# phi: survival probability
# gamma: removal entry probability
# p: capture probability
#--------------------------------------
# States (S):
# 1 not yet entered
# 2 alive
# 3 dead
# Observations (O):
# 1 seen 
# 2 not seen
#--------------------------------------

# Priors and constraints
for (t in 1:(n.occasions-1)){
   phi[t] <- mean.phi
   gamma[t] ~ dunif(0, 1) # Prior for entry probabilities
   p[t] <- mean.p
   }

mean.phi ~ dunif(0, 1)    # Prior for mean survival
mean.p ~ dunif(0, 1)      # Prior for mean capture

# Define state-transition and observation matrices 	
for (i in 1:M){  
   # Define probabilities of state S(t+1) given S(t)
   for (t in 1:(n.occasions-1)){
      ps[1,i,t,1] <- 1-gamma[t]
      ps[1,i,t,2] <- gamma[t]
      ps[1,i,t,3] <- 0
      ps[2,i,t,1] <- 0
      ps[2,i,t,2] <- phi[t]
      ps[2,i,t,3] <- 1-phi[t]
      ps[3,i,t,1] <- 0
      ps[3,i,t,2] <- 0
      ps[3,i,t,3] <- 1
      
      # Define probabilities of O(t) given S(t)
      po[1,i,t,1] <- 0
      po[1,i,t,2] <- 1
      po[2,i,t,1] <- p[t]
      po[2,i,t,2] <- 1-p[t]
      po[3,i,t,1] <- 0
      po[3,i,t,2] <- 1
      } #t
   } #i

# Likelihood 
for (i in 1:M){
   # Define latent state at first occasion
   z[i,1] <- 1   # Make sure that all M individuals are in state 1 at t=1
   for (t in 2:n.occasions){
      # State process: draw S(t) given S(t-1)
      z[i,t] ~ dcat(ps[z[i,t-1], i, t-1,])
      # Observation process: draw O(t) given S(t)
      y[i,t] ~ dcat(po[z[i,t], i, t-1,])
      } #t
   } #i

# Calculate derived population parameters
for (t in 1:(n.occasions-1)){
   qgamma[t] <- 1-gamma[t]
   }
cprob[1] <- gamma[1]
for (t in 2:(n.occasions-1)){
   cprob[t] <- gamma[t] * prod(qgamma[1:(t-1)])
   } #t
psi <- sum(cprob[])            # Inclusion probability
for (t in 1:(n.occasions-1)){
   b[t] <- cprob[t] / psi      # Entry probability
   } #t

for (i in 1:M){
   for (t in 2:n.occasions){
      al[i,t-1] <- equals(z[i,t], 2)
      } #t
   for (t in 1:(n.occasions-1)){
      d[i,t] <- equals(z[i,t]-al[i,t],0)
      } #t   
   alive[i] <- sum(al[i,])
   } #i

for (t in 1:(n.occasions-1)){
   N[t] <- sum(al[,t])        # Actual population size
   B[t] <- sum(d[,t])         # Number of entries
   } #t
for (i in 1:M){
   w[i] <- 1-equals(alive[i],0)
   } #i
Nsuper <- sum(w[])            # Superpopulation size
}
",fill = TRUE)
sink()


# 10.3.3. The superpopulation parameterization
# Specify model in BUGS language
sink("js-super.jags")
cat("
model {
# Priors and constraints
for (i in 1:M){
   for (t in 1:(n.occasions-1)){
      phi[i,t] <- mean.phi
      } #t
   for (t in 1:n.occasions){
      p[i,t] <- mean.p
      } #t
   } #i

mean.phi ~ dunif(0, 1)         # Prior for mean survival
mean.p ~ dunif(0, 1)           # Prior for mean capture
psi ~ dunif(0, 1)              # Prior for inclusion probability

# Dirichlet prior for entry probabilities
for (t in 1:n.occasions){
   beta[t] ~ dgamma(1, 1)
   b[t] <- beta[t] / sum(beta[1:n.occasions])
   }

# Convert entry probs to conditional entry probs
nu[1] <- b[1]
for (t in 2:n.occasions){
   nu[t] <- b[t] / (1-sum(b[1:(t-1)]))
   } #t

# Likelihood
for (i in 1:M){
   # First occasion
   # State process
   w[i] ~ dbern(psi)                  # Draw latent inclusion
   z[i,1] ~ dbern(nu[1])
   # Observation process
   mu1[i] <- z[i,1] * p[i,1] * w[i]
   y[i,1] ~ dbern(mu1[i])

   # Subsequent occasions
   for (t in 2:n.occasions){
      # State process
      q[i,t-1] <- 1-z[i,t-1]
      mu2[i,t] <- phi[i,t-1] * z[i,t-1] + nu[t] * prod(q[i,1:(t-1)]) 
      z[i,t] ~ dbern(mu2[i,t])
      # Observation process
      mu3[i,t] <- z[i,t] * p[i,t] * w[i]
      y[i,t] ~ dbern(mu3[i,t])
      } #t
   } #i 

# Calculate derived population parameters
for (i in 1:M){
   for (t in 1:n.occasions){
      u[i,t] <- z[i,t]*w[i]     # Deflated latent state (u)
      }
   }
for (i in 1:M){
   recruit[i,1] <- u[i,1]
   for (t in 2:n.occasions){
      recruit[i,t] <- (1-u[i,t-1]) * u[i,t]
      } #t
   } #i
for (t in 1:n.occasions){
   N[t] <- sum(u[1:M,t])        # Actual population size
   B[t] <- sum(recruit[1:M,t])  # Number of entries
   } #t
for (i in 1:M){
   Nind[i] <- sum(u[i,1:n.occasions])
   Nalive[i] <- 1-equals(Nind[i], 0)
   } #i
Nsuper <- sum(Nalive[])         # Superpopulation size
}
",fill=TRUE)
sink()


# 10.4. Models with constant survival and time-dependent entry
# Define parameter values
n.occasions <- 7                         # Number of capture occasions
N <- 400                                 # Superpopulation size
phi <- rep(0.7, n.occasions-1)           # Survival probabilities
b <- c(0.34, rep(0.11, n.occasions-1))   # Entry probabilities 
p <- rep(0.5, n.occasions)               # Capture probabilities

PHI <- matrix(rep(phi, (n.occasions-1)*N), ncol = n.occasions-1, nrow = N, byrow = T)
P <- matrix(rep(p, n.occasions*N), ncol = n.occasions, nrow = N, byrow = T)

# Function to simulate capture-recapture data under the JS model
simul.js <- function(PHI, P, b, N){
   B <- rmultinom(1, N, b) # Generate no. of entering ind. per occasion
   n.occasions <- dim(PHI)[2] + 1
   CH.sur <- CH.p <- matrix(0, ncol = n.occasions, nrow = N)
   # Define a vector with the occasion of entering the population
   ent.occ <- numeric()
   for (t in 1:n.occasions){
      ent.occ <- c(ent.occ, rep(t, B[t]))
      }
   # Simulating survival
   for (i in 1:N){
      CH.sur[i, ent.occ[i]] <- 1   # Write 1 when ind. enters the pop.
      if (ent.occ[i] == n.occasions) next
      for (t in (ent.occ[i]+1):n.occasions){
         # Bernoulli trial: has individual survived occasion?
         sur <- rbinom(1, 1, PHI[i,t-1])
         ifelse (sur==1, CH.sur[i,t] <- 1, break)
         } #t
      } #i
   # Simulating capture
   for (i in 1:N){
      CH.p[i,] <- rbinom(n.occasions, 1, P[i,])
      } #i
   # Full capture-recapture matrix
   CH <- CH.sur * CH.p
   
   # Remove individuals never captured
   cap.sum <- rowSums(CH)
   never <- which(cap.sum == 0)
   CH <- CH[-never,]
   Nt <- colSums(CH.sur)    # Actual population size
   return(list(CH=CH, B=B, N=Nt))
   }

# Execute simulation function
sim <- simul.js(PHI, P, b, N)
CH <- sim$CH

# Augment the capture-histories by nz pseudo-individuals
nz <- 500
CH.aug <- rbind(CH, matrix(0, ncol = dim(CH)[2], nrow = nz))

# Bundle data
jags.data <- list(y = CH.aug, n.occasions = dim(CH.aug)[2], M = dim(CH.aug)[1])

# Initial values
# Good initial values for the latent state z are needed to run the model in JAGS. The simplest option that works is to give just a matrix with a 1 at all places.
z.init <- CH.aug
z.init[z.init==0] <- 1
inits <- function(){list(mean.phi = runif(1, 0, 1), mean.p = runif(1, 0, 1), z = z.init)}  

# Parameters monitored
parameters <- c("psi", "mean.p", "mean.phi", "b", "Nsuper", "N", "B", "gamma")

# MCMC settings
ni <- 5000
nt <- 3
nb <- 2000
nc <- 3

# Call JAGS from R (BRT 11 min)
js.occ <- jags(jags.data, inits, parameters, "js-rest.occ.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

print(js.occ, digits = 3)



# 10.4.2 Analysis of the JS model as a multistate model
# Add dummy occasion
CH.du <- cbind(rep(0, dim(CH)[1]), CH)

# Augment data
nz <- 500
CH.ms <- rbind(CH.du, matrix(0, ncol = dim(CH.du)[2], nrow = nz))

# Recode CH matrix: a 0 is not allowed in WinBUGS!
CH.ms[CH.ms==0] <- 2                     # Not seen = 2, seen = 1

# Bundle data
jags.data <- list(y = CH.ms, n.occasions = dim(CH.ms)[2], M = dim(CH.ms)[1])

# Initial values
# As always in JAGS, good initial values need to be specified for the latent state z. They need to correspond to the true state, which is not the same as the observed state. Thus, we have given initial values of "1" for the latend state at all places before an individual was observed, an initial value of "2" at all places when the individual was observed alive or known to be alive and an initial value of "3" at all places after the last observation. The folllowing function creates the initial values.

js.multistate.init <- function(ch, nz){
   ch[ch==2] <- NA
   state <- ch
   for (i in 1:nrow(ch)){
      n1 <- min(which(ch[i,]==1))
      n2 <- max(which(ch[i,]==1))
      state[i,n1:n2] <- 2
      }
   state[state==0] <- NA
   get.first <- function(x) min(which(!is.na(x)))
   get.last <- function(x) max(which(!is.na(x)))   
   f <- apply(state, 1, get.first)
   l <- apply(state, 1, get.last)
   for (i in 1:nrow(ch)){
      state[i,1:f[i]] <- 1
      if(l[i]!=ncol(ch)) state[i, (l[i]+1):ncol(ch)] <- 3
      state[i, f[i]] <- 2
      }   
   state <- rbind(state, matrix(1, ncol = ncol(ch), nrow = nz))
   return(state)
   }

inits <- function(){list(mean.phi = runif(1, 0, 1), mean.p = runif(1, 0, 1), z = js.multistate.init(CH.du, nz))}    

# Parameters monitored
parameters <- c("mean.p", "mean.phi", "b", "Nsuper", "N", "B")

# MCMC settings
ni <- 20000
nt <- 3
nb <- 5000
nc <- 3

# Call JAGS from R (BRT 32 min)
js.ms <- jags(jags.data, inits, parameters, "js-ms.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

print(js.ms, digits = 3)



# 10.4.3 Analysis of the JS model under the superpopulation parameterization
# Augment capture-histories by nz pseudo-individuals
nz <- 500
CH.aug <- rbind(CH, matrix(0, ncol = dim(CH)[2], nrow = nz))

# Bundle data
jags.data <- list(y = CH.aug, n.occasions = dim(CH.aug)[2], M = dim(CH.aug)[1])

# Initial values
# Good initial values for the latent state z are needed to run the model in JAGS. The simplest option that works is to give just a matrix with a 1 at all places. Moreover, initial values for the inclusion parameter w should also be given, the simplest option is to provide a vector with 1's.
z.init <- CH.aug
z.init[z.init==0] <- 1
w.init <- rep(1, nrow(CH.aug))
inits <- function(){list(mean.phi = runif(1, 0, 1), mean.p = runif(1, 0, 1), psi = runif(1, 0, 1), w = w.init, z = z.init)}  

# Parameters monitored
parameters <- c("psi", "mean.p", "mean.phi", "b", "Nsuper", "N", "B", "nu")

# MCMC settings
ni <- 5000
nt <- 3
nb <- 2000
nc <- 3

# Call JAGS from R (BRT 40 min)
js.super <- jags(jags.data, inits, parameters, "js-super.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

print(js.super, digits = 3)


# 10.4.4 Comparison of estimates
# Code to produce Fig. 10-4
par(mfrow = c(1,2), mar = c(5, 6, 2, 1), mgp=c(3.4, 1, 0), las = 1)
plot(density(js.occ$BUGSoutput$sims.list$Nsuper), main = "", xlab = "", ylab = "Density", frame = FALSE, lwd = 2, ylim=c(0, 0.023), col = "blue")
points(density(js.ms$BUGSoutput$sims.list$Nsuper), type = "l", lty = 2, col = "blue", lwd = 2)
points(density(js.super$BUGSoutput$sims.list$Nsuper), type = "l", lty = 3, col = "blue", lwd = 2)
abline(v = N, col = "red", lwd = 2)
mtext("Size of superpopulation", 1, line = 3)
text(x = 470, y = 0.02, "(a)")

b1.lower <- b2.lower <- b3.lower <- b1.upper <- b2.upper <- b3.upper <- numeric()
for (t in 1:n.occasions){
   b1.lower[t] <- quantile(js.occ$BUGSoutput$sims.list$b[,t], 0.025)
   b2.lower[t] <- quantile(js.ms$BUGSoutput$sims.list$b[,t], 0.025)
   b3.lower[t] <- quantile(js.super$BUGSoutput$sims.list$b[,t], 0.025)
   b1.upper[t] <- quantile(js.occ$BUGSoutput$sims.list$b[,t], 0.975)
   b2.upper[t] <- quantile(js.ms$BUGSoutput$sims.list$b[,t], 0.975)
   b3.upper[t] <- quantile(js.super$BUGSoutput$sims.list$b[,t], 0.975)
   }
time <- 1:n.occasions
plot(x = time-0.25, y = js.occ$BUGSoutput$mean$b, xlab = "", ylab = "Entry probability", frame = FALSE, las = 1, xlim = c(0.5, 7.5), pch = 16,  ylim = c(0, max(c(b1.upper, b2.upper))))
segments(time-0.25, b1.lower, time-0.25, b1.upper)
points(x = time, y = js.ms$BUGSoutput$mean$b, pch = 1)
segments(time, b2.lower, time, b2.upper)
points(x = time+0.25, y = js.super$BUGSoutput$mean$b, pch = 17)
segments(time+0.25, b3.lower, time+0.25, b3.upper)
points(x = time, y = b, pch = 18, col = "red")
mtext("Year", 1, line = 3)
text(x = 6, y = 0.39, "(b)")


# 10.5. Models with individual capture heterogeneity
# Define parameter values
n.occasions <- 8                         # Number of capture occasions
N <- 300                                 # Size of the superpopulation
phi <- rep(0.75, n.occasions-1)          # Survival probabilities
b <- c(0.37, rep(0.09, n.occasions-1))   # Entry probabilities 
mean.p <- 0.6                            # Mean capture probability
var.p <- 1                               # Indv. Variance of capture prob.
p <- plogis(rnorm(N, qlogis(mean.p), var.p^0.5))

PHI <- matrix(rep(phi, (n.occasions-1)*N), ncol = n.occasions-1, nrow = N, byrow = T)
P <- matrix(rep(p, n.occasions), ncol = n.occasions, nrow = N, byrow = F)

# Execute simulation function
sim <- simul.js(PHI, P, b, N)
CH <- sim$CH

# Specify model in BUGS language
sink("js-super-indran.jags")
cat("
model {
# Priors and constraints
for (i in 1:M){
   for (t in 1:(n.occasions-1)){
      phi[i,t] <- mean.phi
      } #t
   for (t in 1:n.occasions){
      logit(p[i,t]) <- mean.lp + epsilon[i]
      } #t
   } #i

mean.phi ~ dunif(0, 1)              # Prior for mean survival
mean.lp <- log(mean.p / (1-mean.p))
mean.p ~ dunif(0, 1)                # Prior for mean capture
for (i in 1:M){
   epsilon[i] ~ dnorm(0, tau)
   }
tau <- pow(sigma, -2)
sigma ~ dunif(0, 5)                  # Prior for sd of indv. variation of p
sigma2 <- pow(sigma, 2)
psi ~ dunif(0, 1)                    # Prior for inclusion probability

# Dirichlet prior for entry probabilities
for (t in 1:n.occasions){
   beta[t] ~ dgamma(1, 1)
   b[t] <- beta[t] / sum(beta[1:n.occasions])
   }

# Convert entry probs to conditional entry probs
nu[1] <- b[1]
for (t in 2:n.occasions){
   nu[t] <- b[t] / (1-sum(b[1:(t-1)]))
   } #t

# Likelihood
for (i in 1:M){
   # First occasion
   # State process
   w[i] ~ dbern(psi)                  # Draw latent inclusion
   z[i,1] ~ dbern(nu[1])
   # Observation process
   mu1[i] <- z[i,1] * p[i,1] * w[i]
   y[i,1] ~ dbern(mu1[i])

   # Subsequent occasions
   for (t in 2:n.occasions){
      # State process
      q[i,t-1] <- 1-z[i,t-1]
      mu2[i,t] <- phi[i,t-1] * z[i,t-1] + nu[t] * prod(q[i,1:(t-1)]) 
      z[i,t] ~ dbern(mu2[i,t])
      # Observation process
      mu3[i,t] <- z[i,t] * p[i,t] * w[i]
      y[i,t] ~ dbern(mu3[i,t])
      } #t
   } #i 

# Calculate derived population parameters
for (i in 1:M){
   for (t in 1:n.occasions){
      u[i,t] <- z[i,t]*w[i]     # Deflated latent state (u)
      }
   }
for (i in 1:M){
   recruit[i,1] <- u[i,1]
   for (t in 2:n.occasions){
      recruit[i,t] <- (1-u[i,t-1]) * u[i,t]
      } #t
   } #i
for (t in 1:n.occasions){
   N[t] <- sum(u[1:M,t])        # Actual population size
   B[t] <- sum(recruit[1:M,t])  # Number of entries
   } #t
for (i in 1:M){
   Nind[i] <- sum(u[i,1:n.occasions])
   Nalive[i] <- 1-equals(Nind[i], 0)
   } #i
Nsuper <- sum(Nalive[])         # Superpopulation size
}
",fill=TRUE)
sink()

# Augment the capture-histories by nz pseudo-individuals
nz <- 300
CH.aug <- rbind(CH, matrix(0, ncol = dim(CH)[2], nrow = nz))

# Bundle data
jags.data <- list(y = CH.aug, n.occasions = dim(CH.aug)[2], M = dim(CH.aug)[1])

# Initial values
# Good initial values for the latent state z are needed to run the model in JAGS. The simplest option that works is to give just a matrix with a 1 at all places. Moreover, initial values for the inclusion parameter w should also be given, the simplest option is to provide a vector with 1's.
z.init <- CH.aug
z.init[z.init==0] <- 1
w.init <- rep(1, nrow(CH.aug))
inits <- function(){list(mean.phi = runif(1, 0, 1), mean.p = runif(1, 0, 1), sigma = runif(1, 0, 1), w = w.init, z = z.init)}  

# Parameters monitored
parameters <- c("sigma2","psi", "mean.p", "mean.phi", "N", "Nsuper", "b", "B")

# MCMC settings
ni <- 20000
nt <- 6
nb <- 5000
nc <- 3

# Call JAGS from R (BRT 179 min)
js.ran <- jags(jags.data, inits, parameters, "js-super-indran.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

print(js.ran, digits = 3)


# 10.7. Analysis of a real data set: survival, recruitment and population size of Leisler?s bats
# Specify model in BUGS language
sink("js-tempran.jags")
cat("
model {
# Priors and constraints
for (i in 1:M){
   for (t in 1:(n.occasions-1)){
      logit(phi[i,t]) <- mean.lphi + epsilon[t]
      } #t
   for (t in 1:n.occasions){
      p[i,t] <- mean.p
      } #t
   } #i
mean.p ~ dunif(0, 1)                # Prior for mean capture
mean.phi ~ dunif(0, 1)              # Prior for mean survival
mean.lphi <- log(mean.phi / (1-mean.phi))
for (t in 1:(n.occasions-1)){
   epsilon[t] ~ dnorm(0, tau)
   }
tau <- pow(sigma, -2)
sigma ~ dunif(0, 5)              # Prior for sd of indv. variation of phi
sigma2 <- pow(sigma, 2)

for (t in 1:n.occasions){
   gamma[t] ~ dunif(0, 1)
   } #t

# Likelihood
for (i in 1:M){
   # First occasion
   # State process
   z[i,1] ~ dbern(gamma[1])
   mu1[i] <- z[i,1] * p[i,1]
   # Observation process
   y[i,1] ~ dbern(mu1[i])
   
   # Subsequent occasions
   for (t in 2:n.occasions){
      # State process
      q[i,t-1] <- 1-z[i,t-1]
      mu2[i,t] <- phi[i,t-1] * z[i,t-1] + gamma[t] * prod(q[i,1:(t-1)]) 
      z[i,t] ~ dbern(mu2[i,t])
      # Observation process
      mu3[i,t] <- z[i,t] * p[i,t]
      y[i,t] ~ dbern(mu3[i,t])
      } #t
   } #i

# Calculate derived population parameters
for (t in 1:n.occasions){
   qgamma[t] <- 1-gamma[t]
   }
cprob[1] <- gamma[1]
for (t in 2:n.occasions){
   cprob[t] <- gamma[t] * prod(qgamma[1:(t-1)])
   } #t
psi <- sum(cprob[])            # Inclusion probability
for (t in 1:n.occasions){
   b[t] <- cprob[t] / psi      # Entry probability
   } #t
for (i in 1:M){
   recruit[i,1] <- z[i,1]
   for (t in 2:n.occasions){
      recruit[i,t] <- (1-z[i,t-1]) * z[i,t]
      } #t
   } #i
for (t in 1:n.occasions){
   N[t] <- sum(z[1:M,t])        # Actual population size
   B[t] <- sum(recruit[1:M,t])  # Number of entries
   } #t
for (i in 1:M){
   Nind[i] <- sum(z[i,1:n.occasions])
   Nalive[i] <- 1-equals(Nind[i], 0)
   } #i
Nsuper <- sum(Nalive[])         # Size of superpopulation
}
",fill=TRUE)
sink()

leis <- as.matrix(read.table("leisleri.txt", sep = " ", header = FALSE))
nz <- 300
CH.aug <- rbind(leis, matrix(0, ncol = dim(leis)[2], nrow = nz))

# Bundle data
jags.data <- list(y = CH.aug, n.occasions = dim(CH.aug)[2], M = dim(CH.aug)[1])

# Initial values
# Good initial values for the latent state z are needed to run the model in JAGS. The simplest option that works is to give just a matrix with a 1 at all places.
z.init <- CH.aug
z.init[z.init==0] <- 1
inits <- function(){list(mean.phi = runif(1, 0, 1), mean.p = runif(1, 0, 1), sigma = runif(1, 0, 1), z = z.init)}  

# Parameters monitored
parameters <- c("psi", "mean.p", "sigma2", "mean.phi", "N", "Nsuper", "b", "B")

# MCMC settings
ni <- 10000
nt <- 6
nb <- 5000
nc <- 3

# Call JAGS from R (BRT 127 min)
nl <- jags(jags.data, inits, parameters, "js-tempran.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

print(nl, digits = 3)

# Code to produce Fig. 10-6
# Calculate per-capita recruitment
T <- dim(leis)[2]
f <- matrix(NA, ncol = T, nrow = length(nl$BUGSoutput$sims.list$B[,1]))
for (t in 1:(T-1)){
   f[,t] <- nl$BUGSoutput$sims.list$B[,t+1] / nl$BUGSoutput$sims.list$N[,t+1]
   }
n.lower <- n.upper <- f.lower <- f.upper <- f.mean <- numeric()
for (t in 1:T){
   n.lower[t] <- quantile(nl$BUGSoutput$sims.list$N[,t], 0.025)
   n.upper[t] <- quantile(nl$BUGSoutput$sims.list$N[,t], 0.975)
   }
for (t in 1:(T-1)){
   f.lower[t] <- quantile(f[,t], 0.025)
   f.upper[t] <- quantile(f[,t], 0.975)
   f.mean[t] <- mean(f[,t])
   }

par(mfrow = c(1, 2))
plot(nl$BUGSoutput$mean$N, type = "b", pch = 19, ylab = "Population size", xlab = "", axes = F, cex = 1.5, ylim = c(10, max(n.upper)))
axis(1, at = seq(1, T, 2), labels = seq(1990, 2008, 2))
axis(1, at = 1:T, labels = rep("", T), tcl = -0.25)
axis(2, las = 1)
segments(1:T, n.lower, 1:T, n.upper)

plot(f.mean, type = "b", pch = 19, ylab = "Local per capita recruitment", xlab = "", axes = F, cex = 1.5, ylim = c(0, 0.8))
axis(1, at = seq(1, (T-1), 2), labels = seq(1991, 2008, 2))
axis(1, at = 1:(T-1), labels = rep("", T-1), tcl = -0.25)
axis(2, las = 1)
segments(1:(T-1), f.lower, 1:(T-1), f.upper)




############################################################################
#
# 11. Integrated population models
# 
############################################################################


# 11.3. Example of a simple IPM (counts, capture-recapture, reproduction)
# 11.3.1. Load data
# Population counts (from years 1 to 10)
y <- c(45, 48, 44, 59, 62, 62, 55, 51, 46, 42)

# Capture-recapture data (in m-array format, from years 1 to 10)
m <- matrix(c(11,  0,  0,  0,  0,  0,  0,  0,  0,  70,
               0, 12,  0,  1,  0,  0,  0,  0,  0,  52,
               0,  0, 15,  5,  1,  0,  0,  0,  0,  42,
               0,  0,  0,  8,  3,  0,  0,  0,  0,  51,
               0,  0,  0,  0,  4,  3,  0,  0,  0,  61,
               0,  0,  0,  0,  0, 12,  2,  3,  0,  66,
               0,  0,  0,  0,  0,  0, 16,  5,  0,  44,
               0,  0,  0,  0,  0,  0,  0, 12,  0,  46,
               0,  0,  0,  0,  0,  0,  0,  0, 11,  71,
              10,  2,  0,  0,  0,  0,  0,  0,  0,  13,
               0,  7,  0,  1,  0,  0,  0,  0,  0,  27,
               0,  0, 13,  2,  1,  1,  0,  0,  0,  14,
               0,  0,  0, 12,  2,  0,  0,  0,  0,  20,
               0,  0,  0,  0, 10,  2,  0,  0,  0,  21,
               0,  0,  0,  0,  0, 11,  2,  1,  1,  14,
               0,  0,  0,  0,  0,  0, 12,  0,  0,  18,
               0,  0,  0,  0,  0,  0,  0, 11,  1,  21,
               0,  0,  0,  0,  0,  0,  0,  0, 10,  26), ncol = 10, byrow = TRUE)

# Productivity data (from years 1 to 9)
J <- c(64, 132,  86, 154, 156, 134, 116, 106, 110)
R <- c(21, 28, 26, 38, 35, 33, 31, 30, 33) 


# 11.3.2. Analysis of the model
# Specify model in BUGS language
sink("ipm.jags")
cat("
model {
#-------------------------------------------------
#  Integrated population model
#  - Age structured model with 2 age classes: 
#		1-year old and adult (at least 2 years old)
#  - Age at first breeding = 1 year
#  - Prebreeding census, female-based
#  - All vital rates assumed to be constant
#-------------------------------------------------

#-------------------------------------------------
# 1. Define the priors for the parameters
#-------------------------------------------------
# Observation error
tauy <- pow(sigma.y, -2)
sigma.y ~ dunif(0, 50)
sigma2.y <- pow(sigma.y, 2)

# Initial population sizes
# Note that this part is different than in BUGS:
# 1. JAGS seems to be very sensitive to the choice of the prior distribution for the initial population sizes and of the choice of the observation model (priors)
# 2. Since the initial population sizes are used in binomial distributions, the numbers must be integers, otherwise JAGS does not run.
# The following specification seems to work very well and produces similar results as BUGS:

n1 ~ dnorm(25, tauy)T(0,)     # 1-year
nad ~ dnorm(25, tauy)T(0,)    # Adults
N1[1] <- round(n1)
Nad[1] <- round(nad)

# Survival and recapture probabilities, as well as productivity
for (t in 1:(nyears-1)){
   sjuv[t] <- mean.sjuv
   sad[t] <- mean.sad
   p[t] <- mean.p
   f[t] <- mean.fec
   }

mean.sjuv ~ dunif(0, 1)
mean.sad ~ dunif(0, 1)
mean.p ~ dunif(0, 1)
mean.fec ~ dunif(0, 20)

#-------------------------------------------------
# 2. Derived parameters
#-------------------------------------------------
# Population growth rate
for (t in 1:(nyears-1)){
   lambda[t] <- Ntot[t+1] / Ntot[t]
   }

#-------------------------------------------------
# 3. The likelihoods of the single data sets
#-------------------------------------------------
# 3.1. Likelihood for population population count data (state-space model)
   # 3.1.1 System process
   for (t in 2:nyears){
      mean1[t] <- f[t-1] / 2 * sjuv[t-1] * Ntot[t-1]
      N1[t] ~ dpois(mean1[t])
      Nad[t] ~ dbin(sad[t-1], Ntot[t-1])
      }
   for (t in 1:nyears){
      Ntot[t] <- Nad[t] + N1[t]
      }
   
   # 3.1.2 Observation process
   for (t in 1:nyears){
      y[t] ~ dnorm(Ntot[t], tauy)
      }

# 3.2 Likelihood for capture-recapture data: CJS model (2 age classes)
# Multinomial likelihood
for (t in 1:2*(nyears-1)){
   m[t,1:nyears] ~ dmulti(pr[t,], r[t])
   }

# m-array cell probabilities for juveniles
for (t in 1:(nyears-1)){
   # Main diagonal
   q[t] <- 1-p[t]
   pr[t,t] <- sjuv[t] * p[t]
   # Above main diagonal
   for (j in (t+1):(nyears-1)){
      pr[t,j] <- sjuv[t]*prod(sad[(t+1):j])*prod(q[t:(j-1)])*p[j]
      } #j	
   # Below main diagonal
   for (j in 1:(t-1)){
      pr[t,j] <- 0
      } #j
   # Last column: probability of non-recapture
   pr[t,nyears] <- 1-sum(pr[t,1:(nyears-1)])
   } #t

# m-array cell probabilities for adults
for (t in 1:(nyears-1)){
   # Main diagonal
   pr[t+nyears-1,t] <- sad[t] * p[t]
   # Above main diagonal
   for (j in (t+1):(nyears-1)){
      pr[t+nyears-1,j] <- prod(sad[t:j])*prod(q[t:(j-1)])*p[j]
      } #j
   # Below main diagonal
   for (j in 1:(t-1)){
      pr[t+nyears-1,j] <- 0
      } #j
   # Last column
   pr[t+nyears-1,nyears] <- 1 - sum(pr[t+nyears-1,1:(nyears-1)])
   } #t

# 3.3. Likelihood for productivity data: Poisson regression
for (t in 1:(nyears-1)){
   J[t] ~ dpois(rho[t])
   rho[t] <- R[t]*f[t]
   }
}
",fill = TRUE)
sink()

# Bundle data
jags.data <- list(m = m, y = y, J = J, R = R, nyears = dim(m)[2], r = rowSums(m))

# Initial values
inits <- function(){list(mean.sjuv = runif(1, 0, 1), mean.sad = runif(1, 0, 1), mean.p = runif(1, 0, 1), mean.fec = runif(1, 0, 10), sigma.y = runif(1, 0, 1), n1 = rpois(1, 30), nad = rpois(1, 30))}  

# Parameters monitored
parameters <- c("mean.sjuv", "mean.sad", "mean.p", "mean.fec", "N1", "Nad", "Ntot", "sigma2.y", "lambda")

# MCMC settings
ni <- 2000
nt <- 6
nb <- 500
nc <- 3

# Call JAGS from R (BRT 2 min)
ipm <- jags(jags.data, inits, parameters, "ipm.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

print(ipm, digits = 3)

# Note that sometimes the model does not run, because of the estimation of the population growth rate (lambda). Presumably, awkward numbers such as 0 are produced during the initialising phase and then the ratio cannot be computed. There are two options: 1) just run the model again until it works, or 2) delete the calculation of lambda within the JAGS code and compute it afterwards based on the posterior samples of Ntot.


# Produce Fig. 11-4
par(cex = 1.2)
lower <- upper <- numeric()
for (i in 1:10){
   lower[i] <- quantile(ipm$BUGSoutput$sims.list$Ntot[,i], 0.025)
   upper[i] <- quantile(ipm$BUGSoutput$sims.list$Ntot[,i], 0.975)
   }
plot(ipm$BUGSoutput$mean$Ntot, type = "b", ylim = c(35, 65), ylab = "Population size", xlab = "Year", las = 1, pch = 16, col = "blue", frame = F, cex = 1.5)
segments(1:10, lower, 1:10, upper, col = "blue")
points(y, type = "b", col = "black", pch = 16, lty = 2, cex = 1.5)
legend(x = 1, y = 65, legend = c("Counts", "Estimates"), pch = c(16, 16), col = c("black", "blue"), lty = c(2, 1), bty = "n")
 



# 11.4. Another example of an IPM: Estimating productivity without explicit productivity data
# Specify model in BUGS language
sink("ipm-prod.jags")
cat("
model {
#-------------------------------------------------
#  Integrated population model
#  - Age structured model with 2 age classes: 
#		1-year old and adult (at least 2 years old)
#  - Age at first breeding = 1 year
#  - Prebreeding census, female-based
#  - All vital rates assumed to be constant
#-------------------------------------------------

#-------------------------------------------------
# 1. Define the priors for the parameters
#-------------------------------------------------
# Observation error
tauy <- pow(sigma.y, -2)
sigma.y ~ dunif(0, 50)
sigma2.y <- pow(sigma.y, 2)

# Initial population sizes
n1 ~ dnorm(25, tauy)T(0,)     # 1-year
nad ~ dnorm(25, tauy)T(0,)    # Adults
N1[1] <- round(n1)
Nad[1] <- round(nad)

# Survival and recapture probabilities, as well as productivity
for (t in 1:(nyears-1)){
   sjuv[t] <- mean.sjuv
   sad[t] <- mean.sad
   p[t] <- mean.p
   f[t] <- mean.fec
   }

mean.sjuv ~ dunif(0, 1)
mean.sad ~ dunif(0, 1)
mean.p ~ dunif(0, 1)
mean.fec ~ dunif(0, 20)

#-------------------------------------------------
# 2. Derived parameters
#-------------------------------------------------
# Population growth rate
for (t in 1:(nyears-1)){
   lambda[t] <- Ntot[t+1] / Ntot[t]
   }

#-------------------------------------------------
# 3. The likelihoods of the single data sets
#-------------------------------------------------
# 3.1. Likelihood for population population count data (state-space model)
   # 3.1.1 System process
   for (t in 2:nyears){
      mean1[t] <- f[t-1] / 2 * sjuv[t-1] * Ntot[t-1]
      N1[t] ~ dpois(mean1[t])
      Nad[t] ~ dbin(sad[t-1], Ntot[t-1])
      }
   for (t in 1:nyears){
      Ntot[t] <- Nad[t] + N1[t]
      }

   # 3.1.2 Observation process
   for (t in 1:nyears){
      y[t] ~ dnorm(Ntot[t], tauy)
      }

# 3.2 Likelihood for capture-recapture data: CJS model (2 age classes)
# Multinomial likelihood
for (t in 1:2*(nyears-1)){
   m[t,1:nyears] ~ dmulti(pr[t,], r[t])
   }

# m-array cell probabilities for juveniles
for (t in 1:(nyears-1)){
   # Main diagonal
   q[t] <- 1-p[t]
   pr[t,t] <- sjuv[t] * p[t]
   # Above main diagonal
   for (j in (t+1):(nyears-1)){
      pr[t,j] <- sjuv[t]*prod(sad[(t+1):j])*prod(q[t:(j-1)])*p[j]
      } #j	
   # Below main diagonal
   for (j in 1:(t-1)){
      pr[t,j] <- 0
      } #j
   # Last column: probability of non-recapture
   pr[t,nyears] <- 1-sum(pr[t,1:(nyears-1)])
   } #t

# m-array cell probabilities for adults
for (t in 1:(nyears-1)){
   # Main diagonal
   pr[t+nyears-1,t] <- sad[t] * p[t]
   # Above main diagonal
   for (j in (t+1):(nyears-1)){
      pr[t+nyears-1,j] <- prod(sad[t:j])*prod(q[t:(j-1)])*p[j]
      } #j
   # Below main diagonal
   for (j in 1:(t-1)){
      pr[t+nyears-1,j] <- 0
      } #j
   # Last column
   pr[t+nyears-1,nyears] <- 1 - sum(pr[t+nyears-1,1:(nyears-1)])
   } #t
}
",fill = TRUE)
sink()

# Bundle data
jags.data <- list(m = m, y = y, nyears = dim(m)[2], r = rowSums(m))

# Initial values
inits <- function(){list(mean.sjuv= runif(1, 0, 1), mean.sad = runif(1, 0, 1), mean.p = runif(1, 0, 1), mean.fec = runif(1, 0, 10), n1 = rpois(1, 30), nad = rpois(1, 30), sigma.y = runif(1, 0, 10))}  

# Parameters monitored
parameters <- c("mean.sjuv", "mean.sad", "mean.p", "mean.fec", "N1", "Nad", "Ntot", "sigma2.y", "lambda")

# MCMC settings
ni <- 20000
nt <- 6
nb <- 5000
nc <- 3

# Call JAGS from R (BRT 1 min)
ipm.prod <- jags(jags.data, inits, parameters, "ipm-prod.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

# Summarize posteriors
print(ipm.prod, digits = 3)




# 11.5. IPMs for population viability analysis
# Specify model in BUGS language
sink("ipm-pred.jags")
cat("
model {
#-------------------------------------------------
#  Integrated population model
#  - Age structured model with 2 age classes: 
#		1-year old and adult (at least 2 years old)
#  - Age at first breeding = 1 year
#  - Prebreeding census, female-based
#  - All vital rates assumed to be constant
#-------------------------------------------------

#-------------------------------------------------
# 1. Define the priors for the parameters
#-------------------------------------------------
# Observation error
tauy <- pow(sigma.y, -2)
sigma.y ~ dunif(0, 50)
sigma2.y <- pow(sigma.y, 2)

# Initial population sizes
n1 ~ dnorm(25, tauy)T(0,)     # 1-year
nad ~ dnorm(25, tauy)T(0,)    # Adults
N1[1] <- round(n1)
Nad[1] <- round(nad)

# Survival and recapture probabilities, as well as productivity
for (t in 1:(nyears-1+t.pred)){
   sjuv[t] <- mean.sjuv
   sad[t] <- mean.sad
   p[t] <- mean.p
   f[t] <- mean.fec
   }

mean.sjuv ~ dunif(0, 1)
mean.sad ~ dunif(0, 1)
mean.p ~ dunif(0, 1)
mean.fec ~ dunif(0, 20)

#-------------------------------------------------
# 2. Derived parameters
#-------------------------------------------------
# Population growth rate
for (t in 1:(nyears-1+t.pred)){
   lambda[t] <- Ntot[t+1] / Ntot[t]
   }

#-------------------------------------------------
# 3. The likelihoods of the single data sets
#-------------------------------------------------
# 3.1. Likelihood for population population count data (state-space model)
   # 3.1.1 System process
   for (t in 2:(nyears+t.pred)){
      mean1[t] <- f[t-1] / 2 * sjuv[t-1] * Ntot[t-1]
      N1[t] ~ dpois(mean1[t])
      Nad[t] ~ dbin(sad[t-1], Ntot[t-1])
      }
   for (t in 1:(nyears+t.pred)){
      Ntot[t] <- Nad[t] + N1[t]
      }

   # 3.1.2 Observation process
   for (t in 1:nyears){
      y[t] ~ dnorm(Ntot[t], tauy)
      }

# 3.2 Likelihood for capture-recapture data: CJS model (2 age classes)
# Multinomial likelihood
for (t in 1:2*(nyears-1)){
   m[t,1:nyears] ~ dmulti(pr[t,], r[t])
   }

# m-array cell probabilities for juveniles
for (t in 1:(nyears-1)){
   # Main diagonal
   q[t] <- 1-p[t]
   pr[t,t] <- sjuv[t] * p[t]
   # Above main diagonal
   for (j in (t+1):(nyears-1)){
      pr[t,j] <- sjuv[t]*prod(sad[(t+1):j])*prod(q[t:(j-1)])*p[j]
      } #j	
   # Below main diagonal
   for (j in 1:(t-1)){
      pr[t,j] <- 0
      } #j
   # Last column: probability of non-recapture
   pr[t,nyears] <- 1-sum(pr[t,1:(nyears-1)])
   } #t

# m-array cell probabilities for adults
for (t in 1:(nyears-1)){
   # Main diagonal
   pr[t+nyears-1,t] <- sad[t] * p[t]
   # Above main diagonal
   for (j in (t+1):(nyears-1)){
      pr[t+nyears-1,j] <- prod(sad[t:j])*prod(q[t:(j-1)])*p[j]
      } #j
   # Below main diagonal
   for (j in 1:(t-1)){
      pr[t+nyears-1,j] <- 0
      } #j
   # Last column
   pr[t+nyears-1,nyears] <- 1 - sum(pr[t+nyears-1,1:(nyears-1)])
   } #t

# 3.3. Likelihood for productivity data: Poisson regression
for (t in 1:(nyears-1)){
   J[t] ~ dpois(rho[t])
   rho[t] <- R[t]*f[t]
   }
}
",fill = TRUE)
sink()

# Give the number of future years for which population size shall be estimated
t.pred <- 5

# Bundle data
jags.data <- list(m = m, y = y, J = J, R = R, nyears = dim(m)[2], t.pred = t.pred, r = rowSums(m))

# Initial values
inits <- function(){list(mean.sjuv = runif(1, 0, 1), mean.sad = runif(1, 0, 1), mean.p = runif(1, 0, 1), mean.fec = runif(1, 0, 10), n1 = rpois(1, 30), nad = rpois(1, 30), sigma.y = runif(1, 0, 10))}  

# Parameters monitored
parameters <- c("mean.sjuv", "mean.sad", "mean.p", "mean.fec", "N1", "Nad", "Ntot", "sigma2.y", "lambda")

# MCMC settings
ni <- 20000
nt <- 6
nb <- 5000
nc <- 3

# Call JAGS from R (BRT 1 min)
ipm.pred <- jags(jags.data, inits, parameters, "ipm-pred.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

# Produce Fig. 11-5
par(cex = 1.2)
lower <- upper <- numeric()
for (i in 1:15){
   lower[i] <- quantile(ipm.pred$BUGSoutput$sims.list$Ntot[,i], 0.025)
   upper[i] <- quantile(ipm.pred$BUGSoutput$sims.list$Ntot[,i], 0.975)
   }
plot(ipm.pred$BUGSoutput$mean$Ntot, type = "b", ylim = c(10, max(upper)), ylab = "Population size", xlab = "Year", las = 1, pch = 16, col = "blue", frame = F)
segments(1:15, lower, 1:15, upper, col = "blue")
points(y, type = "b", col = "black", lty = 2, pch = 16)
legend(x = 1, y = 80, legend = c("Counts", "Estimates"), pch = c(16, 16), col = c("black", "blue"), lty = c(2, 1), bty = "n")
 
mean(ipm.pred$BUGSoutput$sims.list$Ntot[,15]<30)



# 11.6. Real data example: hoopoe population dynamics 
# Load data
nyears <- 9	  # Number of years

# Capture-recapture data: m-array of juveniles and adults (these are males and females together)
marray.j <- matrix (c(15, 3, 0, 0, 0, 0, 0, 0, 198, 0, 34, 9, 1, 0, 0, 0, 0, 287, 0, 0, 56, 8, 1, 0, 0, 0, 455, 0, 0, 0, 48, 3, 1, 0, 0, 518, 0, 0, 0, 0, 45, 13, 2, 0, 463, 0, 0, 0, 0, 0, 27, 7, 0, 493, 0, 0, 0, 0, 0, 0, 37, 3, 434, 0, 0, 0, 0, 0, 0, 0, 39, 405), nrow = 8, ncol = 9, byrow = TRUE)
marray.a <- matrix(c(14, 2, 0, 0, 0, 0, 0, 0, 43, 0, 22, 4, 0, 0, 0, 0, 0, 44, 0, 0, 34, 2, 0, 0, 0, 0, 79, 0, 0, 0, 51, 3, 0, 0, 0, 94, 0, 0, 0, 0, 45, 3, 0, 0, 118, 0, 0, 0, 0, 0, 44, 3, 0, 113, 0, 0, 0, 0, 0, 0, 48, 2, 99, 0, 0, 0, 0, 0, 0, 0, 51, 90), nrow = 8, ncol = 9, byrow = TRUE)

# Population count data
popcount <- c(32, 42, 64, 85, 82, 78, 73, 69, 79)

# Productivity data
J <- c(189, 274, 398, 538, 520, 476, 463, 438) # Number of offspring 
R <- c(28, 36, 57, 77, 81, 83, 77, 72)         # Number of surveyed broods 

# Specify model in BUGS language
sink("ipm.hoopoe.jags")
cat("
model {
#------------------------------------------------------------
#  Integrated population model
#  - Age structured model with 2 age classes: 
#		1-year old and adult (at least 2 years old)
#  - Age at first breeding = 1 year
#  - Prebreeding census, female-based
#  - All vital rates are assumed to be time-dependent (random)
#  - Explicit estimation of immigration
#-------------------------------------------------------------

#----------------------------------------
# 1. Define the priors for the parameters
#----------------------------------------
# Initial population sizes
n1 ~ dnorm(100, 0.0001)I(0,)           # 1-year old individuals
nadSurv ~ dnorm(100, 0.0001)I(0,)      # Adults >= 2 years
nadimm ~ dnorm(100, 0.0001)I(0,)       # Immigrants
N1[1] <- round(n1)
NadSurv[1] <- round(nadSurv)
Nadimm[1] <- round(nadimm)

# Mean demographic parameters (on appropriate scale)
l.mphij ~ dnorm(0, 0.0001)
l.mphia ~ dnorm(0, 0.0001)
l.mfec ~ dnorm(0, 0.0001)
l.mim ~ dnorm(0, 0.0001)
l.p ~ dnorm(0, 0.0001)

# Precision of standard deviations of temporal variability
sig.phij ~ dunif(0, 10)
tau.phij <- pow(sig.phij, -2)
sig.phia ~ dunif(0, 10)
tau.phia <- pow(sig.phia, -2)
sig.fec ~ dunif(0, 10)
tau.fec <- pow(sig.fec, -2)
sig.im ~ dunif(0, 10)
tau.im <- pow(sig.im, -2)

# Distribution of error terms (Bounded to help with convergence)
for (t in 1:(nyears-1)){
   epsilon.phij[t] ~ dnorm(0, tau.phij)	
   epsilon.phia[t] ~ dnorm(0, tau.phia)
   epsilon.fec[t] ~ dnorm(0, tau.fec)
   epsilon.im[t] ~ dnorm(0, tau.im)
   }

#-------------------------
# 2. Constrain parameters
#-------------------------
for (t in 1:(nyears-1)){
   logit(phij[t]) <- l.mphij + epsilon.phij[t]  # Juv. apparent survival
   logit(phia[t]) <- l.mphia + epsilon.phia[t]  # Adult apparent survival
   log(f[t]) <- l.mfec + epsilon.fec[t]         # Productivity
   log(omega[t]) <- l.mim + epsilon.im[t]       # Immigration
   logit(p[t]) <- l.p                           # Recapture probability
   }

#-----------------------
# 3. Derived parameters
#-----------------------
mphij <- exp(l.mphij)/(1+exp(l.mphij))   # Mean juvenile survival probability
mphia <- exp(l.mphia)/(1+exp(l.mphia))   # Mean adult survival probability
mfec <- exp(l.mfec)                      # Mean productivity
mim <- exp(l.mim)                        # Mean immigration rate

# Population growth rate
for (t in 1:(nyears-1)){
   lambda[t] <- Ntot[t+1] / Ntot[t]
   logla[t] <- log(lambda[t])
   }
mlam <- exp((1/(nyears-1))*sum(logla[1:(nyears-1)]))   # Geometric mean

#--------------------------------------------
# 4. The likelihoods of the single data sets
#--------------------------------------------
# 4.1. Likelihood for population population count data (state-space model)
   # 4.1.1 System process
   for (t in 2:nyears){
      mean1[t] <- 0.5 * f[t-1] * phij[t-1] * Ntot[t-1]
      N1[t] ~ dpois(mean1[t])
      NadSurv[t] ~ dbin(phia[t-1], Ntot[t-1])
      mpo[t] <- Ntot[t-1] * omega[t-1]
      Nadimm[t] ~ dpois(mpo[t])
      }

   # 4.1.2 Observation process
   for (t in 1:nyears){
      Ntot[t] <- NadSurv[t] + Nadimm[t] + N1[t]
      y[t] ~ dpois(Ntot[t])
      }

# 4.2 Likelihood for capture-recapture data: CJS model (2 age classes)
# Multinomial likelihood
for (t in 1:(nyears-1)){
   marray.j[t,1:nyears] ~ dmulti(pr.j[t,], r.j[t])
   marray.a[t,1:nyears] ~ dmulti(pr.a[t,], r.a[t])
   }

# m-array cell probabilities for juveniles
for (t in 1:(nyears-1)){
   q[t] <- 1-p[t]
   # Main diagonal
   pr.j[t,t] <- phij[t]*p[t]
   # Above main diagonal
   for (j in (t+1):(nyears-1)){
      pr.j[t,j] <- phij[t]*prod(phia[(t+1):j])*prod(q[t:(j-1)])*p[j]
      } #j
   # Below main diagonal
   for (j in 1:(t-1)){
      pr.j[t,j] <- 0
      } #j
   # Last column
   pr.j[t,nyears] <- 1-sum(pr.j[t,1:(nyears-1)])
   } #t

# m-array cell probabilities for adults
for (t in 1:(nyears-1)){
   # Main diagonal
   pr.a[t,t] <- phia[t]*p[t]
   # above main diagonal
   for (j in (t+1):(nyears-1)){
      pr.a[t,j] <- prod(phia[t:j])*prod(q[t:(j-1)])*p[j]
      } #j
   # Below main diagonal
   for (j in 1:(t-1)){
      pr.a[t,j] <- 0
      } #j
   # Last column
   pr.a[t,nyears] <- 1-sum(pr.a[t,1:(nyears-1)])
   } #t

# 4.3. Likelihood for productivity data: Poisson regression
for (t in 1:(nyears-1)){
   J[t] ~ dpois(rho[t])
   rho[t] <- R[t] * f[t]
   }
}
",fill = TRUE)
sink()

# Bundle data
jags.data <- list(nyears = nyears, marray.j = marray.j, marray.a = marray.a, y = popcount, J = J, R = R, r.j = rowSums(marray.j), r.a = rowSums(marray.a))

# Initial values
inits <- function(){list(l.mphij = rnorm(1, 0.2, 0.5), l.mphia = rnorm(1, 0.2, 0.5), l.mfec = rnorm(1, 0.2, 0.5), l.mim = rnorm(1, 0.2, 0.5), l.p = rnorm(1, 0.2, 1), sig.phij = runif(1, 0.1, 10), sig.phia = runif(1, 0.1, 10), sig.fec = runif(1, 0.1, 10), sig.im = runif(1, 0.1, 10), n1 = round(runif(1, 1, 50), 0), nadSurv = round(runif(1, 5, 50), 0), nadimm = round(runif(1, 1, 50), 0))}

# Parameters monitored
parameters <- c("phij", "phia", "f", "omega", "p", "lambda", "mphij", "mphia", "mfec", "mim", "mlam", "sig.phij", "sig.phia", "sig.fec", "sig.im", "N1", "NadSurv", "Nadimm", "Ntot")

# MCMC settings
ni <- 20000
nt <- 6
nb <- 10000
nc <- 3

# Call JAGS from R (BRT 5 min)
ipm.hoopoe <- jags(jags.data, inits, parameters, "ipm.hoopoe.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())


# Code for Fig. 11-7
par(mfrow = c(2, 2), cex.axis = 1.2, cex.lab = 1.2, mar = c(5, 6, 1.5, 2), las = 1)
lower <- upper <- numeric()
year <- 2002:2010
for (i in 1:nyears){
   lower[i] <- quantile(ipm.hoopoe$BUGSoutput$sims.list$Ntot[,i], 0.025)
   upper[i] <- quantile(ipm.hoopoe$BUGSoutput$sims.list$Ntot[,i], 0.975)}
m1 <- min(c(ipm.hoopoe$BUGSoutput$mean$Ntot, popcount, lower), na.rm = T)
m2 <- max(c(ipm.hoopoe$BUGSoutput$mean$Ntot, popcount, upper), na.rm = T)
plot(0, 0, ylim = c(0, m2), xlim = c(1, nyears), ylab = "Population size", xlab = " ", col = "black", type = "l",  axes = F, frame = F)
axis(2)
axis(1, at = 1:nyears, labels = year)
polygon(x = c(1:nyears, nyears:1), y = c(lower, upper[nyears:1]), col = "grey90", border = "grey90")
points(popcount, type = "l", col = "grey30", lwd = 2)
points(ipm.hoopoe$BUGSoutput$mean$Ntot, type = "l", col = "blue", lwd = 2)
legend(x = 2, y = 25, legend = c("Counts", "Estimates"), lty = c(1, 1),lwd = c(2, 2), col = c("grey30", "blue"), bty = "n", cex = 1)

lower <- upper <- numeric()
T <- nyears-1
for (t in 1:T){
   lower[t] <- quantile(ipm.hoopoe$BUGSoutput$sims.list$phij[,t], 0.025)
   upper[t] <- quantile(ipm.hoopoe$BUGSoutput$sims.list$phij[,t], 0.975)}
par(mgp=c(3.8,1,0))
plot(y = ipm.hoopoe$BUGSoutput$mean$phij, x = (1:T)+0.5, xlim= c(1, 9), type = "b", pch = 16, ylim = c(0, 0.6), ylab = "Annual survival probability", xlab = "", axes = F, cex = 1.5, frame = F, lwd = 2)
axis(2)
axis(1, at = 1:(T+1), labels = 2002:2010)
segments((1:T)+0.5, lower, (1:T)+0.5, upper, lwd = 2)
segments(1, ipm.hoopoe$BUGSoutput$mean$mphij, T+1, ipm.hoopoe$BUGSoutput$mean$mphij, lty = 2, lwd = 2, col = "red")
segments(1, quantile(ipm.hoopoe$BUGSoutput$sims.list$mphij, 0.025), T+1, quantile(ipm.hoopoe$BUGSoutput$sims.list$mphij, 0.025), lty = 2, col = "red") 
segments(1, quantile(ipm.hoopoe$BUGSoutput$sims.list$mphij, 0.975), T+1, quantile(ipm.hoopoe$BUGSoutput$sims.list$mphij, 0.975), lty = 2, col = "red")
for (t in 1:T){
   lower[t] <- quantile(ipm.hoopoe$BUGSoutput$sims.list$phia[,t], 0.025)
   upper[t] <- quantile(ipm.hoopoe$BUGSoutput$sims.list$phia[,t], 0.975)}
points(y=ipm.hoopoe$BUGSoutput$mean$phia, x = (1:T)+0.5, type = "b", pch = 1, cex = 1.5, lwd = 2)
segments((1:T)+0.5, lower, (1:T)+0.5, upper, lwd = 2)
segments(1, ipm.hoopoe$BUGSoutput$mean$mphia, T+1, ipm.hoopoe$BUGSoutput$mean$mphia, lty = 2, lwd = 2, col = "red")
segments(1, quantile(ipm.hoopoe$BUGSoutput$sims.list$mphia, 0.025), T+1, quantile(ipm.hoopoe$BUGSoutput$sims.list$mphia, 0.025), lty = 2, col = "red") 
segments(1, quantile(ipm.hoopoe$BUGSoutput$sims.list$mphia, 0.975), T+1, quantile(ipm.hoopoe$BUGSoutput$sims.list$mphia, 0.975), lty = 2, col = "red")
legend(x = 4.5, y = 0.66, legend = c("Adults", "Juveniles"), pch = c(1, 16), bty = "n")

lower <- upper <- numeric()
T <- nyears-1
for (t in 1:T){
   lower[t] <- quantile(ipm.hoopoe$BUGSoutput$sims.list$f[,t], 0.025)
   upper[t] <- quantile(ipm.hoopoe$BUGSoutput$sims.list$f[,t], 0.975)}
plot(y=ipm.hoopoe$BUGSoutput$mean$f, x = (1:T), type = "b", pch = 16, ylim = c(5, 9), ylab = "Fecundity (fledgling / female)", xlab = "", axes = F, cex = 1.5, frame = F, lwd = 2)
axis(2)
axis(1, at = 1:T, labels = 2003:2010)
segments((1:T), lower, (1:T), upper)
segments(1, ipm.hoopoe$BUGSoutput$mean$mfec, T, ipm.hoopoe$BUGSoutput$mean$mfec, lty = 2, lwd = 2, col = "red")
segments(1, quantile(ipm.hoopoe$BUGSoutput$sims.list$mfec, 0.025), T, quantile(ipm.hoopoe$BUGSoutput$sims.list$mfec, 0.025), lty = 2, col = "red") 
segments(1, quantile(ipm.hoopoe$BUGSoutput$sims.list$mfec, 0.975), T, quantile(ipm.hoopoe$BUGSoutput$sims.list$mfec, 0.975), lty = 2, col = "red")

lower <- upper <- numeric()
T <- nyears-1
for (t in 1:T){
   lower[t] <- quantile(ipm.hoopoe$BUGSoutput$sims.list$omega[,t], 0.025)
   upper[t] <- quantile(ipm.hoopoe$BUGSoutput$sims.list$omega[,t], 0.975)}
plot(y = ipm.hoopoe$BUGSoutput$mean$omega, x = (1:T)+0.5, xlim = c(1, 9), type = "b", pch = 16, ylim = c(0, 1.1), ylab = "Immigration rate", xlab = "", axes = F, cex = 1.5, frame = F, lwd = 2)
axis(2)
axis(1, at = 1:(T+1), labels = 2002:2010)
segments((1:T)+0.5, lower, (1:T)+0.5, upper)
segments(1, ipm.hoopoe$BUGSoutput$mean$mim, T+1, ipm.hoopoe$BUGSoutput$mean$mim, lty = 2, lwd = 2, col = "red")
segments(1, quantile(ipm.hoopoe$BUGSoutput$sims.list$mim, 0.025), T+1, quantile(ipm.hoopoe$BUGSoutput$sims.list$mim, 0.025), lty = 2, col = "red") 
segments(1, quantile(ipm.hoopoe$BUGSoutput$sims.list$mim, 0.975), T+1, quantile(ipm.hoopoe$BUGSoutput$sims.list$mim, 0.975), lty = 2, col = "red")


# Code for Fig. 11-8 and some descriptive statistics
nyears <- 9
lambda.h <- lam.lower.h <- lam.upper.h <- numeric()
Fitted.h <- lower.h <- upper.h <- matrix(NA, nrow = nyears-1, ncol = 4)

for (i in 1:(nyears-1)){
   lambda.h[i] <- mean(ipm.hoopoe$BUGSoutput$sims.list$lambda[,i])
   lam.lower.h[i] <- quantile(ipm.hoopoe$BUGSoutput$sims.list$lambda[,i], 0.025)
   lam.upper.h[i] <- quantile(ipm.hoopoe$BUGSoutput$sims.list$lambda[,i], 0.975)
   }

for (i in 1:(nyears-1)){
   Fitted.h[i,1] <- mean(ipm.hoopoe$BUGSoutput$sims.list$phij[,i])
   lower.h[i,1] <- quantile(ipm.hoopoe$BUGSoutput$sims.list$phij[,i], 0.025)
   upper.h[i,1] <- quantile(ipm.hoopoe$BUGSoutput$sims.list$phij[,i], 0.975)
   }

for (i in 1:(nyears-1)){
   Fitted.h[i,2] <- mean(ipm.hoopoe$BUGSoutput$sims.list$phia[,i])
   lower.h[i,2] <- quantile(ipm.hoopoe$BUGSoutput$sims.list$phia[,i], 0.025)
   upper.h[i,2] <- quantile(ipm.hoopoe$BUGSoutput$sims.list$phia[,i], 0.975)
   }

for (i in 1:(nyears-1)){
   Fitted.h[i,3] <- mean(ipm.hoopoe$BUGSoutput$sims.list$f[,i])
   lower.h[i,3] <- quantile(ipm.hoopoe$BUGSoutput$sims.list$f[,i], 0.025)
   upper.h[i,3] <- quantile(ipm.hoopoe$BUGSoutput$sims.list$f[,i], 0.975)
   }

for (i in 1:(nyears-1)){
   Fitted.h[i,4] <- mean(ipm.hoopoe$BUGSoutput$sims.list$omega[,i])
   lower.h[i,4] <- quantile(ipm.hoopoe$BUGSoutput$sims.list$omega[,i], 0.025)
   upper.h[i,4] <- quantile(ipm.hoopoe$BUGSoutput$sims.list$omega[,i], 0.975)
   }


# Calculate some correlation coefficients
correl.h <- matrix(NA, ncol = 4, nrow = 5001)
for (i in 1:5001){
   correl.h[i,1] <- cor(ipm.hoopoe$BUGSoutput$sims.list$lambda[i,], ipm.hoopoe$BUGSoutput$sims.list$phij[i,])
   correl.h[i,2] <- cor(ipm.hoopoe$BUGSoutput$sims.list$lambda[i,], ipm.hoopoe$BUGSoutput$sims.list$phia[i,])
   correl.h[i,3] <- cor(ipm.hoopoe$BUGSoutput$sims.list$lambda[i,], ipm.hoopoe$BUGSoutput$sims.list$f[i,])
   correl.h[i,4] <- cor(ipm.hoopoe$BUGSoutput$sims.list$lambda[i,], ipm.hoopoe$BUGSoutput$sims.list$omega[i,])
   }

# Credible intervals of correlation coefficients
quantile(correl.h[,1], c(0.05, 0.5, 0.95), na.rm = TRUE)
quantile(correl.h[,2], c(0.05, 0.5, 0.95), na.rm = TRUE)
quantile(correl.h[,3], c(0.05, 0.5, 0.95), na.rm = TRUE)
quantile(correl.h[,4], c(0.05, 0.5, 0.95), na.rm = TRUE)


# Compute the posterior modes of correlation coefficients
m <- density(correl.h[,1], na.rm = TRUE)
m$x[which(m$y==max(m$y))]

m <- density(correl.h[,2], na.rm = TRUE)
m$x[which(m$y==max(m$y))]

m <- density(correl.h[,3], na.rm = TRUE)
m$x[which(m$y==max(m$y))]

m <- density(correl.h[,4], na.rm = TRUE)
m$x[which(m$y==max(m$y))]


# Probability that correlation coefficients (r) > 0
sum(correl.h[!is.na(correl.h[,1]),1]>0)/5001
sum(correl.h[!is.na(correl.h[,2]),2]>0)/5001
sum(correl.h[!is.na(correl.h[,3]),3]>0)/5001
sum(correl.h[!is.na(correl.h[,4]),4]>0)/5001


# Plot Fig. 11-8
par(mfrow = c(2, 2), mar = c(5.5, 4, 1.5, 1), mgp=c(3, 1, 0), las = 1, cex = 1)
linecol <- c("grey70")
plot(y = lambda.h, Fitted.h[,1], type = "n", xlim = c(0.05, 0.25), ylim = c(0.6, 1.8), ylab = "Population growth rate", xlab = "Juvenile survival", frame = FALSE, pch = 19)
segments(Fitted.h[,1], lam.lower.h, Fitted.h[,1], lam.upper.h, col = linecol)
segments(lower.h[,1], lambda.h, upper.h[,1], lambda.h, col = linecol)
points(y = lambda.h, Fitted.h[,1], pch = 19, col = "black")
text(x = 0.13, y = 0.75, "r = 0.75 (0.11, 0.90)", pos = 4, font = 3, cex = 0.8)
text(x = 0.13, y = 0.65, "P(r>0) = 0.97", pos = 4, font = 3, cex = 0.8)

par(mar = c( 5.5, 4, 1.5, 1))
plot(y = lambda.h, Fitted.h[,2], type = "n", xlim = c(0.3, 0.51), ylim = c(0.6, 1.8),  ylab = "", xlab = "Adult survival", frame.plot = FALSE, pch = 19)
segments(Fitted.h[,2], lam.lower.h, Fitted.h[,2], lam.upper.h, col = linecol)
segments(lower.h[,2], lambda.h, upper.h[,2], lambda.h, col = linecol)
points(y = lambda.h, Fitted.h[,2], pch = 19, col = "black")
text(x = 0.295, y = 1.75, "r = 0.26 (-0.44, 0.73)", pos = 4, font = 3, cex = 0.8)
text(x = 0.295, y = 1.65, "P(r>0) = 0.71", pos = 4, font = 3, cex = 0.8)

par(mar = c(5, 4, 2, 1))
plot(y = lambda.h, Fitted.h[,3], type = "n", xlim = c(5, 8.5), ylim = c(0.6, 1.8), ylab = "Population growth rate", xlab = "Fecundity", frame.plot = FALSE, pch = 19)
segments(Fitted.h[,3], lam.lower.h, Fitted.h[,3], lam.upper.h, col = linecol)
segments(lower.h[,3], lambda.h, upper.h[,3], lambda.h, col = linecol)
points(y=lambda.h, Fitted.h[,3], pch = 19, col = "black")
text(x = 5, y = 1.75, "r = 0.70 (0.21, 0.87)", pos = 4,  font = 3, cex = 0.8)
text(x = 5, y = 1.65, "P(r>0) = 0.99", pos = 4, font = 3, cex = 0.8)

par(mar = c( 5, 4, 2, 1))
plot(y = lambda.h, Fitted.h[,4], type = "n", xlim = c(0, 0.8), ylim = c(0.6, 1.8),  ylab = "", xlab = "Immigration rate", frame.plot = FALSE, pch = 19)
segments(Fitted.h[,4], lam.lower.h, Fitted.h[,4], lam.upper.h, col = linecol)
segments(lower.h[,4], lambda.h, upper.h[,4], lambda.h, col = linecol)
points(y=lambda.h, Fitted.h[,4], pch = 19, col = "black")
text(x = 0.35, y = 0.8, "r = 0.83 (-0.29, 0.94)", pos = 4, font = 3, cex = 0.8)
text(x = 0.35, y = 0.7, "P(r>0) = 0.86", pos = 4, font = 3, cex = 0.8)





######################################################################################################
# 
# 12. Metapopulation modeling of abundance using hierarchical Poisson regression: binomial mixture models
# 
#######################################################################################################

# 12.2. Generation and analysis of simulated data
# 12.2.1. The simplest case with constant parameters
# Determine sample sizes (spatial and temporal replication)
R <- 200
T <- 3

# Create structure to contain counts
y <- array(dim = c(R, T))

# Sample abundance from a Poisson(lambda = 2)
N <- rpois(n = R, lambda = 2)

# Sample counts from a Binomial(N, p = 0.5)
for (j in 1:T){
   y[,j] <- rbinom(n = R, size = N, prob = 0.5)
   }

# Look at realization of biological and observation processes
cbind(N, y)

# Specify model in BUGS language
sink("model.jags")
cat("
model {

# Priors
lambda ~ dgamma(0.005, 0.005) # Standard vague prior for lambda
# lambda ~ dunif(0, 10)       # Other possibility
p ~ dunif(0, 1)

# Likelihood
# Biological model for true abundance
for (i in 1:R) {
   N[i] ~ dpois(lambda)
   # Observation model for replicated counts
   for (j in 1:T) {
      y[i,j] ~ dbin(p, N[i])
      } # j
   } # i
}
",fill = TRUE)
sink()

# Bundle data
win.data <- list(y = y, R = nrow(y), T = ncol(y))

# Initial values
Nst <- apply(y, 1, max) + 1	# This line is important
inits <- function() list(N = Nst)

# Parameters monitored
params <- c("lambda", "p")

# MCMC settings
ni <- 1200
nt <- 2
nb <- 200
nc <- 3

# Call JAGS from R (BRT 0.1 min)
out <- jags(win.data, inits, params, "model.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

# Summarize posteriors
print(out, dig = 2)


# 12.2.2. Introducing covariates
# Define function for generating binomial-mix model data
data.fn <- function(R = 200, T = 3, xmin = -1, xmax = 1, alpha0 = 1, alpha1 = 3, beta0 = 0, beta1 = -5){

# R: number of sites at which counts were made (= number of spatial reps)
# T: number of times that counts were made at each site 
# (= number of temporal reps)
# xmin, xmax: define range of the covariate X
# alpha0 and alpha1: intercept and slope of log-linear regression 
# relating abundance to the site covariate A
# beta0 and beta1: intercept and slope of logistic-linear regression 
# of detection probability on A

   y <- array(dim = c(R, T))	# Array for counts

   # Ecological process
   # Covariate values: sort for ease of presentation
   X <- sort(runif(n = R, min = xmin, max = xmax))

   # Relationship expected abundance - covariate
   lam <- exp(alpha0 + alpha1 * X)

   # Add Poisson noise: draw N from Poisson(lambda)
   N <- rpois(n = R, lambda = lam)
   table(N)                # Distribution of abundances across sites
   sum(N > 0) / R          # Empirical occupancy
   totalN <- sum(N)  ;  totalN

   # Observation process
   # Relationship detection prob - covariate
   p <- plogis(beta0 + beta1 * X)

   # Make a ?census? (i.e., go out and count things)
   for (i in 1:T){
      y[,i] <- rbinom(n = R, size = N, prob = p)
      }

   # Na?ve regression
   naive.pred <- exp(predict(glm(apply(y, 1, max) ~ X + I(X^2), family =   poisson)))

   # Plot features of the simulated system
   par(mfrow = c(2, 2))
   plot(X, lam, main = "Expected abundance", xlab = "Covariate", ylab = "lambda", las = 1, type = "l", col = "red", lwd = 3, frame.plot = FALSE)
   plot(X, N, main = "Realised abundance", xlab = "Covariate", ylab = "N", las = 1, frame.plot = FALSE, col = "red", cex = 1.2)
   plot(X, p, ylim = c(0, 1), main = "Detection probability", xlab = "Covariate", ylab = "p", type = "l", col = "red", lwd = 3, las = 1, frame.plot = FALSE)
   plot(X, naive.pred, main = "Actual counts \n and na?ve regression", xlab = "Covariate", ylab = "Relative abundance", ylim = c(min(y), max(y)), type = "l", lty = 2, lwd = 4, col = "blue", las = 1, frame.plot = FALSE)
points(rep(X, T), y, col = "black", cex = 1.2)

   # Return stuff
   return(list(R = R, T = T, X = X, alpha0 = alpha0, alpha1 = alpha1, beta0 = beta0, beta1 = beta1, lam = lam, N = N, totalN = totalN, p = p, y = y))
   }

data <- data.fn()
str(data)

# Specify model in BUGS language
sink("model.jags")
cat("
model {

# Priors
alpha0 ~ dunif(-10, 10)
alpha1 ~ dunif(-10, 10)
beta0 ~ dunif(-10, 10)
beta1 ~ dunif(-10, 10)

# Likelihood
# Ecological model for true abundance
for (i in 1:R){
   N[i] ~ dpois(lambda[i])
   log(lambda[i]) <- alpha0 + alpha1 * X[i]

   # Observation model for replicated counts
   for (j in 1:T){
      y[i,j] ~ dbin(p[i,j], N[i])
      p[i,j] <- exp(lp[i,j])/(1+exp(lp[i,j]))
      lp[i,j] <- beta0 + beta1 * X[i]
      } #j
   } #i

# Derived quantities
totalN <- sum(N[])
}
",fill = TRUE)
sink()

# Bundle data
y <- data$y
win.data <- list(y = y, R = nrow(y), T = ncol(y), X = data$X)

# Initial values
Nst <- apply(y, 1, max) + 1	# Important to give good inits for latent N
inits <- function() list(N = Nst, alpha0 = runif(1, -1, 1), alpha1 = runif(1, -1, 1), beta0 = runif(1, -1, 1), beta1 = runif(1, -1, 1))

# Parameters monitored
params <- c("totalN", "alpha0", "alpha1", "beta0", "beta1")

# MCMC settings
ni <- 22000
nt <- 20
nb <- 2000
nc <- 3

# Call JAGS from R (BRT 4 min)
out <- jags(win.data, inits, params, "model.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

# Summarize posteriors
print(out, dig = 3)

# Plot posteriors
par(mfrow = c(2, 2))
hist(out$BUGSoutput$sims.list$alpha0, col = "gray", main = "", xlab = "alpha0", las = 1)
abline(v = data$alpha0, lwd = 3, col = "red")
hist(out$BUGSoutput$sims.list$alpha1, col = "gray", main = "", xlab = "alpha1", las = 1)
abline(v = data$alpha1, lwd = 3, col = "red")
hist(out$BUGSoutput$sims.list$beta0, col = "gray", main = "", xlab = "beta0", las=1)
abline(v = data$beta0, lwd = 3, col = "red")
hist(out$BUGSoutput$sims.list$beta1, col = "gray", main = "", xlab = "beta1", las = 1)
abline(v = data$beta1, lwd = 3, col = "red")

# Plot predicted covariate relationship with abundance
plot(data$X, data$N, main = "", xlab = "Covariate", ylab = "Abundance", las = 1, ylim = c(0, max(data$N)), frame.plot = FALSE)
lines(data$X, data$lam, type = "l", col = "red", lwd = 3)
GLM.pred <- exp(predict(glm(apply(data$y, 1, max) ~ X + I(X^2), family = poisson, data = data)))
lines(data$X, GLM.pred, type = "l", lty = 2, col = "blue", lwd = 3)
Nmix.pred <- exp(out$BUGSoutput$mean$alpha0 + out$BUGSoutput$mean$alpha1 * data$X)
points(data$X, Nmix.pred, type = "l", col = "blue", lwd = 3)
 

# 12.3. Analysis of real data: Open-population binomial-mixture model
# Get the data and put them into 3D array
bdat <- read.table("fritillary.txt", header = TRUE)
y <- array(NA, dim = c(95, 2, 7))	# 95 sites, 2 reps, 7 days

for(k in 1:7){
   sel.rows <- bdat$day == k
   y[,,k] <- as.matrix(bdat)[sel.rows, 3:4]
   }
y                          # Look at data set in 3D layout
str(y)

# Have a look at raw data
day.max <- apply(y, c(1, 3), max, na.rm = TRUE)  # Max count each site and day
day.max
site.max <- apply(day.max, 1, max, na.rm = TRUE) # Max count each site
site.max
table(site.max)         # Frequency distribution of max counts
plot(table(site.max))
table(site.max>0)       # Observed occupancy is only 56%

# Sum of observed max as conventional estimator of total abundance
max1 <- apply(y, c(1, 3), max)
obs.max.sum <- apply(max1, 2, sum, na.rm = TRUE)
obs.max.sum


# 12.3.1. Simple Poisson model
# Specify model in BUGS language
sink("Nmix0.jags")
cat("
model {

# Priors
for (k in 1:7){
   alpha.lam[k] ~ dnorm(0, 0.01)
   p[k] ~ dunif(0, 1)
   }

# Likelihood
# Ecological model for true abundance
for (k in 1:7){                          # Loop over days (7)
   lambda[k] <- exp(alpha.lam[k])
   for (i in 1:R){                       # Loop over R sites (95)
      N[i,k] ~ dpois(lambda[k])          # Abundance

      # Observation model for replicated counts
      for (j in 1:T){                    # Loop over temporal reps (2)
         y[i,j,k] ~ dbin(p[k], N[i,k])   # Detection

         # Assess model fit using Chi-squared discrepancy
         # Compute fit statistic E for observed data
         eval[i,j,k] <- p[k] * N[i,k]   	# Expected values
         E[i,j,k] <- pow((y[i,j,k] - eval[i,j,k]),2) / (eval[i,j,k] + 0.5)
         # Generate replicate data and compute fit stats for them
         y.new[i,j,k] ~ dbin(p[k], N[i,k])
         E.new[i,j,k] <- pow((y.new[i,j,k] - eval[i,j,k]),2) / (eval[i,j,k] + 0.5)

         } #j 
      } #i 
   } #k 

# Derived and other quantities
for (k in 1:7){
   totalN[k] <- sum(N[,k])	# Total pop. size across all sites
   mean.abundance[k] <- exp(alpha.lam[k])
   }
fit <- sum(E[,,])
fit.new <- sum(E.new[,,])
}
",fill = TRUE)
sink()

# Bundle data
R = nrow(y)
T = ncol(y)
win.data <- list(y = y, R = R, T = T)

# Initial values
Nst <- apply(y, c(1, 3), max) + 1
Nst[is.na(Nst)] <- 1
inits <- function(){list(N = Nst, alpha.lam = runif(7, -1, 1))}

# Parameters monitored
params <- c("totalN", "mean.abundance", "alpha.lam", "p", "fit", "fit.new")

# MCMC settings
ni <- 10000
nt <- 8
nb <- 2000
nc <- 3

# Call JAGS from R (BRT 1 min)
out0 <- jags(win.data, inits, params, "Nmix0.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

# Summarize posteriors
print(out0, dig = 3)

# Evaluation of fit
plot(out0$BUGSoutput$sims.list$fit, out0$BUGSoutput$sims.list$fit.new, main = "", xlab = "Discrepancy actual data", ylab = "Discrepancy replicate data", frame.plot = FALSE)
abline(0, 1, lwd = 2, col = "black")
mean(out0$BUGSoutput$sims.list$fit.new > out0$BUGSoutput$sims.list$fit)
mean(out0$BUGSoutput$mean$fit) / mean(out0$BUGSoutput$mean$fit.new)


# 12.3.2. Zero-inflated Poisson binomial-mixture model (ZIP binomial-mixture model)
# Specify model in BUGS language
sink("Nmix1.jags")
cat("
model {

# Priors
omega ~ dunif(0, 1)
for (k in 1:7){
   alpha.lam[k] ~ dnorm(0, 0.01)
   p[k] ~ dunif(0, 1)
   }

# Likelihood
# Ecological model for true abundance
for (i in 1:R){                          # Loop over R sites (95)
   z[i] ~ dbern(omega)                   # Latent suitability state
   for (k in 1:7){                       # Loop over survey periods (seasons)
      N[i,k] ~ dpois(lam.eff[i,k])       # Latent abundance state
      lam.eff[i,k] <- z[i] * lambda[i,k]
      log(lambda[i,k]) <- alpha.lam[k]
      # Observation model for replicated counts
      for (j in 1:T){                    # Loop over temporal reps (2)
         y[i,j,k] ~ dbin(p[k], N[i,k])   # Detection
         # Assess model fit using Chi-squared discrepancy
         # Compute fit statistic for observed data
         eval[i,j,k] <- p[k] * N[i,k]
         E[i,j,k] <- pow((y[i,j,k] - eval[i,j,k]),2) / (eval[i,j,k] + 0.5)
         # Generate replicate data and compute fit stats for them
         y.new[i,j,k] ~ dbin(p[k], N[i,k])
         E.new[i,j,k] <- pow((y.new[i,j,k] - eval[i,j,k]),2) / (eval[i,j,k]+0.5)
         } #j
      } #k
   } #i

# Derived and other quantities
for (k in 1:7){
   totalN[k] <- sum(N[,k])	# Estimate total pop. size across all sites
   mean.abundance[k] <- exp(alpha.lam[k])
   }
fit <- sum(E[,,])
fit.new <- sum(E.new[,,])
}
",fill = TRUE)
sink()

# Bundle data
R = nrow(y)
T = ncol(y)
win.data <- list(y = y, R = R, T = T)

# Initial values
# It is advisable to give initial values for the latent state z, the best option is to provide a vector of 1
Nst <- apply(y, c(1, 3), max) + 1
Nst[is.na(Nst)] <- 1
inits <- function(){list(N = Nst, alpha.lam = runif(7, -1, 1), z = rep(1, 95))}

# Parameters monitored
params <- c("omega", "totalN", "alpha.lam", "p", "mean.abundance", "fit", "fit.new")

# MCMC settings
ni <- 30000
nt <- 15
nb <- 15000
nc <- 3

# Call JAGS from R (BRT 3 min)
out1 <- jags(win.data, inits, params, "Nmix1.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

# Summarize posteriors
print(out1, dig = 3)

# Evaluation of fit
plot(out1$BUGSoutput$sims.list$fit, out1$BUGSoutput$sims.list$fit.new, main = "", xlab = "Discrepancy actual data", ylab = "Discrepancy replicate data", frame.plot = FALSE)
abline(0, 1, lwd = 2, col = "black")
mean(out1$BUGSoutput$sims.list$fit.new > out1$BUGSoutput$sims.list$fit)
mean(out1$BUGSoutput$mean$fit) / mean(out1$BUGSoutput$mean$fit.new)



# 12.3.3. Binomial-mixture model with overdispersion in both abundance and detection
# Specify model in BUGS language
sink("Nmix2.jags")
cat("
model{
# Priors
for (k in 1:7){
    alpha.lam[k] ~ dnorm(0, 0.1)
    beta[k] ~ dnorm(0, 0.1)
   }

# Abundance site and detection site-by-day random effects
for (i in 1:R){
   eps[i] ~ dnorm(0, tau.lam)                    # Abundance noise
   }
tau.lam <- 1 / (sd.lam * sd.lam)
sd.lam ~ dunif(0, 3)
tau.p <- 1 / (sd.p * sd.p)
sd.p ~ dunif(0, 3)

# Likelihood
# Ecological model for true abundance
for (i in 1:R){                                 # Loop over R sites (95)
   for (k in 1:7){                              # Loop over days (7)
      N[i,k] ~ dpois(lambda[i,k])               # Abundance
      log(lambda[i,k]) <- alpha.lam[k] + eps[i]

      # Observation model for replicated counts
      for (j in 1:T){                           # Loop over temporal reps (2)
         y[i,j,k] ~ dbin(p[i,j,k], N[i,k])      # Detection
         p[i,j,k] <- 1 / (1 + exp(-lp[i,j,k])) 
         lp[i,j,k] ~ dnorm(beta[k], tau.p) # random delta defined implicitly

         # Assess model fit using Chi-squared discrepancy
         # Compute fit statistic for observed data
         eval[i,j,k] <- p[i,j,k] * N[i,k]
         E[i,j,k] <- pow((y[i,j,k] - eval[i,j,k]),2) / (eval[i,j,k]+0.5)
         # Generate replicate data and compute fit stats for them
         y.new[i,j,k] ~ dbin(p[i,j,k], N[i,k])
         E.new[i,j,k] <- pow((y.new[i,j,k] - eval[i,j,k]),2) / (eval[i,j,k]+0.5)
         } #j
         ik.p[i,k] <- mean(p[i,,k])
      } #k
   } #i

# Derived and other quantities
for (k in 1:7){
   totalN[k] <- sum(N[,k])   # Estimate total pop. size across all sites
   mean.abundance[k] <- mean(lambda[,k])
   mean.N[k] <- mean(N[,k])
   mean.detection[k] <- mean(ik.p[,k])
   }
fit <- sum(E[,,])
fit.new <- sum(E.new[,,])
}
",fill = TRUE)
sink()

# Bundle data
R = nrow(y)
T = ncol(y)
win.data <- list(y = y, R = R, T = T)

# Initial values
Nst <- apply(y, c(1, 3), max) + 1
Nst[is.na(Nst)] <- 1
inits <- function(){list(N = Nst, alpha.lam = runif(7, -3, 3), beta = runif(7, -3, 3), sd.lam = runif(1, 0, 1), sd.p = runif(1, 0, 1))}

# Parameters monitored
params <- c("totalN", "alpha.lam", "beta", "sd.lam", "sd.p", "mean.abundance", "mean.N", "mean.detection", "fit", "fit.new")

# MCMC settings
ni <- 350000
nt <- 300
nb <- 50000
nc <- 3

# Call JAGS from R (BRT 215 min)
out2 <- jags(win.data, inits, params, "Nmix2.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

# Evaluation of fit
plot(out2$BUGSoutput$sims.list$fit, out2$BUGSoutput$sims.list$fit.new, main = "", xlab = "Discrepancy actual data", ylab = "Discrepancy replicate data", frame.plot = FALSE, xlim = c(50, 200), ylim = c(50, 200))
abline(0, 1, lwd = 2, col = "black")
mean(out2$BUGSoutput$sims.list$fit.new > out2$BUGSoutput$sims.list$fit)
mean(out2$BUGSoutput$mean$fit) / mean(out2$BUGSoutput$mean$fit.new)

# Summarize posteriors
print(out2, dig = 2)

max.day.count <- apply(y, c(1, 3), max, na.rm = TRUE)
max.day.count[max.day.count == "-Inf"] <- NA
mean.max.count <- apply(max.day.count, 2, mean, na.rm = TRUE)
mean.max.count

par(mfrow = c(2, 1))
plot(1:7, mean.max.count, xlab = "Day", ylab = "Mean daily abundance", las = 1, ylim = c(0, 16), type = "b", main = "", frame.plot = FALSE, pch = 16, lwd = 2)
lines(1:7, out2$BUGSoutput$summary[25:31,5], type = "b", pch = 16, col = "blue", lwd = 2)
segments(1:7, out2$BUGSoutput$summary[25:31,3], 1:7, out2$BUGSoutput$summary[25:31,7], col = "blue")

plot(1:7, out2$BUGSoutput$summary[32:38,1], xlab = "Day", ylab = "Detection probability ", las = 1, ylim = c(0, 1), type = "b", col = "blue", pch = 16, frame.plot = FALSE, lwd = 2)
segments(1:7, out2$BUGSoutput$summary[32:38,3], 1:7, out2$BUGSoutput$summary[32:38,7], col = "blue")





###################################################################################################################
# 
# 13. Metapopulation modeling of species distributions using hierarchical logistic regression: Site-occupancy models
# 
####################################################################################################################

# 13.2. What happens when p<1 and constant and p is not accounted for in a species distribution model?
nreps <- 10^5                             # No. replicates
estimates <- array(NA, dim = c(nreps, 2)) # Array to contain the estimates
R <- 250                                  # No. sites

for (i in 1:nreps) {
   cat(i, "\n")  ;  flush.console()
   x <- runif(R, 0, 10) # choose covariate values
   state<-rbinom(n = R, size = 1, prob = plogis(-3 + 1 * x)) # Occ. state
   obs <- rbinom(n = R, size = 1, prob = 0.6) * state # Observations
   fm <- glm(obs~x, family = binomial)
   estimates[i,] <- fm$coef
   }

par(mfrow = c(3, 1))
hist(estimates[,1], col = "gray", nclass = 50, main = "", xlab = "Intercept estimates", las = 1, ylab = "", freq = FALSE)
abline(v = -3, col = "red", lwd = 3)	# Truth
hist(estimates[,2], col = "gray", nclass = 50, main = "", xlab = "Slope estimates", xlim = c(0,1), las = 1, ylab = "", freq = FALSE)
abline(v = 1, col = "red", lwd = 3)		# Truth

plot(1:10, plogis(estimates[1,1] + estimates[1,2] * (1:10)), col = "gray", lwd = 1, ylab = "Occupancy probability", xlab = "Covariate value", type = "l", ylim = c(0, 1), frame.plot = FALSE, las = 1)
samp <- sample(1:nreps, 1000)
for (i in samp){
   lines(1:10, plogis(estimates[i,1] + estimates[i,2] * (1:10)), col = "gray", lwd = 1, type = "l")
   }
lines(1:10, plogis(-3 + 1 * (1:10)), col = "red", lwd = 3, type = "l")


# 13.3. Generation and analysis of simulated data for single-season occupancy
# 13.3.1. The simplest possible site-occupancy model
# Select sample sizes (spatial and temporal replication)
R <- 200
T <- 3

# Determine process parameters
psi <- 0.8    # Occupancy probability
p <- 0.5      # Detection probability

# Create structure to contain counts
y <- matrix(NA, nrow = R, ncol = T)

# Ecological process: Sample true occurrence (z, yes/no) from a Bernoulli (occurrence probability = psi)
z <- rbinom(n = R, size = 1, prob = psi)  # Latent occurrence state

# Observation process: Sample detection/nondetection observations from a Bernoulli(with p) if z=1
for (j in 1:T){
   y[,j] <- rbinom(n = R, size = 1, prob = z * p)
   }

# Look at truth and at our imperfect observations
sum(z)                 # Realized occupancy among 200 surveyed sites
sum(apply(y, 1, max))  # Observed occupancy

# Specify model in BUGS language
sink("model.jags")
cat("
model {

# Priors
psi ~ dunif(0, 1)
p ~ dunif(0, 1)

# Likelihood
# Ecological model for true occurrence
for (i in 1:R) {
   z[i] ~ dbern(psi)
   p.eff[i] <- z[i] * p

   # Observation model for replicated detection/nondetection observations
   for (j in 1:T) {
      y[i,j] ~ dbern(p.eff[i])
      } #j
   } #i

# Derived quantities
occ.fs <- sum(z[])       # Number of occupied sites among the 200
}
",fill = TRUE)
sink()

# Bundle data
win.data <- list(y = y, R = nrow(y), T = ncol(y))

# Initial values
zst <- apply(y, 1, max)		# Observed occurrence as starting values for z
inits <- function() list(z = zst)

# Parameters monitored
params <- c("psi", "p", "occ.fs")

# MCMC settings
ni <- 1200
nt <- 2
nb <- 200
nc <- 3

# Call JAGS from R (BRT < 1 min)
out <- jags(win.data, inits, params, "model.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

# Summarize posteriors
print(out, dig = 2)

# 13.3.2. Site-occupancy models with covariates
# Define function for generating species distribution data
data.fn <- function(R = 200, T = 3, xmin = -1, xmax = 1, alpha.psi = -1, beta.psi = 3, alpha.p = 1, beta.p = -3) {

   y <- array(dim = c(R, T))	# Array for counts

   # Ecological process
   # Covariate values
   X <- sort(runif(n = R, min = xmin, max = xmax))

   # Relationship expected occurrence - covariate
   psi <- plogis(alpha.psi + beta.psi * X)	# Apply inverse logit

   # Add Bernoulli noise: draw occurrence indicator z from Bernoulli(psi)
   z <- rbinom(n = R, size = 1, prob = psi)
   occ.fs <- sum(z)	# Finite-sample occupancy (see Royle and K?ry 2007)

   # Observation process
   # Relationship detection prob - covariate
   p <- plogis(alpha.p + beta.p * X)

   # Make a ?census?
   p.eff <- z * p
   for (i in 1:T){
      y[,i] <- rbinom(n = R, size = 1, prob = p.eff)
      }

   # Na?ve regression
   naive.pred <- plogis(predict(glm(apply(y, 1, max) ~ X + I(X^2), family = binomial)))

   # Plot features of the simulated system
   par(mfrow = c(2, 2))
   plot(X, psi, main = "Expected occurrence", xlab = "Covariate", ylab = "Occupancy probability", las = 1, type = "l", col = "red", lwd = 3, frame.plot = FALSE)
   plot(X, z, main = "Realised (true) occurrence", xlab = "Covariate", ylab = "Occurrence", las = 1, frame.plot = FALSE, col = "red",)
   plot(X, p, ylim = c(0,1), main = "Detection probability", xlab = "Covariate", ylab = "p", type = "l", lwd = 3, col = "red", las = 1, frame.plot = FALSE)
   plot(X, naive.pred, main = "Detection/nondetection observations \n and conventional SDM", xlab = "Covariate", ylab = "Apparent occupancy", ylim = c(min(y), max(y)), type = "l", lwd = 3, lty = 2, col = "blue", las = 1, frame.plot = FALSE)
   points(rep(X, T), y)

   # Return stuff
   return(list(R = R, T = T, X = X, alpha.psi = alpha.psi, beta.psi = beta.psi, alpha.p = alpha.p , beta.p = beta.p, psi = psi, z = z, occ.fs = occ.fs, p = p, y = y))
   }

sodata <- data.fn()
str(sodata)                 # Look at data

summary(glm(apply(y, 1, max) ~ X + I(X^2), family = binomial, data = sodata))

# Specify model in BUGS language
sink("model.jags")
cat("
model {

# Priors
alpha.occ ~ dunif(-10, 10)
beta.occ ~ dunif(-10, 10)
alpha.p ~ dunif(-10, 10)
beta.p ~ dunif(-10, 10)

# Likelihood
for (i in 1:R) {
   # True state model for the partially observed true state
   z[i] ~ dbern(psi[i])             # True occupancy z at site i
   logit(psi[i]) <- alpha.occ + beta.occ * X[i]

   for (j in 1:T) {
      # Observation model for the actual observations
      y[i,j] ~ dbern(p.eff[i,j])    # Detection-nondetection at i and j
      p.eff[i,j] <- z[i] * p[i,j]
      logit(p[i,j]) <- alpha.p + beta.p * X[i]
      } #j
   } #i

# Derived quantities
occ.fs <- sum(z[])       # Number of occupied sites among those studied
}
",fill = TRUE)
sink()

# Bundle data
win.data <- list(y = sodata$y, X = sodata$X, R = nrow(sodata$y), T = ncol(sodata$y))

# Initial values
zst <- apply(sodata$y, 1, max)   #Good inits for latent states essential
inits <- function(){list(z = zst, alpha.occ = runif(1, -3, 3), beta.occ = runif(1, -3, 3), alpha.p = runif(1, -3, 3), beta.p = runif(1, -3, 3))}

# Parameters monitored
params <- c("alpha.occ", "beta.occ", "alpha.p", "beta.p", "occ.fs")

# MCMC settings
ni <- 10000
nt <- 8
nb <- 2000
nc <- 3

# Call JAGS from R (BRT 1 min)
out <- jags(win.data, inits, params, "model.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

TRUTH <- c(sodata$alpha.psi, sodata$beta.psi, sodata$alpha.p, sodata$beta.p, sum(sodata$z))
print(cbind(TRUTH, out$BUGSoutput$summary[c(1:4, 6), c(1,2,3,7)]), dig = 3)
sum(apply(sodata$y, 1, sum) > 0)# Apparent number of occupied sites

naive.pred <- plogis(predict(glm(apply(sodata$y, 1, max) ~ X + I(X^2), family = binomial, data = sodata)))
lin.pred2 <- out$BUGSoutput$mean$alpha.occ + out$BUGSoutput$mean$beta.occ * sodata$X

plot(sodata$X, sodata$psi, ylim = c(0, 1), main = "", ylab = "Occupancy probability", xlab = "Covariate", type = "l", lwd = 3, col = "red", las = 1, frame.plot = FALSE)
lines(sodata$X, naive.pred, ylim = c(0 ,1), type = "l", lty = 2, lwd = 3, col = "blue")
lines(sodata$X, plogis(lin.pred2), ylim = c(0, 1), type = "l", lty = 1, lwd = 2, col = "blue")
 

# 13.4. Analysis of real data set: Single-season occupancy model
# Read in the data
data <- read.table("bluebug.txt", header = TRUE)

# Collect the data into suitable structures
y <- as.matrix(data[,4:9])         # as.matrix essential for WinBUGS
y[y>1] <- 1                        # Reduce counts to 0/1
edge <- data$forest_edge
dates <- as.matrix(data[,10:15])
hours <- as.matrix(data[,16:21])

# Standardize covariates
mean.date <- mean(dates, na.rm = TRUE)
sd.date <- sd(dates[!is.na(dates)])
DATES <- (dates-mean.date)/sd.date     # Standardise date
DATES[is.na(DATES)] <- 0               # Impute zeroes (means)

mean.hour <- mean(hours, na.rm = TRUE)
sd.hour <- sd(hours[!is.na(hours)])
HOURS <- (hours-mean.hour)/sd.hour      # Standardise hour
HOURS[is.na(HOURS)] <- 0                # Impute zeroes (means)

# Specify model in BUGS language
sink("model.jags")
cat("
model {

# Priors
alpha.psi ~ dnorm(0, 0.01)
beta.psi ~ dnorm(0, 0.01)
alpha.p ~ dnorm(0, 0.01)
beta1.p ~ dnorm(0, 0.01)
beta2.p ~ dnorm(0, 0.01)
beta3.p ~ dnorm(0, 0.01)
beta4.p ~ dnorm(0, 0.01)

# Likelihood
# Ecological model for the partially observed true state
for (i in 1:R) {
   z[i] ~ dbern(psi[i])                # True occurrence z at site i
   psi[i] <- 1 / (1 + exp(-lpsi.lim[i]))
   lpsi.lim[i] <- min(999, max(-999, lpsi[i]))
   lpsi[i] <- alpha.psi + beta.psi * edge[i]

   # Observation model for the observations
   for (j in 1:T) {
      y[i,j] ~ dbern(mu.p[i,j])	# Detection-nondetection at i and j
      mu.p[i,j] <- z[i] * p[i,j]
      p[i,j] <- 1 / (1 + exp(-lp.lim[i,j]))
      lp.lim[i,j] <- min(999, max(-999, lp[i,j]))
      lp[i,j] <- alpha.p + beta1.p * DATES[i,j] + beta2.p * pow(DATES[i,j], 2) + beta3.p * HOURS[i,j] + beta4.p * pow(HOURS[i,j], 2)
      } #j
   } #i

# Derived quantities
occ.fs <- sum(z[])                             # Number of occupied sites
mean.p <- exp(alpha.p) / (1 + exp(alpha.p))    # Sort of average detection
}
",fill = TRUE)
sink()

# Bundle data
win.data <- list(y = y, R = nrow(y), T = ncol(y), edge = edge, DATES = DATES, HOURS = HOURS)

# Initial values
zst <- apply(y, 1, max, na.rm = TRUE)	# Good starting values crucial
inits <- function(){list(z = zst, alpha.psi=runif(1, -3, 3), alpha.p = runif(1, -3, 3))}

# Parameters monitored
params <- c("alpha.psi", "beta.psi", "mean.p", "occ.fs", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p")

# MCMC settings
ni <- 30000
nt <- 10
nb <- 20000
nc <- 3

# Call JAGS from R (BRT < 1 min)
out <- jags(win.data, inits, params, "model.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

# Summarize posteriors
print(out, dig = 2)

# Posterior distribution of the number of occupied sites in actual sample
hist(out$BUGSoutput$sims.list$occ.fs, nclass = 30, col = "gray", main = "", xlab = "Number of occupied woodpiles (occ.fs)", xlim = c(9, 27))
abline(v = 10, lwd = 2) # The observed number
 
Pstar <- array(NA, dim = c(out$BUGSoutput$n.sims, 10))
x <- cbind(rep(1, 3000), rep(2, 3000), rep(3, 3000), rep(4, 3000), rep(5, 3000), rep(6, 3000), rep(7, 3000), rep(8, 3000), rep(9, 3000), rep(10, 3000)) 
for (i in 1:out$BUGSoutput$n.sims) {
   for (j in 1:10){
      Pstar[i,j] <- 1 - (1 - out$BUGSoutput$sims.list$mean.p[i])^j
      } #j
   } #i

boxplot(Pstar ~ x, col = "gray", las = 1, ylab = "Pstar", xlab = "Number of surveys", outline = FALSE)
abline(h = 0.95, lty = 2, lwd = 2)
 
par(mfrow = c(2, 1))
hist(plogis(out$BUGSoutput$sims.list$alpha.psi), nclass = 40, col = "gray", main = "Forest interior", xlab = "Occupancy probability", xlim = c(0, 1))
hist(plogis(out$BUGSoutput$sims.list$alpha.psi+ out$BUGSoutput$sims.list$beta.psi), nclass = 40, col = "gray", main = "Forest edge", xlab = "Occupancy probability", xlim = c(0, 1))
 
# Predict effect of time of day with uncertainty
mcmc.sample <- out$BUGSoutput$n.sims

original.date.pred <- seq(0, 60, length.out = 30)
original.hour.pred <- seq(180, 540, length.out = 30)
date.pred <- (original.date.pred - mean.date)/sd.date
hour.pred <- (original.hour.pred - mean.hour)/sd.hour
p.pred.date <- plogis(out$BUGSoutput$mean$alpha.p + out$BUGSoutput$mean$beta1.p * date.pred + out$BUGSoutput$mean$beta2.p * date.pred^2 )
p.pred.hour <- plogis(out$BUGSoutput$mean$alpha.p + out$BUGSoutput$mean$beta3.p * hour.pred + out$BUGSoutput$mean$beta4.p * hour.pred^2 )

array.p.pred.hour <- array.p.pred.date <- array(NA, dim = c(length(hour.pred), mcmc.sample))
for (i in 1:mcmc.sample){
   array.p.pred.date[,i] <- plogis(out$BUGSoutput$sims.list$alpha.p[i] + out$BUGSoutput$sims.list$beta1.p[i] * date.pred + out$BUGSoutput$sims.list$beta2.p[i] * date.pred^2)
   array.p.pred.hour[,i] <- plogis(out$BUGSoutput$sims.list$alpha.p[i] + out$BUGSoutput$sims.list$beta3.p[i] * hour.pred + out$BUGSoutput$sims.list$beta4.p[i] * hour.pred^2)
   }

# Plot for a subsample of MCMC draws
sub.set <- sort(sample(1:mcmc.sample, size = 200))

par(mfrow = c(2, 1))
plot(original.date.pred, p.pred.date, main = "", ylab = "Detection probability", xlab = "Date (1 = 1 July)", ylim = c(0, 1), type = "l", lwd = 3, frame.plot = FALSE)
for (i in sub.set){
   lines(original.date.pred, array.p.pred.date[,i], type = "l", lwd = 1, col = "gray")
   }
lines(original.date.pred, p.pred.date, type = "l", lwd = 3, col = "blue")

plot(original.hour.pred, p.pred.hour, main = "", ylab = "Detection probability", xlab = "Time of day (mins after noon)", ylim = c(0, 1), type = "l", lwd = 3, frame.plot = FALSE)
for (i in sub.set){
   lines(original.hour.pred, array.p.pred.hour[,i], type = "l", lwd = 1, col = "gray")
   }
lines(original.hour.pred, p.pred.hour, type = "l", lwd = 3, col = "blue")

 
# 13.5. Dynamic (multi-season) site-occupancy models
# 13.5.1. Generation and analysis of simulated data
data.fn <- function(R = 250, J = 3, K = 10, psi1 = 0.4, range.p = c(0.2, 0.4), range.phi = c(0.6, 0.8), range.gamma = c(0, 0.1)) {
# Function to simulate detection/nondetection data for dynamic site-occ model
# Annual variation in probabilities of patch survival, colonization and 
# detection is specified by the bounds of a uniform distribution.

# Function arguments:
# R - Number of sites
# J - Number of replicate surveys
# K - Number of years
# psi1 - occupancy probability in first year
# range.p - bounds of uniform distribution from which annual p drawn 
# range.psi and range.gamma - same for survival and colonization probability

   # Set up some required arrays
   site <- 1:R					# Sites
   year <- 1:K					# Years
   psi <- rep(NA, K)				# Occupancy probability
   muZ <- z <- array(dim = c(R, K))	# Expected and realized occurrence
   y <- array(NA, dim = c(R, J, K))	# Detection histories

   # Determine initial occupancy and demographic parameters
   psi[1] <- psi1				# Initial occupancy probability
   p <- runif(n = K, min = range.p[1], max = range.p[2])
   phi <- runif(n = K-1, min = range.phi[1], max = range.phi[2])
   gamma <- runif(n = K-1, min = range.gamma[1], max = range.gamma[2])

   # Generate latent states of occurrence
   # First year
   z[,1] <- rbinom(R, 1, psi[1])		# Initial occupancy state
   # Later years
   for(i in 1:R){				# Loop over sites
      for(k in 2:K){				# Loop over years
         muZ[k] <- z[i, k-1]*phi[k-1] + (1-z[i, k-1])*gamma[k-1] # Prob for occ.
         z[i,k] <- rbinom(1, 1, muZ[k])
         }
      }

   # Plot realised occupancy
   plot(year, apply(z, 2, mean), type = "l", xlab = "Year", ylab = "Occupancy or Detection prob.", col = "red", xlim = c(0,K+1), ylim = c(0,1), lwd = 2, lty = 1, frame.plot = FALSE, las = 1)
   lines(year, p , type = "l", col = "red", lwd = 2, lty = 2)

   # Generate detection/nondetection data
   for(i in 1:R){
      for(k in 1:K){
         prob <- z[i,k] * p[k]
         for(j in 1:J){
            y[i,j,k] <- rbinom(1, 1, prob)
            }
         }
      }

   # Compute annual population occupancy
   for (k in 2:K){
      psi[k] <- psi[k-1]*phi[k-1] + (1-psi[k-1])*gamma[k-1]
      }

   # Plot apparent occupancy
   psi.app <- apply(apply(y, c(1,3), max), 2, mean)
   lines(year, psi.app, type = "l", col = "black", lwd = 2)
   text(0.85*K, 0.06, labels = "red solid - true occupancy\n red dashed - detection\n black - observed occupancy")

   # Return data
   return(list(R = R, J = J, K = K, psi = psi, psi.app = psi.app, z = z, phi = phi, gamma = gamma, p = p, y = y))
}

data <- data.fn(R = 250, J = 3, K = 10, psi1 = 0.6, range.p = c(0.1, 0.9), range.phi = c(0.7, 0.9), range.gamma = c(0.1, 0.5))

attach(data)
str(data)

# Specify model in BUGS language
sink("Dynocc.jags")
cat("
model {

# Specify priors
psi1 ~ dunif(0, 1)
for (k in 1:(nyear-1)){
   phi[k] ~ dunif(0, 1)
   gamma[k] ~ dunif(0, 1)
   p[k] ~ dunif(0, 1) 
   }
p[nyear] ~ dunif(0, 1)

# Ecological submodel: Define state conditional on parameters
for (i in 1:nsite){
   z[i,1] ~ dbern(psi1)
   for (k in 2:nyear){
      muZ[i,k]<- z[i,k-1]*phi[k-1] + (1-z[i,k-1])*gamma[k-1]
      z[i,k] ~ dbern(muZ[i,k])
      } #k
   } #i

# Observation model
for (i in 1:nsite){
   for (j in 1:nrep){
      for (k in 1:nyear){
         muy[i,j,k] <- z[i,k]*p[k]
         y[i,j,k] ~ dbern(muy[i,j,k])
         } #k
      } #j
   } #i

# Derived parameters: Sample and population occupancy, growth rate and turnover
psi[1] <- psi1
n.occ[1]<-sum(z[1:nsite,1])
for (k in 2:nyear){
   psi[k] <- psi[k-1]*phi[k-1] + (1-psi[k-1])*gamma[k-1]
   n.occ[k] <- sum(z[1:nsite,k])
   growthr[k-1] <- psi[k]/psi[k-1]                         # originally we had growthr[k]. JAGS seem to dislike vectoring going from 2..K.
   turnover[k-1] <- (1 - psi[k-1]) * gamma[k-1]/psi[k]
   }
}
",fill = TRUE)
sink()

# Bundle data
win.data <- list(y = y, nsite = dim(y)[1], nrep = dim(y)[2], nyear = dim(y)[3])

# Initial values
zst <- apply(y, c(1, 3), max)	# Observed occurrence as inits for z
inits <- function(){ list(z = zst)}

# Parameters monitored
params <- c("psi", "phi", "gamma", "p", "n.occ", "growthr", "turnover") 


# MCMC settings
ni <- 2500
nt <- 4
nb <- 500
nc <- 3

# Call JAGS from R (BRT 3 min)
out <- jags(win.data, inits, params, "Dynocc.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

# Summarize posteriors
print(out, dig = 2)
psiall <- paste("psi[", 1:K, "]", sep="")
print(cbind(data$psi, out$BUGSoutput$summary[psiall, c(1, 2, 3, 7)]), dig = 3)
phiall <- paste("phi[", 1:(K-1), "]", sep="")
print(cbind(data$phi, out$BUGSoutput$summary[phiall, c(1, 2, 3, 7)]), dig = 3)
gammaall <- paste("gamma[", 1:(K-1), "]", sep="")
print(cbind(data$gamma, out$BUGSoutput$summary[gammaall, c(1, 2, 3, 7)]), dig = 3)
pall <- paste("p[", 1:K, "]", sep="")
print(cbind(data$p, out$BUGSoutput$summary[pall, c(1, 2, 3, 7)]), dig = 3)

plot(1:K, data$psi, type = "l", xlab = "Year", ylab = "Occupancy probability", col = "red", xlim = c(0,K+1), ylim = c(0,1), lwd = 2, lty = 1, frame.plot = FALSE, las = 1)
lines(1:K, data$psi.app, type = "l", col = "black", lwd = 2)
points(1:K, out$BUGSoutput$mean$psi, type = "l", col = "blue", lwd = 2)
segments(1:K, out$BUGSoutput$summary[psiall,3], 1:K, out$BUGSoutput$summary[psiall,7], col = "blue", lwd = 1)

# 13.5.2. Dynamic occupancy modeling in a real data set
# Read in the data and put it into 3D array
bdat <- read.table(file = "burnet.txt", header = T)
str(bdat)

y <- array(NA, dim = c(95, 2, 7))	# 95 sites, 2 reps, 7 days

for (i in 1:7){
   sel.rows <- bdat$day == i
   y[,,i] <- as.matrix(bdat)[sel.rows, 3:4]
   }
str(y)

# Convert counts to detection/nondetection data
y[y>0] <- 1

# Look at the number of sites with detections for each day
tmp <- apply(y, c(1,3), max, na.rm = TRUE)
tmp[tmp == "-Inf"] <- NA
apply(tmp, 2, sum, na.rm = TRUE)

# Bundle data
win.data <- list(y = y, nsite = dim(y)[1], nrep = dim(y)[2], nyear = dim(y)[3])

# Initial values
inits <- function(){ list(z = apply(y, c(1, 3), max))}

# Parameters monitored
params <- c("psi", "phi", "gamma", "p", "n.occ", "growthr", "turnover")  

# MCMC settings
ni <- 5000
nt <- 4
nb <- 1000
nc <- 3

# Call JAGS from R (BRT 1 min)
out1 <- jags(win.data, inits, params, "Dynocc.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

# Summarize posteriors
print(out1, dig = 3)

# Specify model in BUGS language
sink("Dynocc2.jags")
cat("
model {

# Specify priors
psi1 ~ dunif(0, 1)
for (k in 1:(nyear-1)){
   phi[k] ~ dunif(0, 1)
   gamma[k] ~ dunif(0, 1)
   }
p ~ dunif(0, 1)

# Both models at once
for (i in 1:nsite){
   z[i,1] ~ dbern(psi1)     # State model 1: Initial state
   muy[i,1] <- z[i,1]*p
   y[i,1] ~ dbin(muy[i,1], 2)
   for (k in 2:nyear){      # State model 2: State dynamics
      muZ[i,k] <- z[i,k-1]*phi[k-1] + (1-z[i,k-1])*gamma[k-1]
      z[i,k] ~ dbern(muZ[i,k])

      # Observation model
      muy[i,k] <- z[i,k]*p
      y[i,k] ~ dbin(muy[i,k], 2)
      } #k
   } #i

# Derived parameters: Sample and population occupancy, growth rate and turnover
psi[1] <- psi1
n.occ[1] <- sum(z[1:nsite,1])
for (k in 2:nyear){
   psi[k] <- psi[k-1]*phi[k-1] + (1-psi[k-1])*gamma[k-1]
   n.occ[k] <- sum(z[1:nsite,k])
   growthr[k-1] <- psi[k]/psi[k-1]
   turnover[k-1] <- (1 - psi[k-1]) * gamma[k-1]/psi[k]
   }
} # end of model
",fill = TRUE)
sink()

# Aggregate detections over reps within a day and bundle data
yy <- apply(y, c(1, 3), sum, na.rm = TRUE)
win.data <- list(y = yy, nsite = dim(yy)[1], nyear = dim(yy)[2])

# Initial values
inits <- function(){list(z = apply(y, c(1, 3), max))}

# Parameters monitored
params <- c("psi", "phi", "gamma", "p", "n.occ", "growthr", "turnover")  

# MCMC settings
ni <- 2500
nt <- 2
nb <- 500
nc <- 3

# Call JAGS from R (BRT 1 min)
out2 <- jags(win.data, inits, params, "Dynocc2.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

# Summarize posterior
print(out2, dig = 3)

DAY <- cbind(rep(1, out2$BUGSoutput$n.sims), rep(2, out2$BUGSoutput$n.sims), rep(3, out2$BUGSoutput$n.sims), rep(4, out2$BUGSoutput$n.sims), rep(5, out2$BUGSoutput$n.sims), rep(6, out2$BUGSoutput$n.sims), rep(7, out2$BUGSoutput$n.sims))
boxplot(out2$BUGSoutput$sims.list$psi ~ DAY, col = "gray", ylab = "Occupancy probability", xlab = "Day of survey", las = 1, frame.plot = FALSE)
apply(apply(y, c(1, 3), max), 2, function(x){sum(!is.na(x))})


# 13.6. Multistate occupancy models
owls <- read.table("owls.txt", header = TRUE)
str(owls)

# Specify model in BUGS language
sink("model1.jags")
cat("
model  { 

# Priors
p2 ~ dunif(0, 1)
psi ~ dunif(0, 1)
r ~ dunif(0, 1)
for (i in 1:3) {
   beta[i] ~ dgamma(1, 1)   # Induce Dirichlet prior
   p3[i] <- beta[i]/sum(beta[])
   }


# Define state vector
for (s in 1:R){
   phi[s,1] <- 1 - psi            # Prob. of non-occupation
   phi[s,2] <- psi * (1-r)        # Prob. of occupancy without repro
   phi[s,3] <- psi * r            # Prob. of occupancy and repro
   }

# Define observation matrix
# Order of indices: true state, time, observed state
for (t in 1:T){
   p[1,t,1] <- 1
   p[1,t,2] <- 0
   p[1,t,3] <- 0
   p[2,t,1] <- 1-p2
   p[2,t,2] <- p2
   p[2,t,3] <- 0
   p[3,t,1] <- p3[1]
   p[3,t,2] <- p3[2]
   p[3,t,3] <- p3[3]
   }

# State-space likelihood
# State equation: model of true states (z)
for (s in 1:R){
   z[s] ~ dcat(phi[s,])
   }

# Observation equation
for (s in 1:R){
   for (t in 1:T){ 
      y[s,t] ~ dcat(p[z[s],t,])
      } #t
   } #s

# Derived quantities
for (s in 1:R){
   occ1[s] <- equals(z[s], 1)
   occ2[s] <- equals(z[s], 2)
   occ3[s] <- equals(z[s], 3)
   }
n.occ[1] <- sum(occ1[]) # Sites in state 1
n.occ[2] <- sum(occ2[]) # Sites in state 2
n.occ[3] <- sum(occ3[]) # Sites in state 3
}
",fill=TRUE)
sink()

# Bundle data
y <- as.matrix(owls[, 2:6])
y <- y + 1
win.data <- list(y = y, R = dim(y)[1], T = dim(y)[2])

# Initial values
zst <- apply(y, 1, max, na.rm = TRUE)
zst[zst == "-Inf"] <- 1
inits <- function(){list(z = zst)}

# Parameters monitored
params <- c("p2", "p3", "r", "psi", "n.occ") # Might want to add "z"

# MCMC settings
ni <- 2500
nt <- 2
nb <- 500
nc <- 3

# Call JAGS from R (BRT <1 min)
out1 <- jags(win.data, inits, params, "model1.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

# Summarize posteriors
print(out1, dig = 2)

# Specifiy model in BUGS language
sink("model2.jags")
cat("
model  { 

# Priors
psi ~ dunif(0, 1)
r ~ dunif(0,1 )

for (t in 1:T){
   p2[t] ~ dunif(0, 1)
   for (i in 1:3) {
      beta[i,t] ~ dgamma(1, 1)   # Induce Dirichlet prior
      p3[i,t] <- beta[i,t]/sum(beta[,t])
      } #i
   } #t

# Define state vector
for (s in 1:R){
   phi[s,1] <- 1 - psi              # Prob. of non-occupation
   phi[s,2] <- psi * (1-r)          # Prob. of occupancy without repro.
   phi[s,3] <- psi * r              # Prob. of occupancy and repro
   }

# Define observation matrix
# Order of indices: true state, time, observed state
for (t in 1:T){    
   p[1,t,1] <- 1
   p[1,t,2] <- 0
   p[1,t,3] <- 0
   p[2,t,1] <- 1-p2[t]
   p[2,t,2] <- p2[t]
   p[2,t,3] <- 0
   p[3,t,1] <- p3[1,t]
   p[3,t,2] <- p3[2,t]
   p[3,t,3] <- p3[3,t]
   }

# State-space likelihood
# State equation: model of true states (z)
for (s in 1:R){
   z[s] ~ dcat(phi[s,])
   }

# Observation equation
for (s in 1:R){
   for (t in 1:T){ 
      y[s,t] ~ dcat(p[z[s],t,])
      } #t
   } #s

# Derived quantities
for (s in 1:R){
   occ1[s] <- equals(z[s], 1)
   occ2[s] <- equals(z[s], 2)
   occ3[s] <- equals(z[s], 3)
   }
n.occ[1] <- sum(occ1[]) # Sites in state 1
n.occ[2] <- sum(occ2[]) # Sites in state 2
n.occ[3] <- sum(occ3[]) # Sites in state 3
}
",fill=TRUE)
sink()

# Call JAGS from R (BRT 1 min)
out2 <- jags(win.data, inits, params, "model2.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

# Summarize posteriors
print(out2, dig = 2)




###################################################################################################################
# 
# Appendix 2: Two further useful multistate capture-recapture models
# 
####################################################################################################################

# 1. Estimation of age-specific survival probabilities
# 1.2. Generation of simulated data
# Define mean survival, transitions, recapture, as well as number of occasions, states, observations and released individuals
phi.juv <- 0.3
phi.ad <- 0.65
p <- 0.5
n.occasions <- 6  
n.states <- 3
n.obs <- 3
marked <- matrix(NA, ncol = n.states, nrow = n.occasions)
marked[,1] <- rep(200, n.occasions)	# Juveniles
marked[,2] <- rep(30, n.occasions)	# Adults
marked[,3] <- rep(0, n.occasions)	# Dead individuals

# Define matrices with survival, transition and recapture probabilities
# These are 4-dimensional matrices, with 
   # Dimension 1: state of departure
   # Dimension 2: state of arrival
   # Dimension 3: individual
   # Dimension 4: time
# 1. State process matrix
totrel <- sum(marked)*(n.occasions-1)
PSI.STATE <- array(NA, dim=c(n.states, n.states, totrel, n.occasions-1))
for (i in 1:totrel){
   for (t in 1:(n.occasions-1)){
      PSI.STATE[,,i,t] <- matrix(c(
      0, phi.juv, 1-phi.juv,
      0, phi.ad,  1-phi.ad,
      0, 0,       1          ), nrow = n.states, byrow = TRUE)
      } #t
   } #i
# 2.Observation process matrix
PSI.OBS <- array(NA, dim=c(n.states, n.obs, totrel, n.occasions-1))
for (i in 1:totrel){
   for (t in 1:(n.occasions-1)){
      PSI.OBS[,,i,t] <- matrix(c(
      0, 0, 1,
      0, p, 1-p,
      0, 0, 1    ), nrow = n.states, byrow = TRUE)
      } #t
   } #i

# Execute simulation function
sim <- simul.ms(PSI.STATE, PSI.OBS, marked)
CH <- sim$CH

# Compute vector with occasion of first capture
get.first <- function(x) min(which(x!=0))
f <- apply(CH, 1, get.first)

# Recode CH matrix: note, a 0 is not allowed in WinBUGS!
# 1 = seen alive as juvenile, 2 = seen alive as adult, 3 = not seen
rCH <- CH  # Recoded CH
rCH[rCH==0] <- 3

# 1.3. Analysis of the model
# Specify model in BUGS language
sink("age.jags")
cat("
model {
#----------------------------------------
# Parameters:
# phi.juv: juvenile survival probability
# phi.ad: adult survival probability
# p: recapture probability
#----------------------------------------
# States (S):
# 1 alive as juvenile
# 2 alive as adult
# 3 dead
# Observations (O): 
# 1 seen as juvenile 
# 2 seen as adult
# 3 not seen
#----------------------------------------

# Priors and constraints
for (t in 1:(n.occasions-1)){
   phi.juv[t] <- mean.phijuv
   phi.ad[t] <- mean.phiad
   p[t] <- mean.p
   }
mean.phijuv ~ dunif(0, 1)
mean.phiad ~ dunif(0, 1)
mean.p ~ dunif(0, 1)

# Define state-transition and observation matrices 	
for (i in 1:nind){
   # Define probabilities of state S(t+1) given S(t)
      for (t in f[i]:(n.occasions-1)){
      ps[1,i,t,1] <- 0
      ps[1,i,t,2] <- phi.juv[t]
      ps[1,i,t,3] <- 1-phi.juv[t]
      ps[2,i,t,1] <- 0
      ps[2,i,t,2] <- phi.ad[t]
      ps[2,i,t,3] <- 1-phi.ad[t]
      ps[3,i,t,1] <- 0
      ps[3,i,t,2] <- 0
      ps[3,i,t,3] <- 1

      # Define probabilities of O(t) given S(t)
      po[1,i,t,1] <- 0
      po[1,i,t,2] <- 0
      po[1,i,t,3] <- 1
      po[2,i,t,1] <- 0
      po[2,i,t,2] <- p[t]
      po[2,i,t,3] <- 1-p[t]
      po[3,i,t,1] <- 0
      po[3,i,t,2] <- 0
      po[3,i,t,3] <- 1
      } #t
   } #i

# State-space model likelihood
for (i in 1:nind){
   z[i,f[i]] <- y[i,f[i]]
   for (t in (f[i]+1):n.occasions){
      # State equation: draw S(t) given S(t-1)
      z[i,t] ~ dcat(ps[z[i,t-1], i, t-1,])
      # Observation equation: draw O(t) given S(t)
      y[i,t] ~ dcat(po[z[i,t], i, t-1,])
      } #t
   } #i
}
",fill = TRUE)
sink()


# Bundle data
jags.data <- list(y = rCH, f = f, n.occasions = dim(rCH)[2], nind = dim(rCH)[1])

# Initial values
# Good initial values for the latent state matrix need to be provided. The following function is a possibility that works well.
age.ms.init <- function(ch){
   init <- ch
   al <- which(ch==2, arr.ind = T)
   for (i in 1:nrow(al)){
      init[al[i,1], f[al[i,1]]:al[i,2]] <- 2
      }
   for (i in 1:nrow(init)){
      init[i,1:f[i]] <- NA
      }
   return(init)
   }
   
inits <- function(){list(mean.phijuv = runif(1, 0, 1), mean.phiad = runif(1, 0, 1), mean.p = runif(1, 0, 1), z = age.ms.init(rCH))}  

# Parameters monitored
parameters <- c("mean.phijuv", "mean.phiad", "mean.p")

# MCMC settings
ni <- 2000
nt <- 3
nb <- 1000
nc <- 3

# Call JAGS from R (BRT 2 min)
age <- jags(jags.data, inits, parameters, "age.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

print(age, digits = 3)


# 2. Accounting for immediate trap response
# 2.2. Generation of simulated data 
# Define mean survival, transitions, recapture, as well as number of occasions, states, observations and released individuals
phi <- 0.55
pss <- 0.75
pns <- 0.3
n.occasions <- 10  
n.states <- 3
n.obs <- 2
marked <- matrix(NA, ncol = n.states, nrow = n.occasions)
marked[,1] <- rep(100, n.occasions)	# Alive, seen
marked[,2] <- rep(0, n.occasions)	# Alive, not seen
marked[,3] <- rep(0, n.occasions)	# Dead 

# Define matrices with survival, transition and recapture probabilities
# These are 4-dimensional matrices, with 
   # Dimension 1: state of departure
   # Dimension 2: state of arrival
   # Dimension 3: individual
   # Dimension 4: time
# 1. State process matrix
totrel <- sum(marked)*(n.occasions-1)
PSI.STATE <- array(NA, dim=c(n.states, n.states, totrel, n.occasions-1))
for (i in 1:totrel){
   for (t in 1:(n.occasions-1)){
      PSI.STATE[,,i,t] <- matrix(c(
      phi*pss, phi*(1-pss), 1-phi,
      phi*pns, phi*(1-pns), 1-phi,
      0,       0,           1       ), nrow = n.states, byrow = TRUE)
      } #t
   } #i

# 2.Observation process matrix
PSI.OBS <- array(NA, dim=c(n.states, n.obs, totrel, n.occasions-1))
for (i in 1:totrel){
   for (t in 1:(n.occasions-1)){
      PSI.OBS[,,i,t] <- matrix(c(
      1, 0,
      0, 1,
      0, 1  ), nrow = n.states, byrow = TRUE)
      } #t
   } #i

# Execute simulation function
sim <- simul.ms(PSI.STATE, PSI.OBS, marked, unobservable = 2)
CH <- sim$CH

# Compute vector with occasion of first capture
get.first <- function(x) min(which(x!=0))
f <- apply(CH, 1, get.first)

# Recode CH matrix: note, a 0 is not allowed in WinBUGS!
# 1 = seen, 2 = not seen 
rCH <- CH  # Recoded CH
rCH[rCH==0] <- 2

# 2.3. Analysis of the model
# Specify model in BUGS language
sink("immtrap.jags")
cat("
model {
#-----------------------------------------------------------
# Parameters:
# phi: survival probability
# pss: recapture probability at t, given captured at t-1
# pns: recapture probability at t, given not captured at t-1
#-----------------------------------------------------------
# States (S):
# 1 alive, seen at t-1
# 2 alive, not seen at t-1
# 3 dead
# Observations (O): 
# 1 seen 
# 2 not seen
#-----------------------------------------------------------

# Priors and constraints
for (t in 1:(n.occasions-1)){
   phi[t] <- mean.phi
   pss[t] <- mean.pss
   pns[t] <- mean.pns
   }
mean.phi ~ dunif(0, 1)
mean.pss ~ dunif(0, 1)
mean.pns ~ dunif(0, 1)

# Define state-transition and observation matrices 	
for (i in 1:nind){
   # Define probabilities of state S(t+1) given S(t)
   for (t in f[i]:(n.occasions-1)){
      ps[1,i,t,1] <- phi[t]*pss[t]
      ps[1,i,t,2] <- phi[t]*(1-pss[t])
      ps[1,i,t,3] <- 1-phi[t]
      ps[2,i,t,1] <- phi[t]*pns[t]
      ps[2,i,t,2] <- phi[t]*(1-pns[t])
      ps[2,i,t,3] <- 1-phi[t]
      ps[3,i,t,1] <- 0
      ps[3,i,t,2] <- 0
      ps[3,i,t,3] <- 1

      # Define probabilities of O(t) given S(t)
      po[1,i,t,1] <- 1
      po[1,i,t,2] <- 0
      po[2,i,t,1] <- 0
      po[2,i,t,2] <- 1
      po[3,i,t,1] <- 0
      po[3,i,t,2] <- 1
      } #t
   } #i
# State-space model likelihood 
for (i in 1:nind){
   z[i,f[i]] <- y[i,f[i]]
   for (t in (f[i]+1):n.occasions){
      # State equation: draw S(t) given S(t-1)
      z[i,t] ~ dcat(ps[z[i,t-1], i, t-1,])
      # Observation equation: draw O(t) given S(t)
      y[i,t] ~ dcat(po[z[i,t], i, t-1,])
      } #t
   } #i
}
",fill = TRUE)
sink()

# Bundle data
jags.data <- list(y = rCH, f = f, n.occasions = dim(rCH)[2], nind = dim(rCH)[1], z = known.state.ms(rCH, 2))

# Initial values
inits <- function(){list(mean.phi = runif(1, 0, 1), mean.pss = runif(1, 0, 1), mean.pns = runif(1, 0, 1), z = 2*ms.init.z(rCH, f))}  

# Parameters monitored
parameters <- c("mean.phi", "mean.pss", "mean.pns")

# MCMC settings
ni <- 20000
nt <- 6
nb <- 10000
nc <- 3

# Call JAGS from R (BRT 33 min)
immtrap <- jags(jags.data, inits, parameters, "immtrap.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

print(immtrap, digits = 3)
