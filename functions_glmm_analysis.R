#### Functions for transforming data and analyzing model output

# Function to standardize covariates for occupancy models
standardize <- function(covar){
  mean.covar <- mean(covar, na.rm = TRUE)
  sd.covar <- sd(covar[!is.na(covar)])
  return((covar - mean.covar)/sd.covar)
}

# Make data matrices for JAGS GLMM with spatial random effect
makeEcoDataMatrix <- function(var, data = detectData, fill = NA){
  dataf <- data[,c("ymo","cellID",var)]
  dataMatrix <- eval(substitute(spread(dataf, key = cellID, value = y, fill = fill), 
                                list(y = as.name(names(dataf)[3]))))
  return(as.matrix(dataMatrix[,-1]))
}
