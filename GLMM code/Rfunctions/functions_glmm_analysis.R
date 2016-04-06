#### Functions for transforming data and analyzing model output

# Function to standardize covariates for occupancy models
standardize <- function(covar){
  mean.covar <- mean(covar, na.rm = TRUE)
  sd.covar <- sd(covar[!is.na(covar)])
  return((covar - mean.covar)/sd.covar)
}

# Make data matrices for JAGS GLMM with spatial random effect
makeEcoDataMatrix <- function(var, rep = "ymo", reVar = "cellID", data = countData, fill = NA){
  # This function takes a column (vector) of data and creates a matrix,
  # with spatial cells going across the top, and year-months going down
  # var = the variable (either explanatory or response)
  # rep = column identifying replicates
  # reVar = column identifying the random effect variable
  # data = the dataset to transform
  # fill = what should missing cells be filled with?
  dataf <- data[,c(rep, reVar, var)]
  dataMatrix <- eval(substitute(spread(dataf, key = key, value = y, fill = fill), 
                                list(key = as.name(names(dataf)[2]),
                                  y = as.name(names(dataf)[3]))))
  return(as.matrix(dataMatrix[,-1]))
}
