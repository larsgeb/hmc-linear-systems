#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

cat('Estimation of (m)ESS. \r\n\r\n')

# Install packages (do this only once, maybe easier from R-console)
install.packages('coda')
install.packages('mcmcse')

# Load libraries
library(coda)
cat('Loading mcmcse library... \r\n')
library(mcmcse)

# Read samples and drop last column (containing misfits)
cat('Loading samples into dataframe ... \r\n', 'Samples file: ',args[1], '\r\n\r\n')
dataframe_markovchain <- read.table(args[1], header = FALSE)
dataframe_markovchain[ , c(length( names( dataframe_markovchain ) ) -1 ) ] <-NULL

# Create mcmc object for coda package
mcmcObject_markovchain <- mcmc(data=dataframe_markovchain)

# Compute (Multivariate) Effective Sample Sizes
mess <- multiESS(dataframe_markovchain) # The multivariate ESS using coda
ess <- effectiveSize(mcmcObject_markovchain) # The per-parameter ESS using mcmcse

# Compute needed samples
neededSamples <- minESS(ncol(dataframe_markovchain))

# Output results
cat("Number of effective samples (according to coda): \r\n", ess, "\r\n\r\n")
cat("Number of effective multivariate samples (according to mcmcse): \r\n", mess, "\r\n\r\n")
cat("Number of effective multivariate samples needed (according to mcmcse):\r\n ", neededSamples, "\r\n\r\n")
if (mess >= neededSamples){
  cat("Properties are estimated with confidence!\r\n")
}else{
  cat("Chain should be run longer for better confidence!\r\n")
}

