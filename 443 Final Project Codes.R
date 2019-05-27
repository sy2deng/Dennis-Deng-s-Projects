# 443 Assignment

# This code makes use of the package ggplot2
# If this is not installed on your machine you need to run
# install.packages("ggplot2")
library(ggplot2)

# a)

# Number of times process is repeated
numA = 120 

# newARMAa: Runs a new ARMA(p,q) given n, ar, and ma. 
# n = number of ARMA(p,q) steps run
# ar = vector with coefficients for the AR portion
# ma = vector with coefficients for the MA portion
newARMAa <- function(n,ar,ma)
{
  sim.arma <- arima.sim(n=n, list(ar=ar,ma=ma), sd = sqrt(1))
  return(sim.arma)
}

#AICmatrixA: uses min/max inputs of p & q and creates numb matrices of AICs from the models used to fit 
#           ARMA(p,q). Finds the ARMA(p,q) model with the minimum AIC and uses that to predict
#           10 steps ahead using the predict function. Outputs list of matricies, list of (p,q) values with
#           the lowest AICs, forecasted values, and generated values.
# n = number of ARMA(p,q) steps run
# numb = number of times process is repeated
# pstart = minimum p used in ARMA(p,q) models used to fit
# pfinish = maximum p used in ARMA(p,q) models used to fit
# qstart = minimum q used in ARMA(p,q) models used to fit
# qfinish = maximum q used in ARMA(p,q) models used to fit
AICmatrixA <- function(n,numb,pstart,pfinish,qstart,qfinish)
{
  
  output <- list(); simData <- list(); mat <- list(); pqValues <- list(); forecast10 <- list()
  
  for (p in 1:numb)
  {
    # Generating simulated data from model 
    simData[[p]] <- newARMAa(n,c(0.5),c(1))
    mat[[p]] <- matrix(nrow = pfinish + 1,ncol = qfinish + 1)
    dimnames(mat[[p]]) <- list(c(pstart:pfinish),c(qstart:qfinish))
    
    for (i in pstart:pfinish)
    {
      for (l in qstart:qfinish)
      {
        # We pretend that we don't know p and q so we look at the matrix.
        k <- arima(simData[[p]][1:(n-10)], order=c(i,0,l), include.mean = F,method="ML") 
        mat[[p]][i+1,l+1] <- k$aic
      }
    }
    
    # Getting output of min aic values and the corresponding matrices
    output[[p]] <- mat[[p]]
    pqValues[p] <- paste(as.character(c(which(mat[[p]] == min(mat[[p]]), arr.ind = TRUE)) - 1),collapse="")
    pval <- which(mat[[p]] == min(mat[[p]]), arr.ind = TRUE)[1] - 1
    qval <- which(mat[[p]] == min(mat[[p]]), arr.ind = TRUE)[2] - 1
    
    # Predicting simData[[p]] for 10 steps ahead
    forecast10[p] <- predict(arima(simData[[p]][1:(n-10)], order=c(pval,0,qval),method="ML"),n.ahead=10)
  }
  
  out <- list(output,pqValues,forecast10,simData)
  return(out)
}

# Running the function to get output
outA <- AICmatrixA(110,numA,0,2,0,2)
# Extracting Matrices with AIC values
matricesA <- outA[[1]]
# Extracting the p & q values for the model 
pqValsA <- unlist(outA[2])
# Extracting the forecasts 
f10A <- outA[[3]]
# Extracting the simulated data
sDataA <- outA[[4]]

# Renaming the matricesA list with the corresponding p & q values
names(matricesA) <- pqValsA

# Creating table to see how many times it picked the correct p & q values
freqA <- table(pqValsA == "11")
# Getting the proportion 
proportionA <- freqA[[2]]/(freqA[[1]] + freqA[[2]])

# This part is getting the squarred differences between the forecasted and actual data
d10A <- list()

for (i in 1:numA)
{
  d10A[[i]] <- (sDataA[[i]][101:110] - f10A[[i]])^2
}

# This part is adding up all the values
sum10A <- c(0,0,0,0,0,0,0,0,0,0)

for (i in 1:numA)
{
  sum10A <- sum10A + d10A[[i]]
}

# Dividing to get the MSE values
mse10A <- sum10A/numA

# Getting the MSE's for 1,2,5,10 steps ahead forecast
mylistA <- mse10A[c(1,2,5,10)]

mylistA
table(pqValsA)
proportionA

# plot
# Makes the pq plots 
pqValsA.df <- as.data.frame(table(pqValsA))
pqPlot11 <- ggplot(data=pqValsA.df, aes(x=pqValsA, y=Freq)) + geom_bar(stat="identity")
pqPlot11 <- pqPlot11 + labs(title = "Suggested pq Values from ARMA(1,1) Data",x="pq Values",y="Frequency")

# Makes the MSE plots
mseNames <- c(1,2,5,10)
mseA.df <- data.frame(mseNames,mylistA)
msePlot11 <- ggplot(data = mseA.df, aes(x=factor(mseNames), y=mylistA)) + geom_bar(stat="identity")
msePlot11 <- msePlot11 + labs(title = "MSE Values from ARMA(1,1) Data",x="Steps Ahead",y="MSE")

# b)

# Number of times process is repeated
numB = 120 

# newARMAb: Runs a new ARMA(p,q) given n, ar, and ma. 
# n = number of ARMA(p,q) steps run
# ar = vector with coefficients for the AR portion
# ma = vector with coefficients for the MA portion
newARMAb <- function(n,ar)
{
  sim.arma <- arima.sim(n=n, model=list(ar=ar), sd = sqrt(1))
  return(sim.arma)
}

#AICmatrixB: uses min/max inputs of p & q and creates numb matrices of AICs from the models used to fit 
#           ARMA(p,q). Finds the ARMA(p,q) model with the minimum AIC and uses that to predict
#           10 steps ahead using the predict function. Outputs list of matricies, list of (p,q) values with
#           the lowest AICs, forecasted values, and generated values.
# n = number of ARMA(p,q) steps run
# numb = number of times process is repeated
# pstart = minimum p used in ARMA(p,q) models used to fit
# pfinish = maximum p used in ARMA(p,q) models used to fit
# qstart = minimum q used in ARMA(p,q) models used to fit
# qfinish = maximum q used in ARMA(p,q) models used to fit
AICmatrixB <- function(n,numb,pstart,pfinish,qstart,qfinish)
{
  
  output <- list(); simData <- list(); mat <- list(); pqValues <- list(); forecast10 <- list()
  
  for (p in 1:numb)
  {
    # Generating simulated data from model 
    simData[[p]] <- newARMAb(n,c(0.3))
    mat[[p]] <- matrix(nrow = pfinish + 1,ncol = qfinish + 1)
    dimnames(mat[[p]]) <- list(c(pstart:pfinish),c(qstart:qfinish))
    
    for (i in pstart:pfinish)
    {
      for (l in qstart:qfinish)
      {
        # We pretend that we don't know p and q so we look at the matrix.
        k <- arima(simData[[p]][1:(n-10)], order=c(i,0,l), include.mean = F,method="ML") 
        mat[[p]][i+1,l+1] <- k$aic
      }
    }
    
    # Getting output of min aic values and the corresponding matrices
    output[[p]] <- mat[[p]]
    pqValues[p] <- paste(as.character(c(which(mat[[p]] == min(mat[[p]]), arr.ind = TRUE)) - 1),collapse="")
    pval <- which(mat[[p]] == min(mat[[p]]), arr.ind = TRUE)[1] - 1
    qval <- which(mat[[p]] == min(mat[[p]]), arr.ind = TRUE)[2] - 1
    
    # Predicting simData[[p]] for 10 steps ahead
    forecast10[p] <- predict(arima(simData[[p]][1:(n-10)], order=c(pval,0,qval), method="ML"),n.ahead=10)
  }
  
  out <- list(output,pqValues,forecast10,simData)
  return(out)
}

# Running the function to get output
outB <- AICmatrixB(30,numB,0,1,0,1)
# Extracting Matrices with AIC values
matricesB <- outB[[1]]
# Extracting the p & q values for the model 
pqValsB <- unlist(outB[2])
# Extracting the forecasts 
f10B <- outB[[3]]
# Extracting the simulated data
sDataB <- outB[[4]]

# Renaming the matricesB list with the corresponding p & q values
names(matricesB) <- pqValsB

# Creating table to see how many times it picked the correct p & q values
freqB <- table(pqValsB == "10")
# Getting the proportion 
proportionB <- freqB[[2]]/(freqB[[1]] + freqB[[2]])

# This part is getting the squarred differences between the forecasted and actual data
d10B <- list()

for (i in 1:numB)
{
  d10B[[i]] <- (sDataB[[i]][21:30] - f10B[[i]])^2
}

# This part is adding up all the values
sum10B <- c(0,0,0,0,0,0,0,0,0,0)

for (i in 1:numB)
{
  sum10B <- sum10B + d10B[[i]]
}

# Dividing to get the MSE values
mse10B <- sum10B/numB

# Getting the MSE's for 1,2,5,10 steps ahead forecast
mylistB <- mse10B[c(1,2,5,10)]

mylistB
table(pqValsB)
proportionB

# plots
# Makes the pq plots 
pqValsB.df <- as.data.frame(table(pqValsB))
pqPlot10 <- ggplot(data=pqValsB.df, aes(x=pqValsB, y=Freq)) + geom_bar(stat="identity")
pqPlot10 <- pqPlot10 + labs(title = "Suggested pq Values from ARMA(1,0) Data",x="pq Values",y="Frequency")

# Makes the MSE plots
mseNames <- c(1,2,5,10)
mseB.df <- data.frame(mseNames,mylistB)
msePlot10 <- ggplot(data = mseB.df, aes(x=factor(mseNames), y=mylistB)) + geom_bar(stat="identity")
msePlot10 <- msePlot10 + labs(title = "MSE Values from ARMA(1,0) Data",x="Steps Ahead",y="MSE")

# c)

numC = 120 #Number of times process is repeated

# arch.sim: ARCH(1) function copied from slides
# n = number of ARCH(1) steps run
# omega = alpha 0 value for ARCH(1)
# alpha1 = alpha 1 value for ARCH(1)
# sigma = sigma value for ARCH(1)
arch.sim <- function(n, omega, alpha1, sigma)
{
  out <- sqrt(omega)
  for(i in 2:n)
  {
    out[i] <- sqrt(omega + alpha1*out[i-1]^2)*rnorm(1, sd=sigma)
  }
  out
}

# newARCH: Runs a new ARCH(1) given n, omega, and alpha1. We assume sigma = 1
# n = number of ARCH(1) steps run
# om = alpha 0 value for ARCH(1)
# alp = alpha 1 value for ARCH(1)
newARCH <- function(n,om,alp)
{
  sim.arch <- arch.sim(n=n, omega = om, alpha1 = alp, sigma = sqrt(1))
  return(sim.arch)
}

#AICmatrix: uses min/max inputs of p & q and creates numb matrices of AICs from the models used to fit 
#           ARCH(1) to ARMA(p,q). Finds the ARMA(p,q) model with the minimum AIC and uses that to predict
#           10 steps ahead using the predict function. Outputs list of matricies, list of (p,q) values with
#           the lowest AICs, forecasted values, and generated values.
# n = number of ARCH(1) steps run
# numb = number of times process is repeated
# pstart = minimum p used in ARMA models used to fit
# pfinish = maximum p used in ARMA models used to fit
# qstart = minimum q used in ARMA models used to fit
# qfinish = maximum q used in ARMA models used to fit
AICmatrixC <- function(n,numb,pstart,pfinish,qstart,qfinish)
{
  output <- list(); simData <- list(); mat <- list(); pqValues <- list(); forecast10 <- list()
  
  for (p in 1:numb)
  {
    # Generating simulated data from model 
    simData[[p]] <- newARCH(n,0.3,0.65)
    mat[[p]] <- matrix(nrow = pfinish + 1,ncol = qfinish + 1)
    dimnames(mat[[p]]) <- list(c(pstart:pfinish),c(qstart:qfinish))
    
    for (i in pstart:pfinish)
    {
      for (l in qstart:qfinish)
      {
        # We pretend that we don't know p and q so we look at the matrix.
        k <- arima(simData[[p]][1:(n-10)], order=c(i,0,l), include.mean = F, method="ML") 
        mat[[p]][i+1,l+1] <- k$aic
      }
    }
    
    # Getting output of min aic values and the corresponding matrices
    output[[p]] <- mat[[p]]
    pqValues[p] <- paste(as.character(c(which(mat[[p]] == min(mat[[p]]), arr.ind = TRUE)) - 1),collapse="")
    pval <- which(mat[[p]] == min(mat[[p]]), arr.ind = TRUE)[1] - 1
    qval <- which(mat[[p]] == min(mat[[p]]), arr.ind = TRUE)[2] - 1
    
    # Predicting simData[[p]] for 10 steps ahead
    forecast10[p] <- predict(arima(simData[[p]][1:(n-10)], order=c(pval,0,qval), method="ML"),n.ahead=10)
  }
  
  out <- list(output,pqValues,forecast10,simData)
  return(out)
}

# Runs AICmatrixC with 110 data points and with p & q values between 0-2
outC <- AICmatrixC(110,numC,0,2,0,2)
matricesC <- outC[[1]] # matricesC of AICs
pqValsC <- unlist(outC[2]) # pq values with lowest AICs
f10C <- outC[[3]] # forecasted values
sDataC <- outC[[4]] # generated ARCH values
names(matricesC) <- pqValsC

d10C <- list()

# Calculates SEs for each prediction
for (i in 1:numC)
{
  d10C[[i]] <- (sDataC[[i]][101:110] - f10C[[i]])^2
}

sum10C <- c(0,0,0,0,0,0,0,0,0,0)

# Sums SEs at each step
for (i in 1:numC)
{
  sum10C <- sum10C + d10C[[i]]
}

# calculates MSE
mse10C <- sum10C/numC

# MSEs for steps 1,2,5,10 isolated
mylistC <- mse10C[c(1,2,5,10)]

# Output
mylistC
table(pqValsC)

# plots
# Makes the pq plots 
pqValsC.df <- as.data.frame(table(pqValsC))
pqPlot1 <- ggplot(data=pqValsC.df, aes(x=pqValsC, y=Freq)) + geom_bar(stat="identity")
pqPlot1 <- pqPlot1 + labs(title = "Suggested pq Values from ARCH(1) Data",x="pq Values",y="Frequency")

# Makes the MSE plots
mseNames <- c(1,2,5,10)
mseC.df <- data.frame(mseNames,mylistC)
msePlot1 <- ggplot(data = mseC.df, aes(x=factor(mseNames), y=mylistC)) + geom_bar(stat="identity")
msePlot1 <- msePlot1 + labs(title = "MSE Values from ARCH(1) Data",x="Steps Ahead",y="MSE")
