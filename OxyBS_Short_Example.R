#################################################
#  Quick example of OxyBS functionality 
#################################################
library(OxyBS)
data(OxyBSSampleData)

nSpecimens <- 30
nCpGs <- 30

# Calculate Total Signals
signalBS <- exampleMethBS+exampleUnmethBS
signalOxBS <- exampleMethOxBS+exampleUnmethOxBS

# Calculate Beta Values
betaBS <- exampleMethBS/signalBS
betaOxBS <- exampleMethOxBS/signalOxBS

# Create container for results
MethOxy <- array(NA,dim=c(nCpGs,nSpecimens,3))
dimnames(MethOxy) <- list(
  rownames(exampleMethBS)[1:nCpGs],
  colnames(exampleMethBS)[1:nSpecimens],
  c("C","5mC","5hmC"))

# Process results (one array at a time)
for(i in 1:nSpecimens){
  MethOxy[,i,] <- fitOxBS(betaBS[,i],betaOxBS[,i],signalBS[,i],signalOxBS[,i])
}

# Check that results sum to one
table(apply(MethOxy,1:2,sum))

# First specimen
MethOxy[,1,]

# Ranges of values
range(MethOxy[,,1])
range(MethOxy[,,2])
range(