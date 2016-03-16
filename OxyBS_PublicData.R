#################################################
#  EXAMPLE 1 - GSE63179 (4 Cerebellum Samples)
#################################################
# Necessary packages
library(OxyBS)
library(minfi) 
library(IlluminaHumanMethylation450kmanifest) 
library(IlluminaHumanMethylation450kanno.ilmn12.hg19) 
library(GEOquery)
data(OxyBSSampleData)

################################
# 1. Get GEO data
################################

# Load series and platform data from GEO
getGEOSuppFiles("GSE63179", baseDir = "/OxyBS/") 

# Get Phenotype Data
gse <- getGEO("GSE63179")
pheno <- pData(gse[[1]])

# Define a directory 
destination <- "/OxyBS/Data/"
dir.create(destination, recursive = T, showWarnings = F)

# Create a phenotype file for each GEO dataset for downstream analyses
write.table(pheno, "/OxyBS/Data/PhenoData_GSE63179.csv", sep = ",", row.names = T, col.names = NA)

# untar the folder downloaded into an external directory
untar("/OxyBS/GSE63179/GSE63179_RAW.tar", exdir = destination)

# gunzip the files in this newly untarred directory replacing all files
files <- list.files(destination)[grep(".idat.gz", list.files(destination))]
for(i in 1:length(files))
{
  gunzip(paste(destination, files[i], sep = ""), overwrite = T)
}

# Find file names based on GEO "GSM" individual sample identifier
idat_files <- list.files(file.path(destination), pattern="GSM")

# Define the basenames for each IDAT file
base <- unique(lapply(idat_files, function (x) {paste(unlist(strsplit(x, "_"))[1], unlist(strsplit(x, "_"))[2], unlist(strsplit(x, "_"))[3], sep = "_")}))

# Create phenotype basenames to match those listed with each IDAT file
pheno_base <- unique(lapply(as.character(pheno$supplementary_file), function (x) {paste(unlist(strsplit(unlist(strsplit(x, "/"))[11], "_"))[1], unlist(strsplit(x, "_"))[2], unlist(strsplit(x, "_"))[3], sep = "_")}))
pheno$pheno_base <- pheno_base

# Double check that the IDAT files match from the phenotype file and directory of downloaded IDAT files
length(intersect(pheno_base, base)) # 8 total array positions for 4 cerebellum samples

# Order the arrays such that samples match on BS and oxBS 
pheno_ordered <- pheno[order(pheno$title, decreasing=F), ]

# Create 'flags' for sample treatment with BS/oxBS
BS_indices <- grep("brain-BS", pheno_ordered$title)
oxBS_indices <- grep("oxBS", pheno_ordered$title)

# Name the array and slide position for BS treated samples
arraysBS <- intersect(pheno_ordered$pheno_base[BS_indices], base)

# Name the array and slide position for oxBS treated samples
arraysOxBS <- intersect(pheno_ordered$pheno_base[oxBS_indices], base)

####################################################
# 2. Read IDAT files and preprocess Infinium arrays
####################################################

# Input array list to minfi read.450k array function
dataDir <- "/OxyBS/Data/"
datListBS <- read.450k(file.path(dataDir, arraysBS),verbose=TRUE)
datListOxBS <- read.450k(file.path(dataDir, arraysOxBS),verbose=TRUE)

# Create new preprocessFunNorm() function
preprocessFunnormRedGreen <- function (rgSet, nPCs = 2, sex = NULL, verbose = TRUE)
{
  minfi:::.isRG(rgSet)
  rgSet <- updateObject(rgSet)
  if (verbose)
    cat("[preprocessFunnorm] Mapping to genome\n")
  gmSet <- mapToGenome(rgSet)
  subverbose <- max(as.integer(verbose) - 1L, 0)
  if (verbose)
    cat("[preprocessFunnorm] Quantile extraction\n")
  extractedData <- minfi:::.extractFromRGSet450k(rgSet)
  if (is.null(sex)) {
    gmSet <- addSex(gmSet, getSex(gmSet, cutoff = -3))
    sex <- rep(1L, length(gmSet$predictedSex))
    sex[gmSet$predictedSex == "F"] <- 2L
  }
  rm(rgSet)
  if (verbose)
    cat("[preprocessFunnorm] Normalization\n")
  CN <- getCN(gmSet)
  minfi:::.normalizeFunnorm450k(object = gmSet, extractedData = extractedData,
                                sex = sex, nPCs = nPCs, verbose = subverbose)
  
}

# Preprocess using functional normalization separately for both BS and oxBS arrays
rgBS <-  preprocessFunnormRedGreen(datListBS)
rgOxBS <-  preprocessFunnormRedGreen(datListOxBS)


####################################################
# 3. Define normalized signals and calculate Betas
####################################################

# Methylated signals from the BS and oxBS arrays
methBS <- assay(rgBS,"Meth")
methOxBS <- assay(rgOxBS,"Meth")
# Unmethylated signals from the BS and oxBS arrays
unmethBS <- assay(rgBS,"Unmeth")
unmethOxBS <- assay(rgOxBS,"Unmeth")

# Calculate Total Signals
signalBS <- methBS+unmethBS
signalOxBS <- methOxBS+unmethOxBS

# Calculate Beta Values
betaBS <- methBS/signalBS
betaOxBS <- methOxBS/signalOxBS

####################################################
# 4. Apply fitOxBS function to preprocessed values
####################################################

# Select the number of CpGs and Subjects to which the method will be applied 
nCpGs <- dim(unmethOxBS)[1]
nSpecimens <- dim(unmethOxBS)[2]

# Create container for the OxyBS results
MethOxy <- array(NA,dim=c(nCpGs,nSpecimens,3))
dimnames(MethOxy) <- list(
  rownames(methBS)[1:nCpGs],
  colnames(methBS)[1:nSpecimens], c("C","5mC","5hmC"))

# Process results (one array at a time, slow)
# for(i in 1:nSpecimens){
#	MethOxy[,i,] <-fitOxBS(betaBS[,i],betaOxBS[,i],signalBS[,i],signalOxBS[,i])
#	}

# Process the results using parallelization and implementation of foreach function
# Import parallelization packages
library(foreach)
library(doParallel)

# Set-up parallel backend to use all but one available processors
no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)
registerDoParallel(cl)

# Calculate the time of funciton
strt <- Sys.time()

# Parallelized loop for 4 Cerebellum samples across all CpGs on the 450K array 
MethOxy[,1:4,] <- foreach(i = 1:nSpecimens, .combine=cbind, .packages='OxyBS') %dopar% {
  fitOxBS(betaBS[,i],betaOxBS[,i],signalBS[,i],signalOxBS[,i])
}

print(Sys.time()-strt) # Time difference of 31.42 mins on 16GB local machine
stopCluster(cl)

####################################################
# 5. Output - Examples and visualizations
####################################################

# Output is a 3D array for 5C, 5mC, and 5hmC levels
# First specimen
MethOxy[,1,]
# 5-hydroxymethylcytosine values for first few CpGs
head(MethOxy[,,3])

# Ranges for each cytosine modification
range(MethOxy[,,1]) # 5C
range(MethOxy[,,2]) # 5mC
range(MethOxy[,,3]) # 5hmC

# Check that results sum to one
table(apply(MethOxy,1:2,sum))

# Some NaNs may be produced when using OxyBS
any(is.na(MethOxy)) # TRUE

# Compare OxyBS with naive approach
MethOxy0 <- MethOxy
MethOxy0[,,1] <- 1-betaBS
MethOxy0[,,2] <- betaOxBS
MethOxy0[,,3] <- betaBS-betaOxBS