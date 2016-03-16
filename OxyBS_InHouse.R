#################################################
#  EXAMPLE 2 -  GSE73895 (30 Glioblastoma Samples)
#  In-house data set
#################################################
# Necessary packages
library(OxyBS)
library(minfi) 
library(IlluminaHumanMethylation450kmanifest) 
library(IlluminaHumanMethylation450kanno.ilmn12.hg19) 
library(GEOquery)

################################
# 1. Process .IDAT files and match with identifier with phenotype file
################################
# Location of the idat files
dataDir <- "GSE73895/IDAT_Files"
setwd(dataDir)
# Retrieve the name of the microarray chips used in this study
chips <- list.files(file.path(dataDir), pattern="^[[:digit:]]*")
nChips <- length(chips)

#Generate array list that contains named chip location
arrayList <- list()
for(i in 1:nChips){
  ff <- grep("Grn.idat$",list.files(path=chips[i],pattern=chips[i]),value=TRUE)
  arrayList[[i]] <- paste(chips[i],gsub("[_]Grn.idat","",ff),sep="/")
}

#Name each list with the chip names
names(arrayList) <- chips

# The meta file contains the covariate information along with sample bisulfite treatment information
meta <- read.delim("metaData_GSE73895.txt", sep="\t", head=TRUE, 
                   stringsAsFactors=FALSE)
meta$ArrayID <- with(meta, paste(Beadchip,Position,sep="_"))
meta$Label.ID.[meta$OX_BS==1] <- gsub("B$","",meta$ID[meta$OX_BS==1])
meta$Label.ID.[meta$OX_BS==0] <- gsub("A$","",meta$ID[meta$OX_BS==0])
head(meta)

# Ensure that paired samples are matched with one another
with(meta, all(Label.ID.==gsub("[A|B]$","",ID))) # Check consistency of IDs

# Flag each Labelled ID for matched BS and oxBS treatment pair
mrg <- merge(
  data.frame(LabelID=meta$Label.ID.[meta$OX_BS==0],
             ixBS=as.numeric(rownames(meta)[meta$OX_BS==0]),
             stringsAsFactors=FALSE),
  data.frame(LabelID=meta$Label.ID.[meta$OX_BS==1],
             ixOxBS=as.numeric(rownames(meta)[meta$OX_BS==1]),
             stringsAsFactors=FALSE)
)

# Create index to the row location of matched BS-oxBS treatment
meta$Pair[mrg$ixOxBS] <- mrg$ixBS 
meta$Pair[mrg$ixBS] <- mrg$ixOxBS 

# New rownames are the location on the given array
rownames(meta) <- meta$ArrayID

# Name the array and slide position for BS treated samples
arraysBS <- intersect(unlist(arrayList),
                      with(meta, paste(Beadchip,ArrayID,sep="/")[OX_BS==0]))

# Name the array and slide position for oxBS treated samples
arraysOxBS <- intersect(unlist(arrayList),
                        with(meta, paste(Beadchip,ArrayID,sep="/")[OX_BS==1]))

####################################################
# 2. Read IDAT files and preprocess Infinium arrays
####################################################
datListBS <- read.450k(file.path(dataDir, arraysBS),verbose=TRUE)
datListOxBS <- read.450k(file.path(dataDir, arraysOxBS),verbose=TRUE)

# Preprocess using functional normalization separately for both BS and oxBS arrays
rgBS <-  preprocessFunnormRedGreen(datListBS,
                                   sex=1*(meta[gsub("^[[:digit:]]*/","",arraysBS),"GENDER"]=="M"))
rgOxBS <-  preprocessFunnormRedGreen(datListOxBS,
                                    sex=1*(meta[gsub("^[[:digit:]]*/","",arraysOxBS),"GENDER"]=="M"))

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


# Process results (one array at a time, VERY slow for >10 samples. Consider parallelization)
for(i in 1:nSpecimens){
  MethOxy[,i,] <-fitOxBS(betaBS[,i],betaOxBS[,i],signalBS[,i],signalOxBS[,i])
}


