# OxyBS
Vignette for data processing and use of the fitOxBS function from the OxyBS R package

OxyBS is a R package that provides a maximum-likelihood approach to the estimation of 5mC, 5hmC, and 5C in the analysis of paired oxidative bisulfite and bisulfite Infinium DNA methylation microarrays. Supported assays include the Infinium 450K, InfiniumEpic (850K), and sequencing data with read counts for the four input data types (BS-methylated, BS-unmethylated, oxBS-methylated, and oxBS-unmethylated).

# Contents
* Getting Started
  * Datasets
  * GEO download and processing of in-house dataset
* Example use of fitOxBS
* Advanced Usage 
  * Visualization

# 1. Getting Started
## Datasets
In this vignette, we describe how to process data generated by the Infinium 450k methylation array using the Cambridge Epigenetix TrueMethyl kit. We will display the workhorse function of the OxyBS package (i.e., fitOxBS) using a study of publicly available data on GEO as well as “in-house” dataset of thirty glioblastomas. The publicly available data comprises of eight DNA samples from the same DNA source that have each been treated with oxBS-BS and hybridized to the Infinium 450K array. This dataset is currently published (Field et al 2015 PLoS One) and can be downloaded from the Gene Expression Omnibus website (http://www.ncbi.nlm.nih.gov/geo/, GSE63179) or accessed with the R-package GEOquery. The in-house dataset is composed of sixty DNA samples from thirty glioblastoma samples treated with BS and oxBS (True Methyl Array, GSE73895). This data set was under review at the time of this publication.

## Data Download
In the example of publicly available data, it will be necessary to download the .IDAT (IDAT = intensity data) files for each of the study samples. The .IDAT files contain the 5-methylcytosine and 5-hydroxymethylcytosine data. Here, we describe use of the GEOquery package so that data procurement from GEO can be performed exclusively using R. 

Alternatively, if the user is working with an in-house dataset (i.e., experiments done by the laboratory group) we have provided example data sheets and utilities to assist users unfamiliar with the minfi preprocessing pipeline. 

# 2. Example Estimation
OxyBS provides functions to estimate cytosine modifications when BS-methylated, BS-unmethylated, oxBS-methylated, and oxBS-unmethylated signals are provided. OxyBS doesn’t not currently provide tools for comprehensive analysis of cytosine modifications nor does the package provide methods to directly import data. Instead, this manual presents approaches to processing raw .IDAT files when available on GEO as well as for experiments produced in-house. Figure 1 provides a general flow diagram to assist with the processing of .IDAT files and to apply the maximum likelihood estimation approach (i.e., fitOxBS function).

### 2A. Example with Publicly Available Data (setting up the analysis environment)
Step 0: Before obtaining estimates from the paired oxBS-BS data, we need to prepare the appropriate input files and have essential libraries loaded. To this end, the OxyBS, minfi (plus associated annotation packages), and GEOquery packages need to be installed and loaded in the user’s workspace. 

Step 1: For this example, we will use a small dataset of 4 samples in GSE63179. It is important to note that the user will have to adapt the directories corresponding to their own file system and operating system. For the present example, the directory “/OxyBS/” will contain the downloaded GEO files.

Step 2: Once the .IDAT files have  been procured from GEO and we have separated the arrays into groups either treated with sodium bisulfite only (i.e. BS) or a selective oxidation step before sodium bisulfite (i.e., oxBS). Next, we apply the ‘read.450k’ function from the minfi package and preprocess the data with functional normalization. We refer users’ to the minfi help page for issues/questions with this step.

Step 3: To return methylation signals and beta values for the BS-treated and oxBS-treated arrays we need to access the preprocessed samples (i.e., rgBS object) for methylated/unmethylated signals, tally the total signals, and calculate the Beta values for both cytosine modifications (i.e., object 'betaBS').

Step 4: Finally, to estimate cytosine modification from paired BS-oxBS arrays it is necessary to apply the “fitOxBS” function to each CpG for each sample. This process is computationally demanding making use of this algorithm for a large number of samples inefficient when the function is fit sequentially. Therefore, it is recommended the user implement parallelization of their code on their local computer or computing cluster.

Step 5: The output from the ‘fitOxBS’ function can be placed in a three-dimensional array with the structure of CpGs x Subjects x cytosine modifications (5C, 5mC, 5hmC).

### 2B. Example with an in-house dataset (Setting up the analysis environment)
To apply OxyBS to a dataset produced in-house requires a phenotype and chip information meta file in addition to .IDAT files. Below we present an example of such a meta-file and progress with the analysis (Steps: 2-4) in the same manner as was presented for section 2A.

# 3. Visualization
To quickly assess the frequency and magnitude of negative 5hmC values when a naive (BS-oxBS=5hmC) estimation
