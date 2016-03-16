#################################################
#  EXAMPLE 3 - GSE63179 (4 Cerebellum Samples)
#  Example visualizations of naive estimations of 5hmC
#################################################
# Files processed in OxyBS pipeline
load("annotation-amended.RData") # 'longAnnot' object
load("OxyBrainSignals-FunNorm.RData")
load("OxyBrainMethOxy-FunNorm.RData")

# Keep the probes that are not cross-hybridizing, SNP-related, poor performing as indicated in Chen et al 2013 Epigenetics 
geoUCSC <- gsub("^[N|S]_","", longAnnot$RelToIslandUCSC)
geoUCSC.good.I <- geoUCSC[lmIndex[["I:0:0:0:0"]]]
geoUCSC.good.II <- geoUCSC[lmIndex[["II:0:0:0:0"]]]

# Stratify CpG probes into Illumina Type I and Type II
setI.strat <- split(lmIndex[["I:0:0:0:0"]], geoUCSC.good.I)
setII.strat <- split(lmIndex[["II:0:0:0:0"]], geoUCSC.good.II)

# Use the signals produced in the OxyBS_PublicData.R file
signalBS <- unmethBS+ methBS
signalOxBS <- unmethOxBS+ methOxBS

# Calculate the beta values from BS treated and oxBS treated DNA
betaBS <- methBS/signalBS # All Modifications
betaOxBS <- methOxBS/signalOxBS # Naive estimated 'True' 5mC 
betaDiff <- betaBS-betaOxBS # Naive estimated 5hmC levels

####################################################
# Generate summary statistics to assess probes with values less than 0 
# for 5hmC and define the magnitude of loss 
####################################################
# Define % of CpGs < 0 for 5hmC and quartiles for negative 5hmC values 
# Type I Infinium probes
stats.I <- sapply(setI.strat, function(u){
  v <- as.vector(betaDiff[u,])
  flag <- (v < 0)
  c(mean(flag),quantile(v[flag],c(0,.25,0.5,0.75,1)))
})

# Type II Infinium probes
stats.II <- sapply(setII.strat, function(u){
  v <- as.vector(betaDiff[u,])
  flag <- (v < 0)
  c(mean(flag, na.rm=TRUE),quantile(v[flag],c(0,.25,0.5,0.75,1), na.rm=TRUE))
})

###########################
# Produce sample figure 
###########################
# Positions of the graphs
pos <- c(1:4,6:9)
# Width of the bar graphs on the x-axis
xwd <- 0.35
# Color for separate genomic regions (i.e., CpG Island, etc)
clr <- c(topo.colors(4)[4:1],topo.colors(4)[4:1])

# Separate figure into two portions
layout(1:2, height=c(4,5))
# Convert fraction of negative 5hmC values to percentage
pcts <- 100*cbind(stats.I[1,c(1,4,3,2)],stats.II[1,c(1,4,3,2)])

# Top Portion: % of CpGs per probe-type less than 0 5hmC
par(mar=c(0.1, 4.1, 2.1, 2.1))
plot(c(0,9.96), c(1,max(pcts)), type="n", xlab="", bty="n",
     ylab="% less than 0",xaxt="n",cex.axis=0.9, las=2)
for(i in 1:8) {
  rect(pos[i]-xwd,0,pos[i]+xwd,pcts[i],col=clr[i])
}
lines(pos[c(1,4)]+2*xwd*c(-1,1), c(0,0))
lines(pos[c(5,8)]+2*xwd*c(-1,1), c(0,0))

# Bottom Portion: Naive 5hmC %
par(mar=c(5.1, 4.1, 1.1, 2.1))
pcts2 <- 100*cbind(stats.I[-1,c(1,4,3,2)],stats.II[-1,c(1,4,3,2)])
plot(c(0,9.96), range(pcts2), type="n", xlab="", ylab="Naive 5hmC (%)", xaxt="n", 
     cex.axis=0.9, las=2)
for(i in 1:8) {
  rect(pos[i]-xwd,pcts2[2,i],pos[i]+xwd,pcts2[4,i],col=clr[i])
  lines(pos[i]+xwd*c(-1,1), rep(pcts2[3,i],2),lwd=2)
  lines(rep(pos[i],2),pcts2[4:5,i],lwd=2)
  lines(rep(pos[i],2),pcts2[1:2,i],lwd=2)
}
axis(1, pos, gsub("Open","", colnames(pcts2)), cex.axis=0.8)
axis(1, mean(pos[1:4]), "Type I", line=1.25, tick=FALSE, cex.axis=0.8)
axis(1, mean(pos[5:8]), "Type II", line=1.25, tick=FALSE, cex.axis=0.8)


