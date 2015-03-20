args <- commandArgs(TRUE)

options(stringsAsFactors=FALSE)


inFile <- args[1]  #TITAN
inFile2 <- args[2] #HMMcopy
inFile3 <- args[3] #APOLLOH 
outFile <- args[4]

## load in spike in results ##
loh <- read.delim(inFile,header=T,sep="\t")
chrs <- c(1,2,9,18)
loh <- loh[loh$Chr %in% chrs,1:12]

## load in original APOLLOH1 results ##
cna2 <- read.delim(inFile2,header=T,sep="\t")
loh2 <- read.delim(inFile3,header=F,sep="\t")
chrs <- c(1,2,8,9,16,18)
loh2 <- loh2[loh2[,1] %in% chrs,]
cna2 <- cna2[cna2[,2] %in% chrs,]

## load normal content of spikein sample ##
paramFile <- gsub("loh.txt","params.txt",inFile)
params <- read.delim(paramFile,header=F,sep="\t")
norm <- as.numeric(params[params[,1]=="Normal contamination estimate:",2])

png(outFile,width=1000,height=1500,res=100)
#par(mfrow=c(3,2))
layout(matrix(c(1,2,3,4,5),5,1,byrow=F))

source('/share/data/gha/software/LOH-HMM_0.1.0/plotting/plotSNPallelicRatio.R')
source('/share/data/gha/software/HMMcopy/scripts/plotCNlogRByChr.R')
plotCNlogRByChr(cna2, chr=NULL, geneAnnot=genes, spacing=4, cytoBand=F, alphaVal=1)
plotSNPallelicRatio(loh2, chr=NULL, geneAnnot=genes, spacing=4, cytoBand=F, alphaVal=1)

#plot(0,type="n",bty="n",xlab="",ylab="",xaxt="n",yaxt="n")

source('/share/data/gha/software/APOLLOH_2.0.1/scripts/spikeInExpt/plotting.R')
plotCNlogRByChr(loh, segments=NULL, chr=NULL, geneAnnot=NULL, ylim=c(-4,6), spacing=4, alphaVal=1, xlim=NULL, cex=0.5)
plotAllelicRatio(loh, chr=NULL, geneAnnot=NULL, spacing=4, ylim=c(0,1), xlim=NULL, cex=0.5)
plotClonalFrequency(loh, chr=NULL, normal=norm, geneAnnot=NULL, spacing=4, ylim=c(0,1),xlim=NULL, cex=0.5)

dev.off()