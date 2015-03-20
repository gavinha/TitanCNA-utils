library(IRanges)

args <- commandArgs(TRUE)

options(stringsAsFactors=FALSE)


inFile <- args[1]  #TITAN
inFile2 <- args[2] #HMMcopy
inFile3 <- args[3] #APOLLOH 
regionsFile <- args[4] #spike-in truth file
truthFile1 <- args[5] #SpikeIn_TumCounts_LogR.txt
truthFile2 <- args[6] #SpikeIn_truth_posns.txt
outDir <- args[7] #output dir


## load in spike in results ##
loh <- read.delim(inFile,header=T,sep="\t")
chrs <- c(1,2,9,18)
loh <- loh[loh$Chr %in% chrs,1:12]

## load in regions ##
regions <- read.delim(regionsFile,header=T,sep="\t")
regions <- cbind(1:nrow(regions),regions[,c("chr","start","end")])

## load truth positions ##
truthData <- read.delim(truthFile1,header=T,sep="\t")
truthData <- truthData[truthData[,1] %in% chrs,]

truthPosn <- read.delim(truthFile2,header=T,sep="\t")
## compute clonal cluster, cellular prevalence, apollohCalls
cellPrev <- rep(NA,nrow(truthPosn))
clust <- rep(NA,nrow(truthPosn))
calls <- rep("HET",nrow(truthPosn))
cellPrev[grepl("tum100",truthPosn$id)] <- 0.001; clust[grepl("tum100",truthPosn$id)] <- 1
cellPrev[grepl("tum80-norm20",truthPosn$id)] <- 0.20; clust[grepl("tum80-norm20",truthPosn$id)] <- 2
cellPrev[grepl("tum60-norm40",truthPosn$id)] <- 0.40; clust[grepl("tum60-norm40",truthPosn$id)] <- 3
calls <- truthPosn$type.type; calls[calls=="DEL"] <- "DLOH"; calls[calls=="AMP"] <- "ASCNA"
truthPosn <- cbind(truthPosn,APOLLOHcall=calls,ClonalCluster=clust,CellularFrequency=cellPrev)

## setup truth prevalence ##
truthDataIR <- RangedData(ranges=IRanges(start=truthData$start,end=truthData$start),space=truthData$chr)
truthPosnIR <- RangedData(ranges=IRanges(start=truthPosn$start,end=truthPosn$start),space=truthPosn$chr,
													truthPosn[,c("APOLLOHcall","ClonalCluster","CellularFrequency")])
ind <- match(truthDataIR,truthPosnIR) #indices for truthPosnIR
truthPosnDF <- as.data.frame(truthPosnIR)
truthDataDF <- as.data.frame(truthDataIR)

## assign the truth calls, clonal clust, cellular prev to truthData variable ##
truthDataDF <- cbind(truthDataDF,APOLLOHcall="HET",ClonalCluster=NA,CellularFrequency=NA)
truthDataDF[which(!is.na(ind)),c("APOLLOHcall","ClonalCluster","CellularFrequency")] <- truthPosnDF[ind[!is.na(ind)],c("APOLLOHcall","ClonalCluster","CellularFrequency")]

colnames(truthDataDF)[1:2] <- c("Chr","Position")

## load in original APOLLOH1 results ##
cna2 <- read.delim(inFile2,header=T,sep="\t")
chrsHMMcopy <- c(1,2,8,9,16,18)
cna2 <- cna2[cna2[,2] %in% chrsHMMcopy,]
loh2 <- read.delim(inFile3,header=F,sep="\t")
loh2 <- loh2[loh2[,1] %in% chrsHMMcopy,]


## load normal content of spikein sample ##
paramFile <- gsub("loh.txt","params.txt",inFile)
params <- read.delim(paramFile,header=F,sep="\t")
norm <- as.numeric(params[params[,1]=="Normal contamination estimate:",2])

library(SNPchip)
outFile <- paste(outDir,"/Spikein_plot_truth.png",sep="")
png(outFile,width=2400,height=1200,res=100)
#pdf(outFile,width=14,height=8,useDingbats=FALSE)
par(mfrow=c(2,1))
#layout(matrix(c(1:8),4,2,byrow=F))
source('/share/data/gha/software/LOH-HMM_0.1.0/plotting/plotSNPallelicRatio.R')
source('/share/data/gha/software/HMMcopy/scripts/plotCNlogRByChr.R')
plotCNlogRByChr(cna2, chr=NULL, geneAnnot=regions, spacing=4, yrange=c(-2,2), cytoBand=F, alphaVal=1)
plotSNPallelicRatio(loh2, chr=NULL, geneAnnot=regions, spacing=4, cytoBand=F, alphaVal=1)
dev.off()

for (i in chrs){
	outFile <- paste(outDir,"/Spikein_plot_chr",i,".png",sep="")
	png(outFile,width=1200,height=800,res=100)
	par(mfrow=c(3,1))
	source('/share/data/gha/software/APOLLOH_2.0.1/scripts/spikeInExpt/plotting.R')
	plotCNlogRByChr(loh, segments=NULL, chr=i, geneAnnot=regions, ylim=c(-2,2), spacing=4, alphaVal=1, xlim=NULL, cex=0.5)
	plotClonalFrequency(loh, chr=i, normal=norm, geneAnnot=NULL, spacing=4, ylim=c(0,1),xlim=NULL, cex=0.5)
	plotClonalFrequency(truthDataDF, chr=i, normal=0.45, geneAnnot=NULL, spacing=4, ylim=c(0,1),xlim=NULL, cex=0.5)
	pI <- plotIdiogram(i,build="hg19",unit="bp",label.y=-0.35,new=FALSE,ylim=c(-0.2,-0.1),label.cytoband=FALSE)
	dev.off()
}


