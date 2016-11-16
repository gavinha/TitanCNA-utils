library(TitanCNA)
library(optparse)
sessionInfo()
args <- commandArgs(TRUE)
options(bitmapType='cairo', scipen=0)

option_list <- list(
	make_option(c("-l", "--libdir"), type = "character", help = "Directory containing code and reference wig files"),
	make_option(c("-id", "--id"), type = "character", help = "Sample ID"),
	make_option(c("-t", "--hetFile"), type = "character", help = "File containing allelic read counts at HET sites."),
	make_option(c("-c", "--cnFile"), type = "character", help = "File containing normalized log2 ratios"),
	make_option(c("-numC", "--numClusters"), type = "integer", default = 1, help = "Number of clonal clusters"),
	make_option(c("-cores", "--numCores"), type = "integer", default = 1, help = "Number of cores to use"),
	make_option(c("-p0", "--ploidy_0"), type = "numeric", default = 2, help = "Initial ploidy value"),
	make_option(c("-estP", "--estimatePloidy"), type = "logical", default = TRUE, help = "Estimate ploidy (TRUE)"),
	make_option(c("-n0", "--normal_0"), type = "numeric", default = 0.5, help = "Initial normal contamination (1-purity)"),
	make_option(c("-estN", "--estimateNormal"), type = "character", default = "map", help = "Estimate normal contamination method {'map', 'fixed'}"),
	make_option(c("-CN", "--maxCN"), type = "integer", default = 8, help = "Maximum number of copies to model"),
	make_option(c("-aK", "--alphaK"), type = "numeric", default = 1000, help = "Hyperparameter on Gaussian variance; defaults 1000 (WES) and 10000 (WGS)"),
	make_option(c("-aKH", "--alphaKHigh"), type = "numeric", default = 1000, help = "Hyperparameter on Gaussian variance for extreme copy number states"),
	make_option(c("-txnL", "--txnExpLen"), type = "numeric", default = 1e15, help = "Expected length of segments; higher leads to longer (less sensitive) segments"),
	make_option(c("-txnZ", "--txnZStrength"), type = "numeric", default = 1, help = "Expected length of clonal cluster segmentation (factor of txnExpLen)"),
	make_option(c("-minD", "--minDepth"), type = "numeric", default = 10, help = "Minimum read depth of a HET site to include in analysis"),
	make_option(c("-maxD", "--maxDepth"), type = "numeric", default = 1000, help = "Maximum read depth of a HET site to include in analysis"),
	make_option(c("--skew"), type = "numeric", default=0, help = "Allelic reference skew for all states except heterozygous states (e.g. 1:1, 2:2, 3:3). Value is additive to baseline allelic ratios."),
	make_option(c("--hetBaselineSkew"), type="numeric", default=NULL, help="Allelic reference skew for heterozygous states (e.g. 1:1, 2:2, 3:3). Value is the additive to baseline allelic ratios."), 
	make_option(c("-chrs", "--chrs"), type = "character", default = "c(1:22, 'X')", help = "Chromosomes to analyze"),
	make_option(c("-map", "--mapWig"), type = "character", default = NULL, help = "Mappability score file for bin sizes matching cnfile"),
	make_option(c("-mapThres", "--mapThres"), type = "character", default = 0.9, help = "Minimum mappability score threshold to use"),
	make_option(c("-centromere", "--centromere"), type = "character", help = "Centromere gap file"),
	make_option(c("-of", "--outFile"), type = "character", help = "Output file to write position-level file"),
	make_option(c("-os", "--outSeg"), type = "character", help = "Output file to write detailed segments"),
	make_option(c("-oi", "--outIGV"), type = "character", help = "Output file to write segments for loading into IGV"),
	make_option(c("-op", "--outParam"), type = "character", help = "Output file to write parameters"),
	make_option(c("--outPlotDir"), type = "character", help = "Output directory to save plots.")
)

parseobj <- OptionParser(option_list=option_list)
opt <- parse_args(parseobj)
print(opt)

libdir <- opt$libdir
source(paste0(libdir, "/v1.10.1/R/plotting.R"))
source(paste0(libdir, "/v1.10.1/R/utils.R"))
createSeg <- paste0(libdir, "/v1.10.1/createTITANsegmentfiles.pl")

id <- opt$id
hetfile <- opt$hetFile
cnfile <- opt$cnFile
numClusters <- opt$numClusters
numCores <- opt$numCores
ploidy_0 <- opt$ploidy_0
boolEstPloidy <- opt$estimatePloidy
norm_0 <- opt$normal_0
normEstMeth <- opt$estimateNormal
maxCN <- opt$maxCN
alphaK <- opt$alphaK
alphaHigh <- opt$alphaKHigh
txn_exp_len <- opt$txnExpLen
txn_z_strength <- opt$txnZStrength
mapThres <- opt$mapThres
minDepth <- opt$minDepth
maxDepth <- opt$maxDepth
skew <- opt$skew
hetBaselineSkew <- opt$hetBaselineSkew
chrs <- eval(parse(text = opt$chrs))
mapWig <- opt$mapWig
centromere <- opt$centromere 
outfile <- opt$outFile
outparam <- opt$outParam
outseg <- opt$outSeg
outigv <- opt$outIGV
outplot <- opt$outPlotDir

pseudo_counts <- 1e-300
centromereFlank <- 100000
maxI <- 50

message('Running TITAN...')

#### LOAD DATA ####
data <- loadAlleleCounts(hetfile, header=T)
 
#### REMOVE CENTROMERES ####
centromere <- read.delim(centromere,header=T,stringsAsFactors=F,sep="\t")

#### GC AND MAPPABILITY CORRECTION ####
message('titan: Loading GC content and mappability corrected log2 ratios...')
cnData <- read.delim(cnfile, header=T, stringsAsFactors=F, sep="\t")

#### READ COPY NUMBER FROM HMMCOPY FILE ####
message('titan: Extracting read depth...')
logR <- getPositionOverlap(data$chr,data$posn,cnData)
data$logR <- log(2^logR)
rm(logR,cnData)

#### FILTER DATA FOR DEPTH, MAPPABILITY, NA, etc ####
if (!is.null(mapWig)){
	mScore <- as.data.frame(wigToRangedData(mapWig))
	mScore <- getPositionOverlap(data$chr,data$posn,mScore[,-4])
	data <- filterData(data,chrs,minDepth=minDepth,maxDepth=maxDepth,map=mScore,mapThres=mapThres, centromeres = centromere)
	rm(mScore)
}else{
	data <- filterData(data,chrs,minDepth=minDepth,maxDepth=maxDepth,centromeres = centromere)
}

#### LOAD PARAMETERS ####
message('titan: Loading default parameters')
params <- loadDefaultParameters(copyNumber=maxCN,numberClonalClusters=numClusters, 
																skew=skew, hetBaselineSkew=hetBaselineSkew, data=data)

#### MODEL SELECTION USING EM (FWD-BACK) TO SELECT NUMBER OF CLUSTERS ####
library(doMC)
registerDoMC()
options(cores=numCores)
message("Using ",getDoParWorkers()," cores.")
K <- length(params$genotypeParams$rt)
params$genotypeParams$alphaKHyper <- rep(alphaK,K)
#params$genotypeParams$alphaKHyper[c(1,7:K)] <- alphaHigh 
params$ploidyParams$phi_0 <- ploidy_0
params$normalParams$n_0 <- norm_0
#params$genotypeParams$rt[c(4, 9)] <- hetAR

convergeParams <- runEMclonalCN(data,gParams=params$genotypeParams,nParams=params$normalParams,
                                pParams=params$ploidyParams,sParams=params$cellPrevParams,
                                maxiter=maxI,maxiterUpdate=1500,
                                txnExpLen=txn_exp_len,txnZstrength=txn_z_strength,
                                useOutlierState=FALSE,
                                normalEstimateMethod=normEstMeth,estimateS=TRUE,
                                estimatePloidy=boolEstPloidy, pseudoCounts=pseudo_counts)
    
#### COMPUTE OPTIMAL STATE PATH USING VITERBI ####
options(cores=1)
message("Using ",getDoParWorkers(),"cores.")
optimalPath <- viterbiClonalCN(data,convergeParams)

#### PRINT RESULTS TO FILES ####
results <- outputTitanResults(data,convergeParams,optimalPath,
			filename=NULL,posteriorProbs=F,subcloneProfiles=TRUE)
corrResults <- removeEmptyClusters(convergeParams, results, proportionThreshold = 0.02, proportionThresholdClonal = 0.05)
convergeParams <- corrResults$convergeParams
results <- corrResults$results
numClustersToPlot <- nrow(convergeParams$s)
write.table(results, file = outfile, col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
outputModelParameters(convergeParams, results, outparam)


if (numClusters < 10) { numClusters <- paste("0",numClusters,sep="") }
# save specific objects to a file
# if you don't specify the path, the cwd is assumed 
convergeParams$rhoG <- NULL; convergeParams$rhoZ <- NULL
outImage <- gsub(".titan.txt", ".RData", outfile)
save(convergeParams, file=outImage)

#### OUTPUT SEGMENTS ####
message("Writing segments to ", outseg)
segCmd <- paste0(createSeg, " -id ", id, " -i ", outfile, " -o ",  outseg, " -igv ", outigv)
system(segCmd)
segs <- read.delim(outseg, header = TRUE, stringsAsFactors = FALSE, sep = "\t")

#### PLOT RESULTS ####
dir.create(outplot)
library(SNPchip)  ## use this library to plot chromosome idiogram (optional)
norm <- tail(convergeParams$n,1)
ploidy <- tail(convergeParams$phi,1)
for (chr in chrs){
	outfig <- paste(outplot,"/",id,"_cluster",numClusters,"_chr",chr,".png",sep="")
	png(outfig,width=1200,height=1000,res=100)
	if (as.numeric(numClusters) <= 2){
		par(mfrow=c(4,1))
	}else{
		par(mfrow=c(3,1))  
	}
	plotCNlogRByChr(results, chr, segs = segs, ploidy=ploidy, normal = norm, geneAnnot=NULL,  cex.axis=1.5, 
					ylim=c(-2,2), cex=0.5, xlab="", main=paste("Chr ",chr,sep=""))
	plotAllelicRatio(results, chr, geneAnnot=NULL, spacing=4, cex.axis=1.5,
					ylim=c(0,1), xlab="", cex=0.5, main=paste("Chr ",chr,sep=""))
	plotClonalFrequency(results, chr, normal=norm, geneAnnot=NULL, spacing=4, 
					cex.axis=1.5, ylim=c(0,1), xlab="", cex=0.5, main=paste("Chr ",chr,sep=""))
                    
	if (as.numeric(numClustersToPlot) <= 2 && as.numeric(numClusters) <= 2){
		plotSubcloneProfiles(results, chr, cex = 2, spacing=6, main=paste("Chr ",chr,sep=""), cex.axis=1.5)
		pI <- plotIdiogram(chr, build="hg19", unit="bp", label.y=-4.25, new=FALSE, ylim=c(-2,-1))
	}else{
		pI <- plotIdiogram(chr, build="hg19", unit="bp", label.y=-0.35, new=FALSE, ylim=c(-0.2,-0.1))
	}
	
	dev.off()
}

################################################
############## GENOME WIDE PLOTS ###############
################################################
outFile <- paste(outplot,"/",id,"_cluster",numClusters,"_CNA.pdf",sep="")
#png(outFile,width=1000,height=300)
pdf(outFile,width=20,height=6)
plotCNlogRByChr(dataIn=results, chr=NULL, segs = segs, ploidy=ploidy,  normal = norm, geneAnnot=genes, spacing=4, main=id, xlab="", ylim=c(-2,2), cex=0.5, cex.axis=1.5, cex.lab=1.5, cex.main=1.5)
dev.off()

outFile <- paste(outplot,"/",id,"_cluster",numClusters,"_LOH.pdf",sep="")
#png(outFile,width=1000,height=300)
pdf(outFile,width=20,height=6)
plotAllelicRatio(dataIn=results, chr=NULL, geneAnnot=genes, spacing=4, main=id, xlab="", ylim=c(0,1), cex=0.5, cex.axis=1.5, cex.lab=1.5, cex.main=1.5)	
dev.off()

outFile <- paste(outplot,"/",id,"_cluster",numClusters,"_CF.pdf",sep="")
#png(outFile,width=1000,height=300)
pdf(outFile,width=20,height=6)
plotClonalFrequency(dataIn=results, chr=NULL, norm, geneAnnot=genes, spacing=4, main=id, xlab="", ylim=c(0,1), cex.axis=1.5, cex.lab=1.5, cex.main=1.5)
dev.off()

if (as.numeric(numClusters) <= 2){
	outFile <- paste(outplot,"/",id,"_cluster",numClusters,"_subclone.pdf",sep="")
	#png(outFile,width=1000,height=300)
	pdf(outFile,width=20,height=6)
	plotSubcloneProfiles(dataIn=results, chr=NULL, cex = 0.5, spacing=4, main=id, cex.axis=1.5, xlab="")
	dev.off()
}

