## R script ##
## Search the space of the S_Dbw scale factor (S_Dbw=scale*dens.bw+scat)
##  until we find the inflection (bend) in the curve (plateau)
## scale is selected based on diminishing effects on distribution of samples across number of clonal clusters
library(inflection)

args <- commandArgs(TRUE)
options(stringsAsFactors=FALSE)

inDir <- args[1]
maxScale <- as.numeric(args[2])
type <- "Both" #args[3]
outFile <- args[3]
outPlot <- args[4]

files <- list.files(inDir, pattern="params.txt", full.names=T)
idClusts <- as.data.frame(do.call(rbind,
			strsplit(gsub("_params.txt","",basename(files)),split="_cluster")))
ids <- unique(idClusts[,1])
clusts <- unique(as.numeric(idClusts[,2]))

matDens <- matrix(NA, nrow=length(ids), ncol=length(clusts), dimnames=list(ids, clusts))
matScat <- matDens

for (i in 1:length(files)){
	param <- read.delim(files[i], header=F, sep="\t")
	
	idClust <- strsplit(sub("_params.txt","",
					basename(files[i])),split="_cluster")[[1]]
	
	dens.bw <- as.numeric(param[grep(paste("dens.bw \\(",type,"\\)",sep=""), 
							param[,1]), 2])
	matDens[idClust[1],as.numeric(idClust[2])] <- dens.bw
	scat <-  as.numeric(param[grep(paste("scat \\(",type,"\\)",sep=""), 
							param[,1]), 2])
	matScat[idClust[1],as.numeric(idClust[2])] <- scat

}

## compute S_Dbw for each possible value of scale
## S_Dbw = scale * dens + scat
scale <- 1:maxScale
ratio <- rep(NA, length(scale))
for (i in 1:length(scale)){
	SDbw <- scale[i] * matDens + matScat
	counts <- table(apply(SDbw, 1, which.min))
	ratio[i] <- counts[1] / sum(counts[-1])  # clust1 / sum_C (clustC)
	#ratio[i] <- counts[2] / sum(counts[-c(1,2)])  #clust2 / sum(rest of clust, minus 1 & 2)
}

## plot ratio curve ##
pdf(outPlot)
plot(scale, ratio, type="o", pch=19, xaxt="n", las=2, ylab="Clonal:Subclonal Ratio", 
		xlab="Scale", main="S_Dbw = scale * dens + scat")
axis(1)

## find and plot inflection plot using EDE (R package inflection)
inflectPt <- findiplist(x=as.matrix(scale), y=ratio, index=1)
#abline(v=inflectPt[1,1], col="red"); mtext("left",side=3,at=inflectPt[1,1])
#abline(v=inflectPt[1,3], col="red"); mtext("ESE",side=3,at=inflectPt[1,3])
#abline(v=inflectPt[1,2], col="red"); mtext("right",side=3,at=inflectPt[1,2])
abline(v=inflectPt[2,1], col="red")
mtext(text=scale[inflectPt[2,1]],side=3,at=inflectPt[2,1])

dev.off()

## output if outFile is given (not "0")
if (outFile != "0"){
	SDbw <- scale[inflectPt[2,1]] * matDens + matScat
	optClust <- apply(SDbw, 1, which.min)
	outStr <- paste(ids,"_cluster0",optClust,sep="")
	write.table(outStr,file=outFile,row.names=F,col.names=F,quote=F,sep="\t")
}

message("Using scale ", scale[inflectPt[2,1]],".")
