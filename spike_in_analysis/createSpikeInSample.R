args <- commandArgs(TRUE)

options(stringsAsFactors=FALSE)

inFile <- args[1]
regionsFile <- args[2]
numEvents <- as.numeric(args[3])
outDir <- args[4]

library(TitanCNA)


############################################################
############## LOAD FILES AND REGIONS ######################
############################################################
files <- read.delim(inFile,header=F,sep="\t")
numFiles <- nrow(files)
regions <- read.delim(regionsFile,header=T,sep="\t")
regions <- RangedData(ranges=IRanges(start=regions$Start,end=regions$Stop),type=regions$Type,space=regions$Chr)


########### setup new sample objects #######################
# assign the first files because they are the 100% tumour
# not RangedData object because can't change values in those objects
newAR <- read.delim(files[1,3],header=F,sep="\t")
newAR[newAR[,1]=="X",1] <- 23; newAR[newAR[,1]=="Y",1] <- 24; newAR <- newAR[newAR[,1]!="MT",]
colnames(newAR) <- c("chr","start","refBase","refCount","nrefBase","nrefCount")
newARIR <- RangedData(ranges=IRanges(start=newAR$start,end=newAR$start),
					space=as.numeric(newAR$chr),newAR[,3:6])

newLogR <- read.delim(files[1,2],header=T,sep="\t")
newLogR[newLogR[,1]=="X",1] <- 23; newLogR[newLogR[,1]=="Y",1] <- 24; newLogR <- newLogR[newLogR[,1]!="MT",]
## need RangedData object to help us find where to replace logR values
newLogRIR <- RangedData(ranges=IRanges(start=newLogR$start,end=newLogR$end),
						space=as.numeric(newLogR$chr),logR=newLogR$logR)
						
## find the indicies for which LogR overlaps AR
indAR.LogR <- as.matrix(findOverlaps(newARIR,newLogRIR))
indAR.LogR <- indAR.LogR[sort(indAR.LogR[,1]),]  #sort to make sure that ar indices are consecutive

## add logR info to AR RangedData object ##
newARIR$logR <- newLogRIR[indAR.LogR[,2],][["logR"]]

## convert newARIR to newAR (data.frame) to work with later as the main output 
newAR <- as.data.frame(newARIR)[,c(-3,-4)]
colnames(newAR)[1] <- "chr"

############################################################
###### FUNCTION TO SEARCH FOR REGIONS TO REPLACE ###########
############################################################
## repeatedly sample the starting point until found neutral region with original values
## and selected start point does not overlap used Regions
## global variables: "newLogRIR"
replaceRegionsWithSpikeIn <- function(s,chrs,sampleAR,ar,usedRegions,id){
	repeat{
		chr <- sample(chrs,size=1) #randomly choose a chromosome from chrs
		chrInd <- which(ar$chr==chr) 
		ind <- sort(sample(chrInd,size=1))
		ind <- ind:(ind+s-1)
		tmpIR <- RangedData(ranges=IRanges(start=ar$start[ind],end=ar$start[ind]),
							space=ar$chr[ind],refCount=ar$refCount[ind],
							nrefCount=ar$nrefCount[ind],logR=ar$logR[ind])
		if (sum(unlist(as.list(overlapsAny(tmpIR,usedRegions))))==0){  #no overlap of existing positions
			if (length(unique(ar$chr[ind]))==1){#does not extend past the chromosome.
				if (tail(ar$start[ind],n=1) <= max(ar$start[chrInd],na.rm=T)){ # this new region is s-number of positions and fits 
					break
				}
			}
		}
		message("Searching for replacement region of size=",s," in chr",chr)
	}
	## add replaced regions to used regions list
	usedRegions <- rbind(usedRegions,tmpIR)
	## replace AR values 
	newAR <- ar
	newAR[ind,c(4,6,7)] <- cbind(sampleAR$refCount,sampleAR$nrefCount,sampleAR$logR)
	## keep track of the new coordinates and the new segment
	newCoords <- cbind(newAR[ind,1:2],id,as.data.frame(sampleAR)[,c(-3,-4)])
	
	## setup single row output of the segment with median allelic ratio and logR
	tmpDF <- as.data.frame(tmpIR)
	newMedianAR <- median(apply(cbind(sampleAR$refCount,sampleAR$nrefCount),1,max)
						/(sampleAR$refCount+sampleAR$nrefCount),na.rm=T)
	newMedianLogR <- median(sampleAR$logR,na.rm=T)
	newSegs <- cbind(id,chr=tmpDF[1,1],start=tmpDF[1,2],end=tail(tmpDF,n=1)[,3],
					length=tail(tmpDF,n=1)[,3]-tmpDF[1,2]+1,numPosns=nrow(newCoords),
					AllelicRatio=newMedianAR,LogR=newMedianLogR)
	
	## output the objects to a list
	out <- list()
	out$usedRegions <- tmpIR
	out$newAR <- newAR
	out$newCoords <- newCoords
	out$newSegs <- newSegs
	return(out)
}

############################################################
########### SAMPLE POSITIONS FROM EACH MIXTURE #############
############################################################

sampleSize <- c(rep(10,numEvents),rep(100,numEvents),rep(1000,numEvents),rep(10000,1))
chrs <- c(1,2,9,18)
## keep track of regions used for replacement (in neutral regions of chr9,18,1,2)
usedRegions <- RangedData()
newCoords <- NULL
newSegs <- NULL

for (i in 1:numFiles){
	message("Using file ",files[i,1])
	########### load depth and allelic ratios ##############
	id <- files[i,1]
	ar <- read.delim(files[i,3],header=F,sep="\t")
	ar[ar[,1]=="X",1] <- 23; ar[ar[,1]=="Y",1] <- 24; ar <- ar[ar[,1]!="MT",]
	ar <- RangedData(ranges=IRanges(start=ar[,2],end=ar[,2]),
					space=as.numeric(ar[,1]),refCount=ar[,4],nrefCount=ar[,6])
	logR <- read.delim(files[i,2],header=T,sep="\t")
	logR[logR[,1]=="X",1] <- 23; logR[logR[,1]=="Y",1] <- 24; logR <- logR[logR[,1]!="MT",]
	logR <- RangedData(ranges=IRanges(start=logR$start,end=logR$end),
						space=as.numeric(logR$chr),logR=logR$logR)	
	## find the indicies for which LogR overlaps AR
	indAR.LogR <- as.matrix(findOverlaps(ar,logR))
	indAR.LogR <- indAR.LogR[sort(indAR.LogR[,1]),]  #sort to make sure that ar indices are consecutive
	## add logR info to AR RangedData object ##
	ar$logR <- logR[indAR.LogR[,2],][["logR"]]

	########### sample the positions from regions ##############
	for (r in 1:nrow(regions)){
		for (s in sampleSize){	
			###### sample the region positions ######
			ind <- which(unlist(as.list(overlapsAny(ar,regions[r,],type="within"))))
			ARind <- sort(sample(ind,size=s,replace=FALSE))
			sampleAR <- ar[ARind,]
			
			###### replace existing values ######
			message("Replacing data for segment size=",s,
							" from chr",as.character(regions[r,][["space"]]),
							" ",paste(as.data.frame(regions[r,])[,c(-1,-4)],collapse=" "))			
			#if (s==10){ ## for smaller variants (sized 10), place in chr18	
				results <- replaceRegionsWithSpikeIn(s,chrs,sampleAR,newAR,usedRegions,id)
			#}else if (s==100){ ## for medium variants (sized 100), place in chr9
			#	results <- replaceRegionsWithSpikeIn(s,chr="9",sampleAR,newAR,usedRegions)
			#}else if (s==1000){ ## for large variants (sized 1000), place in chr2
			#	results <- replaceRegionsWithSpikeIn(s,chr="2",sampleAR,newAR,usedRegions)
			#}else if (s==10000){ ## for x-large variants (sized 10000), place in chr1
			#	results <- replaceRegionsWithSpikeIn(s,chr="1",sampleAR,newAR,usedRegions)
			#}
			newAR <- results$newAR  #reassign because some positions have replaced values
			usedRegions <- rbind(usedRegions,results$usedRegions)
			newCoords <- rbind(newCoords,cbind(results$newCoords,type=as.data.frame(regions[r,])[,-4]))
			newSegs <- rbind(newSegs,cbind(results$newSegs,as.data.frame(regions[r,])[,-4]))
		}
	}
}


############################################################
################## OUTPUT THE RESULTS ######################
############################################################

############## output new TITAN input file #################
outFile1 <- paste(outDir,"/SpikeIn_TumCounts_LogR.txt",sep="")
write.table(newAR,file=outFile1,col.names=T,row.names=F,quote=F,sep="\t")

############ output truth coordinates positions file #######
outFile2 <- paste(outDir,"/SpikeIn_truth_posns.txt",sep="")
write.table(newCoords,file=outFile2,col.names=T,row.names=F,quote=F,sep="\t")

############## output truth segment file ###################
outFile3 <- paste(outDir,"/SpikeIn_truth_segs.txt",sep="")
write.table(newSegs,file=outFile3,col.names=T,row.names=F,quote=F,sep="\t")






