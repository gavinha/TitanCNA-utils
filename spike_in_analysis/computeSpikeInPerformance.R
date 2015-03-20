args <- commandArgs(TRUE)

options(stringsAsFactors=FALSE)


inFile <- args[1]
truthSegFile <- args[2]
truthPosnFile <- gsub("segs.txt","posns.txt",truthSegFile)
truthCellPrevFile <- gsub("segs.txt","cellPrev.txt",truthSegFile)
outDir <- args[3]
cp.buffer <- 0.10
id <- gsub("_loh.txt","",basename(inFile))
#############################################################################################
################################# READ FILES AND PREPARE OBJECTS#############################
#############################################################################################
library(IRanges)
#### READ IN SEG AND POSN GROUND TRUTH FILES ####
message("Loading ground truth data...")
#id	chr	start	end	length	numPosns	AllelicRatio	LogR	space	start	end	type
truthSeg <- read.delim(truthSegFile,header=T,sep="\t")
#chr	start	id	space	start	refCount	nrefCount	logR	type.space	type.start	type.end	type.type
truthPosn <- read.delim(truthPosnFile,header=T,sep="\t")
truthCellPrev <- read.delim(truthCellPrevFile,row.names=1,header=F,sep="\t")

truthSegIR <- RangedData(ranges=IRanges(start=as.numeric(truthSeg$start),end=as.numeric(truthSeg$end)),
													space=as.numeric(truthSeg$chr),numPosns=as.numeric(truthSeg$numPosns),
													id=truthSeg$id,eventType=truthSeg$type)
truthSegIR <- truthSegIR[order(space(truthSegIR),start(truthSegIR)),]
truthPosnIR <- RangedData(ranges=IRanges(start=as.numeric(truthPosn$start),end=as.numeric(truthPosn$start)),
													space=as.numeric(truthPosn$chr),eventType=truthPosn$type.type)
truthPosnIR <- truthPosnIR[order(space(truthPosnIR),start(truthPosnIR)),]
truthPosn <- as.data.frame(truthPosnIR)
#### READ IN TITAN PREDICTION POSITION-BASED FIL ####
cat("Loading ",inFile, "...\n")
loh <- read.delim(inFile,header=T,sep="\t")
loh <- loh[,1:12]
numClust <- max(loh[,"ClonalCluster"],na.rm=T)

#read in param file
paramFile <- gsub("loh.txt","params.txt",inFile)
params <- read.delim(paramFile,row.names=1,header=F,sep="\t")
norm <- as.numeric(params["Normal contamination estimate:",1])
#recompute cellular frequencies with normal (% of cells WITH event)
loh[,"CellularFrequency"] <- (1-norm)*(1-loh[,"CellularFrequency"])

#############################################################################################
########################## FUNCTION TO COMPUTE PERFORMANCE  ################################
############################################################################################
## assumes global variables: truthCellPrev
## x is the merge(loh predictions, truth positions) with columns:
##   space   start     end width CopyNumber APOLLOHcall CellularFrequency eventType
## info is a row from truthSegLen with columns:
##  space   start     end width numPosn  id  eventType 
## cp.buffer= cellular prevalence +/- buffer for true positive match
computePerformance <- function(x,info,cp.buffer=0.05){
	## count number of true positives based on whether it is DEL or AMP
	if (info$eventType=="DEL"){
		indEvent <- x$CopyNumber < 2 & x$eventType=="DEL"
	}else if (info$eventType=="AMP"){
		indEvent <- x$CopyNumber > 2 & x$eventType=="AMP"
	}
	
	## count number of true positives for cellular prevalence
	indCP <- x$CellularFrequency>(truthCellPrev[info$id,1]-cp.buffer) & x$CellularFrequency<(truthCellPrev[info$id,1]+cp.buffer)
	
	outMat <- cbind(numPosnUsed=nrow(x),eventTP=sum(indEvent,na.rm=T),cpTP=sum(indCP,na.rm=T))
	return(as.data.frame(outMat))
}

#############################################################################################
#################### FUNCTION TO COMPUTE SUMMARY PERFORMANCE  ###############################
############################################################################################
## summarizes performance for events of the same length
getOverallPerformance <- function(perf,type=NULL,TPthres=0.9){
	if (!is.null(type)){
		ind <- perf$eventType==type
	}else{
		ind <- !logical(length=nrow(perf))
	}
	numPosns <- sum(perf[ind,"numPosns"],na.rm=T)
	medianLen <- median(perf[ind,"width"])
	eventTP <- sum(perf[ind,"eventTP"],na.rm=T)
	eventNumPosns <- sum(perf[ind,"numPosnUsed"],na.rm=T)
	eventTPR <- eventTP / eventNumPosns
	cpTP <- sum(perf[ind,"cpTP"],na.rm=T)
	cpNumPosns <- sum(perf[ind,"numPosnUsed"],na.rm=T)
	cpTPR <- cpTP / cpNumPosns
	numTPevents <- sum(perf[ind,"eventTPR"]>TPthres)
	numTPcp <- sum(perf[ind,"cpTPR"]>TPthres)
	outMat <- cbind(numPosns=numPosns,medianLength=medianLen,eventNumPosns,eventTP,eventTPR,numTPevents,cpNumPosns,cpTP,cpTPR,numTPcp)
	return(as.data.frame(outMat))
}

#############################################################################################
#################### FUNCTION TO COMPUTE SUMMARY PERFORMANCE  ###############################
############################################################################################
## summarizes performance for all events for all sizes
getAllEventSummary <- function(perf){
	allEventNumPosns <- sum(perf$eventNumPosns)
	allMedianLength <- NA
	allEventTP <- sum(perf$eventTP)
	allCPNumPosns <- sum(perf$cpNumPosns)
	allCPTP <- sum(perf$cpTP)
	totalNumPosns <- sum(perf$numPosns)
	totalNumTPevents <- sum(perf$numTPevents)
	totalNumTPcp <- sum(perf$numTPcp)
	outMat <- cbind(numSNPs="total",numPosns=totalNumPosns,medianLength=allMedianLength,eventNumPosns=allEventNumPosns,
												eventTP=allEventTP,eventTPR=allEventTP/allEventNumPosns,numTPevents=totalNumTPevents,
												cpNumPosns=allCPNumPosns,cpTP=allCPTP,cpTPR=allCPTP/allCPNumPosns,numTPcp=totalNumTPcp)
	return(as.data.frame(outMat))												
}


#############################################################################################
########################## COMPUTE PERFORMANCE FOR EACH EVENT SIZE ##########################
#############################################################################################
eventLengths <- unique(truthSeg$numPosns)
chrList <- sort(as.numeric(unique(truthSeg$chr)))
#### Reduce LOH to only chromosomes of interest ####
loh <- loh[loh$Chr %in% chrList,]  
lohIR <- RangedData(space=as.numeric(loh$Chr),ranges=IRanges(start=as.numeric(loh$Position),end=as.numeric(loh$Position)),
										loh[,c("CopyNumber","APOLLOHcall","CellularFrequency")])
lohIR <- lohIR[order(space(lohIR),start(lohIR)),]
## reassign LOH so that the rows match the lohIR object
loh <- as.data.frame(lohIR)

#### FIND TRUE NEGATIVE POSITIONS ####
ind <- unlist(as.list(overlapsAny(lohIR,truthPosnIR,type="within")))
truthNeg <- as.data.frame(lohIR[!ind,])
FP <- sum(truthNeg$CopyNumber!=2,na.rm=T)
numNeg <- sum(!is.na(truthNeg$CopyNumber))
fpr <- FP/numNeg  #FALSE POSITIVE RATE

perfSegResult <- NULL
perfEventResult <- NULL
perfAMPResult <- NULL
perfDELResult <- NULL
#### COMPUTE PERFORMANCE FOR EACH EVENT SIZE ####
for (i in 1:length(eventLengths)){
	message("Analyzing length ",eventLengths[i])
	## find truth positions corresponding to event length (from truthSegs)
	truthSegLenIR <- truthSegIR[truthSegIR$numPosns==eventLengths[i],]
	truthSegLen <- as.data.frame(truthSegLenIR)
	## find predicted positions corresponding to event segment
	perfSegResultByEvent <- NULL
	for (j in 1:nrow(truthSegLenIR)){		
		message("Computing performance for segment ",truthSegLen[j,])
		truthPosnLenIR <- subsetByOverlaps(truthPosnIR,truthSegLenIR[j,],type="within")
		lohLenIR <- subsetByOverlaps(lohIR,truthPosnLenIR,type="within")
		lohTruthJoin <- merge(as.data.frame(lohLenIR),as.data.frame(truthPosnLenIR),sort=TRUE)
		
		perf <- computePerformance(lohTruthJoin,truthSegLen[j,],cp.buffer=cp.buffer)
		perfSegResultByEvent <- rbind(perfSegResultByEvent,cbind(truthSegLen[j,],perf,eventTPR=perf$eventTP/perf$numPosn,cpTPR=perf$cpTP/perf$numPosn))
	}
	perfSegResult <- rbind(perfSegResult,perfSegResultByEvent) ##add to global list of segs
	
	##compute performance for all segments of the same length
	perfAMPResult <- rbind(perfAMPResult,cbind(numSNPs=eventLengths[i],
											getOverallPerformance(perfSegResultByEvent,type="AMP")))
	perfDELResult <- rbind(perfDELResult,cbind(numSNPs=eventLengths[i],
											getOverallPerformance(perfSegResultByEvent,type="DEL")))
	perfEventResult <- rbind(perfEventResult,cbind(numSNPs=eventLengths[i],
											getOverallPerformance(perfSegResultByEvent,type=NULL)))
}

#### compute overall performance for all event sizes combined ####
perfAMPResult <- rbind(perfAMPResult,getAllEventSummary(perfAMPResult))
perfDELResult <- rbind(perfDELResult,getAllEventSummary(perfDELResult))
perfEventResult <- rbind(perfEventResult,getAllEventSummary(perfEventResult))

perfEventResult <- cbind(perfEventResult,DEL=perfDELResult,AMP=perfAMPResult,numNeg=numNeg,allFPR=fpr)

outFile <- paste(outDir,"/",id,"_eventPerformance.txt",sep="")
write.table(perfEventResult,file=outFile,col.names=T,row.names=F,quote=F,sep="\t")
outFile <- paste(outDir,"/",id,"_segmentPerformance.txt",sep="")
write.table(perfSegResult,file=outFile,col.names=T,row.names=F,quote=F,sep="\t")

#############################################################################################
########################## PLOT PERFORMANCE FOR ACROSS EVENT SIZE ###########################
#############################################################################################
## one plot for each of overall, loss, and gain
## one for each of numTPevents, eventTPR and numTPcp, cellPrev TPR
outFile <- paste(outDir,"/",id,"_perf_plot.pdf",sep="")
pdf(outFile,width=11,height=8,useDingbats=FALSE)

par(mfrow=c(3,4),xpd=NA) 
totalNumEvents <- c(24,24,24,6)

################## ALL EVENT ####################
## plot proportion of numTPevents ##
plot(as.numeric(perfEventResult$numTPevents[-5])/totalNumEvents, type="o", pch=19, cex=2, ylim=c(0,1), ylab="Number of Events", xlab="Event Length (bp)", xaxt="n", las=2, cex.axis=1.5, cex.lab=1.5, cex.main=1.5, main="Prediction for\nAll Spike-In Events")
axis(1,labels=c(10,100,1000,10000), at=1:4, cex.axis=1.5)
text(x=c(1.2,2,3,3.8),y=0.05,labels=paste("n=",totalNumEvents,sep=""),cex=1.5)

## plot event TPR ##
plot(as.numeric(perfEventResult$eventTPR[-5]), type="o", pch=19, cex=2, ylim=c(0,1), ylab="True Positive Rate (TPR)", xlab="Event Length (bp)", xaxt="n", las=2, cex.axis=1.5, cex.lab=1.5, cex.main=1.5, main="Prediction for\nAll Spike-In Events")
axis(1,labels=c(10,100,1000,10000), at=1:4, cex.axis=1.5)

## plot proportion of cpTP ##
plot(as.numeric(perfEventResult$numTPcp[-5])/totalNumEvents, type="o", pch=19, cex=2, ylim=c(0,1), ylab="Number of Events", xlab="Event Length (bp)", xaxt="n", las=2, cex.axis=1.5, cex.lab=1.5, cex.main=1.5, main="Cellular Prevalence for\nAll Spike-In Events")
axis(1,labels=c(10,100,1000,10000), at=1:4, cex.axis=1.5)
text(x=c(1.2,2,3,3.8),y=0.05,labels=paste("n=",totalNumEvents,sep=""),cex=1.5)

## plot cp TPR ##
plot(as.numeric(perfEventResult$cpTPR[-5]), type="o", pch=19, cex=2, ylim=c(0,1), ylab="True Positive Rate (TPR)", xlab="Event Length (bp)", xaxt="n", las=2, cex.axis=1.5, cex.lab=1.5, cex.main=1.5, main="Cellular Prevalence for\nAll Spike-In Events")
axis(1,labels=c(10,100,1000,10000), at=1:4, cex.axis=1.5)

################## DELETIONS ####################
totalNumEvents <- c(24,24,24,6)/2
## plot proportion of numTPevents ##
plot(as.numeric(perfEventResult$DEL.numTPevents[-5])/totalNumEvents, type="o", pch=19, cex=2, ylim=c(0,1), ylab="Number of Events", xlab="Event Length (bp)", xaxt="n", las=2, cex.axis=1.5, cex.lab=1.5, cex.main=1.5, main="Prediction for\nDeletion Spike-In Events")
axis(1,labels=c(10,100,1000,10000), at=1:4, cex.axis=1.5)
text(x=c(1.2,2,3,3.8),y=0.05,labels=paste("n=",totalNumEvents,sep=""),cex=1.5)

## plot event TPR ##
plot(as.numeric(perfEventResult$DEL.eventTPR[-5]), type="o", pch=19, cex=2, ylim=c(0,1), ylab="True Positive Rate (TPR)", xlab="Event Length (bp)", xaxt="n", las=2, cex.axis=1.5, cex.lab=1.5, cex.main=1.5, main="Prediction for\nDeletion Spike-In Events")
axis(1,labels=c(10,100,1000,10000), at=1:4, cex.axis=1.5)

## plot proportion of cpTP ##
plot(as.numeric(perfEventResult$DEL.numTPcp[-5])/totalNumEvents, type="o", pch=19, cex=2, ylim=c(0,1), ylab="Number of Events", xlab="Event Length (bp)", xaxt="n", las=2, cex.axis=1.5, cex.lab=1.5, cex.main=1.5, main="Cellular Prevalence for\nDeletion Spike-In Events")
axis(1,labels=c(10,100,1000,10000), at=1:4, cex.axis=1.5)
text(x=c(1.2,2,3,3.8),y=0.05,labels=paste("n=",totalNumEvents,sep=""),cex=1.5)

## plot cp TPR ##
plot(as.numeric(perfEventResult$DEL.cpTPR[-5]), type="o", pch=19, cex=2, ylim=c(0,1), ylab="True Positive Rate (TPR)", xlab="Event Length (bp)", xaxt="n", las=2, cex.axis=1.5, cex.lab=1.5, cex.main=1.5, main="Cellular Prevalence for\nDeletion Spike-In Events")
axis(1,labels=c(10,100,1000,10000), at=1:4, cex.axis=1.5)

################## AMP ####################
totalNumEvents <- c(24,24,24,6)/2
## plot proportion of numTPevents ##
plot(as.numeric(perfEventResult$AMP.numTPevents[-5])/totalNumEvents/2, type="o", pch=19, cex=2, ylim=c(0,1), ylab="Number of Events", xlab="Event Length (bp)", xaxt="n", las=2, cex.axis=1.5, cex.lab=1.5, cex.main=1.5, main="Prediction for\nAmplification Spike-In Events")
axis(1,labels=c(10,100,1000,10000), at=1:4, cex.axis=1.5)
text(x=c(1.2,2,3,3.8),y=0.05,labels=paste("n=",totalNumEvents,sep=""),cex=1.5)

## plot event TPR ##
plot(as.numeric(perfEventResult$AMP.eventTPR[-5]), type="o", pch=19, cex=2, ylim=c(0,1), ylab="True Positive Rate (TPR)", xlab="Event Length (bp)", xaxt="n", las=2, cex.axis=1.5, cex.lab=1.5, cex.main=1.5, main="Prediction for\nAmplification Spike-In Events")
axis(1,labels=c(10,100,1000,10000), at=1:4, cex.axis=1.5)

## plot proportion of cpTP ##
plot(as.numeric(perfEventResult$AMP.numTPcp[-5])/totalNumEvents/2, type="o", pch=19, cex=2, ylim=c(0,1), ylab="Number of Events", xlab="Event Length (bp)", xaxt="n", las=2, cex.axis=1.5, cex.lab=1.5, cex.main=1.5, main="Cellular Prevalence for\nAmplification Spike-In Events")
axis(1,labels=c(10,100,1000,10000), at=1:4, cex.axis=1.5)
text(x=c(1.2,2,3,3.8),y=0.05,labels=paste("n=",totalNumEvents,sep=""),cex=1.5)

## plot cp TPR ##
plot(as.numeric(perfEventResult$AMP.cpTPR[-5]), type="o", pch=19, cex=2, ylim=c(0,1), ylab="True Positive Rate (TPR)", xlab="Event Length (bp)", xaxt="n", las=2, cex.axis=1.5, cex.lab=1.5, cex.main=1.5, main="Cellular Prevalence for\nAmplification Spike-In Events")
axis(1,labels=c(10,100,1000,10000), at=1:4, cex.axis=1.5)


dev.off()
