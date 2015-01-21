args <- commandArgs(TRUE)

options(stringsAsFactors=FALSE)

inFile <- args[1] #list of files for APOLLOH LOH
geneList <- args[2] #/share/data/gha/code2/scripts/GRCh37.p3_genes/GRCh37.p3_protein-coding.txt
method <- args[3] #{'common','severity','complete'}; default is 'severity'
headerType <- args[4] #{'Genes','chrPosn'}
filterLen <- as.numeric(args[5])
outDir <- args[6]



genes <- read.delim(geneList,header=T,sep="\t")
#colnames(genes)[1:6] <- c("EnsgId","Gene","Chr","Start","Stop","Strand")#,"Band","TxnCount","Description","Status")
colnames(genes)[1:3] <- c("Chr","Start","Stop")
genes[genes[,"Chr"]==23,"Chr"] <- "X"
genes[genes[,"Chr"]==24,"Chr"] <- "Y"
genes[genes[,"Chr"]==25,"Chr"] <- "MT"
chrPosn <- paste(genes$Chr,":",genes$Start,"-",genes$Stop,sep="")
genes <- cbind(genes, chrPosn)

stateKey <- array(0:9); names(stateKey) <- c('HOMD','DLOH','HET','NLOH','GAIN','ALOH','BCNA','UBCNA','ASCNA','OUT')
severity <- array(c(6,5,4,5,1,2,3,3,4,0)); names(severity) <- c('HOMD','DLOH','NLOH','ALOH','HET','GAIN','BCNA','UBCNA','ASCNA','OUT')


#########################################################################################
####################### FUNCTION: FIND OVERLAP OF SEGMENT AND GENE ######################
#########################################################################################
#Input: gene is an array with chr,position,start,stop
#Uses global variable "loh"
#Uses global variable "stateKey"
#Returns overlapping values for logR, AR, cell freq, clonal clust.
getLOH <- function(gene,loh){
	#print(as.character(gene["Gene"]))
	ind <- which(as.character(loh[,"Chr"])==as.character(gene["Chr"]) & 
		((as.numeric(loh[,"Start"])>=as.numeric(gene["Start"]) & as.numeric(loh[,"Start"])<=as.numeric(gene["Stop"])) |
		(as.numeric(loh[,"Stop"])>=as.numeric(gene["Start"]) & as.numeric(loh[,"Stop"])<=as.numeric(gene["Stop"])) |
		(as.numeric(loh[,"Start"])<=as.numeric(gene["Start"]) & as.numeric(loh[,"Stop"])>=as.numeric(gene["Stop"]))))
	

	geneCN <- NaN
	geneAR <- NaN
	geneLogR <- NaN
	geneCF <- NaN
	geneCC <- NaN
	stateInt <- NaN
	state <- NaN

	#if we found an overlap
	if (length(ind)>0){
	
		########## FIND STATE OF SEGMENT BASED ON METHOD ##################
		calls <- table(loh[ind,"TITAN_call"])
		lengths <- getOverlapLength(loh[ind,], gene["Start"], gene["Stop"])
		if (method=="common"){ 
			#use call of most frequent call of overlapping positions; in case of tie, it takes the first call based on the order from the table command
			#state <- names(which.max(calls))
			# for segments, will use largest segment
			state <- loh[ind[which.max(lengths)], "TITAN_call"]
		}else if(method=="complete"){
			#only accept call if all overlapping positions is same state; otherwise, leave as NaN
			if(length(calls)==1){  
				state <- names(calls)
			}
		}else{
			#use call of most severe state
			state <- getMostSevereState(names(calls))				
		}
		stateInt <- as.numeric(stateKey[state])
		
		########## USE DETERMINED STATE TO EXTRACT VALUES ##################
		if (state!="OUT"){
		   segInd <- which(loh[ind,"TITAN_call"]==state)
		   #stateInt <- loh[ind[segInd[1]],"APOLLOH_state"]
		   geneCN <- loh[ind[segInd[1]],"Copy_Number"]
		   geneAR <- loh[ind,"Median_Ratio"] %*% lengths / sum(lengths)
		   geneLogR <- loh[ind,"Median_logR"] %*% lengths / sum(lengths)
		   #find most common clonal freq/cluster for the segments with determined state
		   if (state!="HET"){
			   geneCF <- names(which.max(table(loh[ind[segInd],"Clonal_Frequency"])))
			   geneCC <- names(which.max(table(loh[ind[segInd],"Clonal_Cluster"])))
		   }
		}
	}
	outmat <- as.data.frame(cbind(geneCN,geneAR,geneLogR,geneCF,geneCC,stateInt))
	colnames(outmat) <- c("geneCN","geneAR","geneLogR","geneCF","geneCC","stateInt")
	return(outmat)
}	

## input: lohHits = loh rows that overlap region of interest
# start = start coordinate of region of interest
# end = end coordinate of region of interest
getOverlapLength <- function(lohHits, start, end){
	coords <- cbind(lohHits[, c("Start","Stop")], as.numeric(start), as.numeric(end))
	coordsSort <- t(apply(coords, 1, sort))
	dist <- coordsSort[, 3] - coordsSort[, 2] + 1
	return(dist)
}

#Input: states is an array of state names (e.g. (DLOH,NLOH,...,))
#Uses global variable "severity"
getMostSevereState <- function(states){	
	severityValue <- 0
	severeState <- states[1]
	for (i in states){
		if (severity[i] > severityValue){
			severeState <- i
			severityValue <- severity[i]
		}
	}
	return(severeState)
}


# output the matrix to file
#Input: Matrix to output; output file name
writeMatrixToFile <- function(mat,outfile){
	outMat <- cbind(rownames(mat),mat)
	if (!is.null(colnames(outMat))){
		colnames(outMat)[1] <- "Sample"
	}
	write.table(outMat,file=outfile,row.names=F,col.names=T,quote=F,na="NaN",sep="\t")
}


files <- read.delim(inFile,header=F,sep="\t")
numFiles <- nrow(files)

numSamples <- numFiles
numGenes <- dim(genes)[1]
geneCallmat <- matrix(NaN,nrow=numSamples,ncol=numGenes); colnames(geneCallmat) <- genes[,headerType]; 
geneCNmat <- matrix(NaN,nrow=numSamples,ncol=numGenes); colnames(geneCNmat) <- genes[,headerType]; 
geneARmat <- matrix(NaN,nrow=numSamples,ncol=numGenes); colnames(geneARmat) <- genes[,headerType]; 
geneLogRmat <- matrix(NaN,nrow=numSamples,ncol=numGenes); colnames(geneLogRmat) <- genes[,headerType]; 
geneCFmat <- matrix(NaN,nrow=numSamples,ncol=numGenes); colnames(geneCFmat) <- genes[,headerType]; 
geneCCmat <- matrix(NaN,nrow=numSamples,ncol=numGenes); colnames(geneCCmat) <- genes[,headerType]; 
rownames(geneCallmat) <- files[,1]
rownames(geneCNmat) <- files[,1]
rownames(geneARmat) <- files[,1]
rownames(geneLogRmat) <- files[,1]
rownames(geneCFmat) <- files[,1]
rownames(geneCCmat) <- files[,1]

for (i in 1:numFiles){
	caseId <- files[i,1]
	cat("Analyzing sample:\t",caseId," for file:",files[i,2],"...\n")

	#LOH
	loh <- read.delim(files[i,2],header=T,sep="\t")
	colnames(loh)[c(2,3,4)] <- c("Chr","Start","Stop")
	
	## filter by length threshold ##
	loh <- loh[loh[,"Length.bp."]>=filterLen,]

	#convert to gene scaffold
	overlaps <- apply(genes,1,getLOH,loh)
	geneTITAN <- do.call(rbind,overlaps)

	#build matrices
	geneCallmat[caseId,] <- geneTITAN[,"stateInt"]  #call matrix
	geneCNmat[caseId,] <- geneTITAN[,"geneCN"]  #call matrix
	geneARmat[caseId,] <- geneTITAN[,"geneAR"]   #AR matrix
	geneLogRmat[caseId,] <- geneTITAN[,"geneLogR"]   #logR matrix
	geneCFmat[caseId,] <- geneTITAN[,"geneCF"]   #cellular prevalence matrix
	geneCCmat[caseId,] <- geneTITAN[,"geneCC"]   #clonal cluster matrix

}

#geneCallmat[is.na(geneCallmat)] <- NaN
#geneARmat[is.na(geneARmat)] <- NaN
#geneLogRmat[is.na(geneLogRmat)] <- NaN
#geneCFmat[is.na(geneCFmat)] <- NaN
#geneCCmat[is.na(geneCCmat)] <- NaN

outCall <- paste(outDir,"_geneCalls.txt",sep="")
writeMatrixToFile(geneCallmat,outCall)
outCN <- paste(outDir,"_geneCN.txt",sep="")
writeMatrixToFile(geneCNmat,outCN)
outAR <- paste(outDir,"_geneAR.txt",sep="")
writeMatrixToFile(geneARmat,outAR)
outLogR <- paste(outDir,"_geneLogR.txt",sep="")
writeMatrixToFile(geneLogRmat,outLogR)
outCF <- paste(outDir,"_geneCF.txt",sep="")
writeMatrixToFile(geneCFmat,outCF)
outCC <- paste(outDir,"_geneCC.txt",sep="")
writeMatrixToFile(geneCCmat,outCC)


