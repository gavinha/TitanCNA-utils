# file  : plotting.R
# author: Gavin Ha <gha@bccrc.ca>
#         Dept of Molecular Oncolgy
#         British Columbia Cancer Agency
#         University of British Columbia
# date  : November 7, 2012

#data is the output format of APOLLOH
#cytoBand = {T, F}
#alphaVal = [0,1]
#geneAnnot is a dataframe with 4 columns: geneSymbol, chr, start, stop
#spacing is the distance between each track
plotAllelicRatio <- function(dataIn, chr=NULL, geneAnnot=NULL, spacing=4, ylim=c(0,1), xlim=NULL, cex=0.5){
  #color coding
  #alphaVal <- ceiling(alphaVal * 255); class(alphaVal) = "hexmode"
  lohCol <- c("#00FF00","#006400","#0000FF","#8B0000","#006400","#BEBEBE","#FF0000","#BEBEBE","#FF0000")
  #lohCol <- paste(lohCol,alphaVal,sep="")
  #lohCol <- col2rgb(c("green","darkgreen","blue","darkgreen","grey","red"))
  names(lohCol) <- c("HOMD","DLOH","NLOH","GAIN","ALOH","HET","ASCNA","BCNA","UBCNA")
  
  if (!is.null(chr)){
  	for (i in chr){
    	dataByChr <- dataIn[dataIn[,"Chr"]==i,]
    	dataByChr <- dataByChr[dataByChr[,"APOLLOHcall"]!="OUT",]
    	#plot the data
    	#if (outfile!=""){ pdf(outfile,width=10,height=6) }
    	par(mar=c(spacing,8,2,2),xpd=FALSE)
    	xaxt <- "s"
    	if (is.null(xlim)){
    		xlim <- as.numeric(c(1,dataByChr[dim(dataByChr)[1],"Position"]))
    		xaxt <- "n"
    	}
    	plot(dataByChr[,"Position"],dataByChr[,"AllelicRatio"],col=lohCol[dataByChr[,"APOLLOHcall"]],pch=16,xaxt=xaxt,ylim=ylim,xlim=xlim,xlab="",ylab="Allelic Ratio",cex.lab=1.5,cex.axis=1.5,main=paste("Chromosome ",i,sep=""),cex=cex,las=1)
    	lines(as.numeric(c(1,dataByChr[dim(dataByChr)[1],"Position"])),rep(0.5,2),type="l",col="grey",lwd=0.75)
    
    	if (!is.null(geneAnnot)){
     	plotGeneAnnotation(geneAnnot,i,cex=1)      
    	}
  	}
  }else{  #plot for all chromosomes
    	coord <- getGenomeWidePositions(dataIn[,"Chr"],dataIn[,"Position"])
    	plot(coord$posns,as.numeric(dataIn[,"AllelicRatio"]),col=lohCol[dataIn[,"APOLLOHcall"]],pch=16,xaxt="n",ylim=c(0,1),xlim=c(1,as.numeric(coord$posns[length(coord$posns)])),xlab="Chromosome",ylab="Allelic Ratio",cex.lab=1.5,cex.axis=1.5,cex=cex,las=1,bty="n")
    	lines(as.numeric(c(1,coord$posns[length(coord$posns)])),rep(0.5,2),type="l",col="grey",lwd=3)
    	plotChrLines(dataIn[,"Chr"],coord$chrBkpt,c(-0.1,1.1))
    
  }
  par(xpd=NA)
}

#data is the output format of APOLLOH
#alphaVal = [0,1]
#geneAnnot is a dataframe with 4 columns: geneSymbol, chr, start, stop
#spacing is the distance between each track
plotClonalFrequency <- function(dataIn, chr=NULL, normal=NULL, geneAnnot=NULL, spacing=4, ylim=c(0,1),xlim=NULL, cex=0.5){
  #color coding
  lohCol <- c("#00FF00","#006400","#0000FF","#8B0000","#006400","#BEBEBE","#FF0000","#FF0000","#FF0000")
  names(lohCol) <- c("HOMD","DLOH","NLOH","GAIN","ALOH","HET","ASCNA","BCNA","UBCNA")
  
  #get unique set of cluster and estimates
  #table: 1st column is cluster number, 2nd column is clonal freq
  clusters <- unique(dataIn[,c("ClonalCluster","CellularFrequency")])
  clusters <- clusters[!is.na(clusters[,1]),,drop=F] #exclude NA
  if (!is.null(normal)){ clusters[,2] <- (1-as.numeric(clusters[,2])) * (1-as.numeric(normal)) }
  
	dataToUse <- dataIn[dataIn[,"APOLLOHcall"]!="OUT",]
    dataToUse[dataToUse[,"CellularFrequency"]=="NA"|is.na(dataToUse[,"CellularFrequency"]),c("ClonalCluster","CellularFrequency")] <- c(NA,NA) 
    #extract clonal info
    clonalFreq <- cbind(as.numeric(dataToUse[,"ClonalCluster"]),as.numeric(dataToUse[,"CellularFrequency"]))
    #mode(clonalFreq) <- "numeric"
    clonalFreq[,2] <- 1 - clonalFreq[,2]
    if (!is.null(normal)){ clonalFreq[,2] <- clonalFreq[,2] * (1-normal) }
    clonalFreq[is.na(clonalFreq[,2]) | clonalFreq[,2]=="0" | clonalFreq[,2]=="NA",2] <- 0;
    
  #plot per chromosome
  if (!is.null(chr)){
  for (i in chr){
    ind <- dataToUse[,"Chr"]==as.character(i)
    dataByChr <- dataToUse[ind,]
    clonalFreq <- clonalFreq[ind,]
    #plot the data
    par(mar=c(spacing,8,2,2),xpd=FALSE)
    xaxt <- "s"
    if (is.null(xlim)){
    	xlim <- c(1,as.numeric(dataByChr[dim(dataByChr)[1],"Position"]))
    	xaxt <- "n"
    }
    
    #PLOT CLONAL FREQUENCIES
    plot(dataByChr[,"Position"],clonalFreq[,2],type="h",col=lohCol[dataByChr[,"APOLLOHcall"]],ylim=ylim,pch=16,xaxt=xaxt,xlim=xlim,xlab="",ylab="Cellular Prevalence",cex.lab=1.5,cex.axis=1.5,main=paste("Chromosome ",i,sep=""),cex=cex,las=1)
    
    #plot cluster lines and labels
    for (j in 1:length(clusters[,1])){
      chrLen <- as.numeric(dataByChr[dim(dataByChr)[1],"Position"])
      lines(c(1-chrLen*0.02,chrLen*1.02),rep(clusters[j,2],2),type="l",col="grey",lwd=3)
      #text(x=1-dataByChr[dim(dataByChr)[1],"Position"]*0.02,y=clusters[j,2],label=paste("Z=",clusters[j,1],"",sep=""),pos=2,cex=1.2)
      #text(x=dataByChr[dim(dataByChr)[1],"Position"]*1.02,y=clusters[j,2],label=paste("Z=",clusters[j,1],"",sep=""),pos=4,cex=1.2)
      mtext(side=4,at=clusters[j,2],text=paste("Z",clusters[j,1],"",sep=""),cex=1,padj=0.5,adj=1,las=2,outer=FALSE)
      mtext(side=2,at=clusters[j,2],text=paste("Z",clusters[j,1],"",sep=""),cex=1,padj=0.5,adj=0,las=2,outer=FALSE)
    }
    #points(dataByChr[,"Position"],clonalFreq[,2],col=lohCol[dataByChr[,"APOLLOHcall"]],ylim=c(0,1),pch=16,xaxt="n",xlim=c(1,as.numeric(dataByChr[dim(dataByChr)[1],"Position"])),xlab="",ylab="Cellular Frequency",cex.lab=1,cex.axis=1,main=paste("Chromosome ",i,sep=""),cex=0.5,las=1)
        
    if (!is.null(normal)){
      chrLen <- as.numeric(dataByChr[dim(dataByChr)[1],"Position"])
      lines(c(1-chrLen*0.02,chrLen*1.02),rep((1-normal),2),type="l",col="#000000",lwd=3)
      #text(x=1,y=(1-normal),label=paste("--T--",sep=""),pos=2,cex=1.2)
      #text(x=dataByChr[dim(dataByChr)[1],"Position"],y=(1-normal),label=paste("--T--",sep=""),pos=4,cex=1.2)
      #mtext(side=4,at=(1-normal),text=paste("-T-",sep=""),padj=0.5,adj=1,cex=1,las=2,outer=FALSE)
      #mtext(side=2,at=(1-normal),text=paste("-T-",sep=""),padj=0.5,adj=0,cex=1,las=2,outer=FALSE)
    }

    if (!is.null(geneAnnot)){
		plotGeneAnnotation(geneAnnot,i,cex=1)
    }
   }   
  }  else{  #plot genome-wide
	coord <- getGenomeWidePositions(dataToUse[,"Chr"],dataToUse[,"Position"])
  	plot(coord$posns,clonalFreq[,2],type="h",col=lohCol[dataToUse[,"APOLLOHcall"]],ylim=c(0,1),pch=16,xaxt="n",xlim=c(1,as.numeric(coord$posns[length(coord$posns)])),xlab="Chromosome",ylab="Cellular Prevalence",cex.lab=1.5,cex.axis=1.5,cex=cex,las=1,bty="n")    	
    	plotChrLines(dataIn[,"Chr"],coord$chrBkpt,c(0,1))
    	
    	 #plot cluster lines and labels
    for (j in 1:length(clusters[,1])){
      chrLen <- as.numeric(coord$posns[length(coord$posns)])
      lines(c(1-chrLen*0.02,chrLen*1.02),rep(clusters[j,2],2),type="l",col="grey",lwd=3)
      mtext(side=4,at=clusters[j,2],text=paste("Z",clusters[j,1],"",sep=""),cex=1,padj=0.5,adj=1,las=2,outer=FALSE)
      mtext(side=2,at=clusters[j,2],text=paste("Z",clusters[j,1],"",sep=""),cex=1,padj=0.5,adj=0,las=2,outer=FALSE)
    }
     if (!is.null(normal)){
      chrLen <- as.numeric(coord$posns[length(coord$posns)])
      lines(c(1-chrLen*0.02,chrLen*1.02),rep((1-normal),2),type="l",col="#000000",lwd=3)
      #mtext(side=4,at=(1-normal),text=paste("-T-",sep=""),padj=0.5,adj=1,cex=1,las=2,outer=FALSE)
      #mtext(side=2,at=(1-normal),text=paste("-T-",sep=""),padj=0.5,adj=0,cex=1,las=2,outer=FALSE)
    }
    
  }
  par(xpd=NA)
}



#data is the output format of APOLLOH2.0 (*loh.txt)
#alphaVal = [0,1]
#geneAnnot is a dataframe with 4 columns: geneSymbol, chr, start, stop
#spacing is the distance between each track
plotCNlogRByChr <- function(dataIn, segments=NULL, chr=NULL, geneAnnot=NULL, ylim=c(-4,6), spacing=4, alphaVal=1, xlim=NULL, cex=0.5){
  #color coding
  alphaVal <- ceiling(alphaVal * 255); class(alphaVal) = "hexmode"
  cnCol <- c("#00FF00","#006400","#0000FF","#8B0000","#FF0000","#FF0000")
  cnCol <- paste(cnCol,alphaVal,sep="")
  #cnCol <- col2rgb(c("green","darkgreen","blue","darkred","red","brightred"))
  names(cnCol) <- c("0","1","2","3","4","5")

  if (!is.null(chr)){
  	for (i in chr){
    	dataByChr <- dataIn[dataIn[,"Chr"]==i,]
    	dataByChr <- dataByChr[dataByChr[,"APOLLOHcall"]!="OUT",]
    	#plot the data
   		#if (outfile!=""){ pdf(outfile,width=10,height=6) }
    	par(mar=c(spacing,8,2,2),xpd=FALSE)
    	coord <- as.numeric(dataByChr[,"Position"])
    	xaxt <- "s"
    	if (is.null(xlim)){
    		xlim <- c(1,as.numeric(dataByChr[dim(dataByChr)[1],"Position"]))
    		xaxt <- "n"
    	}
    	plot(coord,as.numeric(dataByChr[,"LogRatio"]),col=cnCol[as.character(dataByChr[,"CopyNumber"])],pch=16,xaxt=xaxt,ylim=ylim,xlim=xlim,xlab="",ylab="Copy Number (log ratio)",cex.lab=1.5,cex.axis=1.5,main=paste("Chromosome ",i,sep=""),cex=cex,las=1)
    	lines(c(1,as.numeric(dataByChr[dim(dataByChr)[1],"Position"])),rep(0,2),type="l",col="grey",lwd=0.75)
    	
		if (!is.null(segments)){
			segments <- segments[segments[,"Chromosome"]==as.character(i),]
			for (j in 1:nrow(segments)){
			lines(c(segments[j,"Start_Position.bp."],segments[j,"End_Position.bp."]),rep(segments[j,"Median_logR"],2),col=cnCol[as.character(segments[,"Copy_Number"])])
			}
		}		    	
  
    	if (!is.null(geneAnnot)){
			plotGeneAnnotation(geneAnnot,i,cex=1)
    	}
     }
    }else{  #plot for all chromosomes
    	coord <- getGenomeWidePositions(dataIn[,"Chr"],dataIn[,"Position"])
    	plot(coord$posns,as.numeric(dataIn[,"LogRatio"]),col=cnCol[as.character(dataIn[,"CopyNumber"])],pch=16,xaxt="n",ylim=ylim,xlim=c(1,as.numeric(coord$posns[length(coord$posns)])),xlab="Chromosome",ylab="Copy Number (log ratio)",cex.lab=1.5,cex.axis=1.5,cex=cex,las=1,bty="n")
    	lines(as.numeric(c(1,coord$posns[length(coord$posns)])),rep(0,2),type="l",col="grey",lwd=2)
    	plotChrLines(dataIn[,"Chr"],coord$chrBkpt,ylim)
    }
    par(xpd=NA)
}  


plotGeneAnnotation <- function(geneAnnot,chr=1,...){
    colnames(geneAnnot) <- c("Gene","Chr","Start","Stop")
    geneAnnot <- geneAnnot[geneAnnot[,"Chr"]==as.character(chr),]
    if (nrow(geneAnnot) != 0){
		for (g in 1:dim(geneAnnot)[1]){
			#print(geneAnnot[g,"Gene"])
			abline(v=as.numeric(geneAnnot[g,"Start"]),col="black",lty=3,xpd=F)
			abline(v=as.numeric(geneAnnot[g,"Stop"]),col="black",lty=3,xpd=F)			
			atP <- (as.numeric(geneAnnot[g,"Stop"]) - as.numeric(geneAnnot[g,"Start"]))/2 + as.numeric(geneAnnot[g,"Start"])
			#if (atP < dataByChr[1,2]){ 
			#	atP <- dataByChr[1,2] 
			#}else if (atP > dataByChr[dim(dataByChr)[1],2]){ 
			#	atP <- dataByChr[dim(dataByChr)[1],2] 
			#}
			mtext(geneAnnot[g,"Gene"],side=3,line=0,at=atP,...)
		}
	}
}

plotChrLines <- function(chrs,chrBkpt,yrange){
	#plot vertical chromosome lines
    	for (j in 1:length(chrBkpt)){
    		lines(rep(chrBkpt[j],2),yrange,type="l",lty=2,col="black",lwd=0.75)
    	}
    	numLines <- length(chrBkpt)
    	mid <- (chrBkpt[1:(numLines-1)]+chrBkpt[2:numLines])/2
    	chrs[chrs=="X"] <- 23; chrs[chrs=="Y"] <- 24;
      chrsToShow <- sort(unique(as.numeric(chrs)))
      chrsToShow[chrsToShow==23] <- "X"; chrsToShow[chrsToShow==24] <- "Y";
    	axis(side=1,at=mid,labels=c(chrsToShow),cex.axis=1.5,tick=FALSE)
}

getGenomeWidePositions <- function(chrs,posns){  
  #create genome coordinate scaffold
	positions <- as.numeric(posns)
	chrsNum <- unique(chrs)
	chrBkpt <- rep(0,length(chrsNum)+1)
	for (i in 2:length(chrsNum)){
    	chrInd <- which(chrs==chrsNum[i])
    	prevChrPos <- positions[chrInd[1]-1]      
    	chrBkpt[i] = prevChrPos
    	positions[chrInd] = positions[chrInd] + prevChrPos
	}
	chrBkpt[i+1] <- positions[length(positions)]
	return(list(posns=positions,chrBkpt=chrBkpt))
}