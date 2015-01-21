#library(TitanCNA)
source("~/software/code/R_SVN/TitanCNA/R/utils.R")
args <- commandArgs(TRUE)

options(stringsAsFactors=FALSE)

titanfile <- args[1]
paramObj <- args[2]
scale <- as.numeric(args[3])
data.type <- args[4]
outfile <- args[5]
symmetric <- TRUE

titan <- read.delim(titanfile,header=T,sep="\t")
load(paramObj)

outputModelParameters(convergeParams, titan, outfile, S_Dbw.scale = scale)
		
