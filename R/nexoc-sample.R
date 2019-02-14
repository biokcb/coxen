### Example script on how to use the NEXOC function... which I guess is really just a coxen pipeline function, but whatever

#set working directory
setwd("/Users/Kristen/Desktop/COXEN")

#input directory where the coxen set is located
input.dir <- "./Laval"

#input directory where the independent test set is located
indep.dir <- "./BladderFive"

#load necessary librarys
library(affy)
library(frma)
library(hgu133plus2frmavecs) # for processing expression data

library(MASS)
library(MiPP)
library(samr)
library(qvalue) #t-test
library(SuperLearner)
library(randomForest) #for deg selection, model building

#source nexoc.R script, which contains copies of all the functions we'll need that aren't in packages
source("./nexoc.R")

#load up gene expression data for NCI60 
load("./NCI60 Batch Effects/NCI60.fRMA.collapsed.rda")
genematrix <- frozen.rma.collapse
rm(frozen.rma.collapse)

# read in .CEL files or processed expression matrix
filenames <- list.files(path=input.dir, pattern=".CEL")
eset <- ReadAffy(filenames=paste(input.dir,'/',filenames,sep=''), phenoData=NULL)
slotNames(eset)
dates <- eset@protocolData
slotNames(dates)
date.info <- dates@data
date.info # all were scanned on the same date.
coxenmat <- frma(eset)
coxenmat <- exprs(coxenmat)
rm(eset)

# alternatively, load the processed laval set
load("./Laval Batch/Laval.FrozenRMA.rda")
coxenmat <- exprs(laval.frozen)
rm(laval.frozen)

# read in .CEL files or processed expression matrix
filenames <- list.files(path=indep.dir, pattern=".CEL")
eset <- ReadAffy(filenames=paste(indep.dir,'/',filenames,sep=''), phenoData=NULL)
slotNames(eset)
dates <- eset@protocolData
slotNames(dates)
date.info <- dates@data
date.info # all were scanned on the same date.
indepmat <- frma(eset)
indepmat <- exprs(indepmat)
rm(eset)

# read in drug data - here we are just loading the FDA approved list I had previously saved. 
load("./meanac_REG0.1alldrugs.rda") # the object here is named "meanac"
identical(colnames(meanac), colnames(genematrix))

# for the 76 FDA approved drug list, I also have group sizes that were selected 
# based on maximizing DEGs, so we'll use that too. 
num.group <- read.csv("./NCI60frma-FDA76groups.csv", header=TRUE, row.names=1)
num.group <- num.group[order(num.group[,1]),]
identical(as.character(num.group[,1]), rownames(meanac))

# regression coefficients calculated originally
reg.coef <- read.csv("./Regression.Coefficients.Cutoff_0.1.csv", header=TRUE, row.names=1)
reg.coef <- reg.coef[order(reg.coef[,1]),]
identical(as.character(reg.coef[,1]), rownames(meanac))

#randomForests using the ivDrug method requires tissue type information
info <- read.csv("./NCI60 Batch Effects/filename.match.csv", header=TRUE, row.names=1)
tissue <- info[match(colnames(meanac), info[,4]),c(4,1)]
rownames(tissue) <- tissue[,1]

#### applying the nexoc function. ####

#create empty master matrix to store predictions and other information
master.pred <- matrix(data=NA, nrow=nrow(meanac), ncol=(ncol(coxenmat) + ncol(indepmat) + 10))
colnames(master.pred) <- c("NSC","Reg.Coef","Average.Score","Num.GeneModels", "Average.NumGenes","T.Test", "SAM","SAM-FDR", "CorTest", "COXEN", colnames(coxenmat), colnames(indepmat))
master.pred[,1] <- rownames(meanac)
master.pred[,2] <- reg.coef[,2]

# only indep set predictions
master.pred <- matrix(data=NA, nrow=nrow(meanac), ncol=(ncol(indepmat) + 10))
colnames(master.pred) <- c("NSC","Reg.Coef","Average.Score","Num.GeneModels", "Average.NumGenes","T.Test", "SAM","SAM-FDR", "CorTest", "COXEN", colnames(indepmat))
master.pred[,1] <- rownames(meanac)
master.pred[,2] <- reg.coef[,2]


for (k in 1:nrow(meanac)){
  master.pred[k,3:ncol(master.pred)] <- nexoc(genematrix, coxenmat, indepmat, scale=TRUE, pred.all=TRUE, rownames(meanac[k,]), meanac[k,], rev=TRUE, tissue, num=num.group[k,2:3], deg=c("t.test"), min.gene=10, max.gene=1000, method=c("randomForests"))
}

write.csv(master.pred, file="./BladderFDAPredictions-mipp.csv")


#### group extremes
for (i in 1:nrow(meanac)){
  extremes <- select.extremes(genematrix, order(meanac[i,]), 0.1, min.num.sensitive=8, min.num.resistant=8, max.num.resistant=15, max.num.sensitive=15) 
  num.group[i,2] <- length(extremes$drug.sensitive) 
  num.group[i,3] <- length(extremes$drug.resistant)
  num.group[i,1] <- rownames(meanac[i,])
}