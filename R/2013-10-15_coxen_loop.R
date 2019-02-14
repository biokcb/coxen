## This script was modified for an automated, general run of Jared's COXEN runs. The goal is to look at many COXEN models for..... 1 - a set of various drugs, 2 - different numbers of dogs in sensitive/resistant groups in the ACC29, 3 - different correlation score cutoffs (90,95,98), 4 - Using the ACC29 to fine tune the model & using the NCI60 panel only.

#set the appropriate working directory where all the required files are & load ALL the libraries to start. 
setwd("/Users/Kristen/Documents/GustafsonLab/Jared/")
library(affy)
library(MiPP)
library(e1071)
library(MASS)

#create a vector with the different drug shorthand names, needs to be a matrix to work in the following steps..
druglist<-c("vbl")
druglist<-as.matrix(druglist)
probability<-c(0.98,0.95,0.9)
rulelist<-c("lda","qda","svmlin","svmrbf","logistic") #to run all the models

#creating lists with file names to load based on drug name. 
x0<-1:length(druglist)  
NCIRMA_file<-sprintf("%s_%s.%s","NCI60RMA",druglist[x0,],"rda") # NCI60RMA_dox.rda
HumProbeTable<-sprintf("%s_%s%s.%s","hum",druglist[x0,],"psid","txt") # hum_doxpsid.txt
geneTable<-sprintf("%s%s.%s",druglist[x0,],"gensym","txt") # doxgenesym.txt
DogProbeTable<-sprintf("%s_%s%s.%s","dog",druglist[x0,],"psid","txt") # dog_doxpsid.txt
ACCdrugTable<-sprintf("%s%s%s.%s","ACC29",druglist[x0,],"rank","txt") #ACC29doxrank.txt

### Beginning of the long loop. First we start by loading the appropriate files for each drug, then running the individual analysis for each drug one at a time. i will indicate we are working with one drug at a time per iteration. In a more general code, we would load the NCI60 RMA data then sort it based on each drug sensitivity profile instead of loading pre-sorted data each time.

for (i in 1:length(druglist)){
	load("ACC29RMA.rda")
	load(NCIRMA_file[i])
	drugs<-druglist[i,]
	print(drugs)
	A <- read.table(file=HumProbeTable[i], sep="\t", header=FALSE)
	A <- A[,1]
	A <-as.vector(A)
	HumIDs <- NCI60RMA[A,]

##to replace probesetIDs with gene symbols for correlation purposes##
GenSym <- read.table(file=geneTable[i],sep="\t",header=FALSE)
GenSym <- GenSym[,1]
GenSym <- as.vector(GenSym)
HumIDs <- t(HumIDs)
colnames(HumIDs) <- GenSym

##to acquire ACCRMA values from probesets of interest##
##make sure the canine probeset ID names are in correct format.  ex:  "CfaAffx.20329.1.S1_at" ##

yourProbeSetIDs <- read.table(file=DogProbeTable[i],sep="\t",header=FALSE)
A <- yourProbeSetIDs[,1]
A <- as.vector(A)
	if (any(drugs==c("lom","ptx","vbl"))==TRUE){
		ACC29RMA<-ACC29RMA[,-c(7,24)]
	}
K9IDs <- ACC29RMA[A,]
k9length<-ncol(K9IDs)

##to replace probesetIDs with gene symbols for correlation purposes##
K9IDs <- t(K9IDs)
colnames(K9IDs) <- GenSym

##to reorder ACC29RMA based on drug sensitivity##

ACCdrugdata <- read.table(file=ACCdrugTable[i],header=TRUE, dec=".") 
#numACCdrug<-ACCdrugdata[-(charmatch("resistant",ACCdrugdata, nomatch=0)),]
#numACCdrug<-as.vector(numACCdrug)
#maxnum<-max(numACCdrug)
#("resistant" %in% ACCdrugdata)<-max(ACCdrugdata[is.numeric(ACCdrugdata)])
#ACCdrugdata[charmatch("resistant",ACCdrugdata)]<-250000
ACCdrugdata <- as.matrix(ACCdrugdata[,1])
K9IDs <- cbind(ACCdrugdata,K9IDs)
K9IDs <- K9IDs[order(K9IDs[,1]),]
K9IDs <- K9IDs[,2:ncol(K9IDs)] #remove drug data

##scale the NCI60RMA and ACC29RMA data. Columns should be gene symbols.

K9IDs <- scale(K9IDs)
HumIDs <- scale(HumIDs)
cutoff<-nrow(HumIDs)-12 + 1
HumIDs.GH <- HumIDs[c(1:12,cutoff:nrow(HumIDs)),] #split into sen/res groups

###########Step5 of the COXEN analysis##############

##generate original correlation matrices in NCI60RMA and ACCRMA##

CMat1 <- cor(HumIDs.GH)
CMat2 <- cor(K9IDs)

##generate 2nd correlation matrix between datasets##

CMat1 <- as.matrix(CMat1)
CMat2 <- as.matrix(CMat2)
CMat1m <- CMat1
is.na(diag(CMat1m)) <- TRUE
CMat2m <- CMat2
is.na(diag(CMat2m)) <- TRUE
COR <- cor(CMat1m,CMat2m,use="pairwise.complete")
corfile<-sprintf("%s_%s.%s","COR",drugs,"csv")
write.csv(COR,file=corfile)
RCout <- diag(cor(CMat1m,CMat2m,use="pairwise.complete"))
rcoutfile<-sprintf("%s_%s.%s","RCout",drugs,".rda")
save(RCout,file=rcoutfile)
PermRC <- c()
for(i in 1:1000)
{
order <- sample(seq(1:length(A)), length(A), replace = FALSE)
Perm1 <- CMat1m[order,order]
RCO <- diag(cor(Perm1,CMat2m,use="pairwise.complete"))
PermRC <- cbind(PermRC,RCO)
rm(order,Perm1,RCO)
}
HST<-hist(PermRC,main="Histogram of PermRC for COR")
#begin the loop for the different correlation score cutoffs. 
for (p in 1:3){
	probs<-probability[p]
	PermRC_prob<-quantile(PermRC,prob=probs)
	RCout_prob<- RCout[RCout>PermRC_prob]
	ModelGenes<-names(RCout_prob)
	modelgenefile<-sprintf("%s.%s.%s_%s","ModelGenes",probs,drugs,"NCI60_ACC30.csv")
	write.csv(ModelGenes,file=modelgenefile)
if (length(ModelGenes)>5){
##################### Step 6 of COXEN analysis #################################

##First need to assemble RMA expression data for both datasets with your selected genes for modeling##
B <- ModelGenes
B <- as.vector(B)
HumMiPP <- HumIDs[,B]

#### subsetting the NCI60 panel into the 12 most sensitive and resistant groups#########
HumMiPP_split <- HumMiPP[c(1:12,cutoff:nrow(HumMiPP)),]
HumMiPP_split <- t(HumMiPP_split) #genes x cells

######extract model genes from K9IDs object.  Use all the cell lines for model building for this round
K9MiPP <- K9IDs[,B]
K9MiPP <-t(K9MiPP)

#####subset K9MiPP as sen/res groups. Using j as the # of cell lines in the cutoff, then including anything with a GI50 value within 10% of that cutoff. Using groups of 5-12
ACCdrugdata<-ACCdrugdata[order(ACCdrugdata)]

for (w in 5:12){
	m <- w #set initial values for top/bottom cutoff
	n <- k9length - w + 1
	totlen<-k9length+1
	
	while (ACCdrugdata[m+1] <= ACCdrugdata[w] + ACCdrugdata[w]*0.1){
		m = m + 1
	}
	
	while (ACCdrugdata[n-1] >= ACCdrugdata[totlen-w]-ACCdrugdata[totlen-w]*0.1){
		n = n - 1
	}

	K9MiPP_split <- K9MiPP[,c(1:m,n:k9length)]

for (k in 1:length(rulelist)){ #within this loop: using the ACC as a test set for fine tuning
	rule=rulelist[k]
	print(c("NCI60 & ACC29",drugs,rule,probs))
	MiPP <- cbind(HumMiPP_split,K9MiPP_split)
	MiPP <-as.matrix(MiPP)

	x.train<-MiPP[,1:24]
	y.train<-(c(rep("Sen",12),rep("Res",12)))
	x.test<-MiPP[,25:ncol(MiPP)]
	y.test<-(c(rep("Sen",m),rep("Res",(totlen-n)))) 

	tryCatch({
		out<-c()
		out<-mipp.seq(x=x.train,y=y.train,x.test=x.test,y.test=y.test,n.fold=5,percent.cut=1,rule=rule)}, error=function(e){})

	if (is.null(out)==FALSE){ #if there's an output from the model, save it! then pull out individual predictions
		
		out.model <- out$model
		modelfile<-sprintf("%s_%s.%s.%s.%s_%s.%s", "outmodel",drugs,probs,w,"NCI60_ACC30",rule,"csv")
		write.csv(out.model,file=modelfile)

# picking out the most parsimonious model with the lowest error rate out of all the sequences mipp runs through (0-3) 
comparemodels<- out.model[out.model$Select=="**",]
comparemodels <- comparemodels[order(comparemodels$Te.ER),]
sequence <- comparemodels[1,"Seq"]
#pulling out the sequence of genes for said model
seq <- out.model[out.model$Seq==sequence,]
l <- charmatch("**",seq$Select)
seq <- seq[1:l,]
	geneseq <- seq$Gene
	geneseq <-as.vector(geneseq)
	
#pulling out individual predictions based on which rule was used to make the model	
	if (rule=="lda"){
		colnames(x.train)<-c(1:ncol(x.train))
		colnames(x.test)<-c(1:ncol(x.test))
		x.train<-t(x.train)
		x.train<-as.matrix(x.train[,geneseq])
		x.test<-t(x.test)
		x.test<-as.matrix(x.test[,geneseq])
		fit<-lda(x.train,y.train)
		output<-predict(fit,x.test)	
	}
	else if (rule=="qda"){
		colnames(x.train)<-c(1:ncol(x.train))
		colnames(x.test)<-c(1:ncol(x.test))
		x.train<-t(x.train)
		x.train<-as.matrix(x.train[,geneseq])
		x.test<-t(x.test)
		x.test<-as.matrix(x.test[,geneseq])
		fit<-qda(x.train,y.train)
		output<-predict(fit,x.test)	
	}
	else if (rule=="logistic"){
		x.train<-t(x.train)
		x.test<-t(x.test)
		y.train <- factor(y.train); levels(y.train) <- c("1","0")
		y.test <- factor(y.test); levels(y.test) <- c("1","0")
		if(is.data.frame(x.train)) x.train <- as.matrix(x.train)
		if(is.data.frame(x.test)) x.test <- as.matrix(x.test)
		fit <- glm(y.train ~ x.train[,geneseq], family="binomial")
		predx <- cbind(1,x.test[,geneseq])%*%t(matrix(fit$coef, nrow=1))
		prob <- 1/(1+exp(-predx))
		postdf <- data.frame(prob, y.test)
		post.prob <- ifelse(postdf$y.test=="1", 1-postdf$prob, postdf$prob)
		ind <- ifelse(post.prob > .5, 1, 0)
		output <-cbind(postdf,post.prob,ind)
	}
	else if (rule=="svmlin"){
		x.train <- t(x.train)
		x.test <- t(x.test)
		if(is.data.frame(x.train)) x.train <- as.matrix(x.train)
		if(is.data.frame(x.test)) x.test <- as.matrix(x.test)
		y.train <- factor(y.train)
		y.test <- factor(y.test)
		fit <- svm(x.train[,geneseq],y.train, kernel="linear")
		True.class <- y.test
		Pred.class <- predict(fit, x.test[,geneseq])
		fofx <- numeric(length(y.test))
		for(h in 1:length(y.test)){
			xin <- x.test[h,]
			fofx[h] <- linearkernel.decision.function(xin, x.train, fit)
			}
		c <-1 #optimal parameter
		prob <- 1/(1+c*exp(-fofx))
		postdf <- data.frame(prob, True.class)
		post.prob <- ifelse(postdf$True.class==Pred.class, 1-postdf$prob,postdf$prob)
		N <- length(y.test)
		nMiss <- N - sum(True.class==Pred.class)
		Er <- nMiss/N
		MiPP <- sum(post.prob)-nMiss
		sMiPP <- MiPP/N
		output<-cbind(Pred.class,postdf,post.prob)
	}
	else { #the last possible case, svmrbf
		x.train <- t(x.train)
		x.test <- t(x.test)
		if(is.data.frame(x.train)) x.train <- as.matrix(x.train)
		if(is.data.frame(x.test)) x.test <- as.matrix(x.test)
		y.train <- factor(y.train)
		y.test <- factor(y.test)
		gammap <- 1/length(ncol(x.train))
		fit <- svm(x.train[,geneseq], y.train, kernel="radial", gamma=gammap)
		True.class <- y.test
		Pred.class <- predict(fit, x.test[,geneseq])
		fofx <- numeric(length(y.test))
		for(h in 1:length(y.test)){
			xin <- x.test[h,]
			fofx[h] <- rbfkernel.decision.function(xin, x.train, fit)
		}
		c <-100 #optimal parameter
		prob <- 1/(1+c*exp(-fofx))
		postdf <- data.frame(prob, True.class)
		post.prob <- ifelse(postdf$True.class==Pred.class, 1-postdf$prob,postdf$prob)
		N <- length(y.test)
		nMiss <- N - sum(True.class==Pred.class)
		Er <- nMiss/N
		MiPP <- sum(post.prob)-nMiss
		sMiPP <- MiPP/N
		output <-cbind(Pred.class,postdf,post.prob)
	}
		outfile<-sprintf("%s_%s_%s.%s_%s_%s.%s","predict",drugs,probs,w,"NCI60_ACC30",rule,"csv")
		write.csv(output,file=outfile)	
	} 
		rm(output)
		rm(outfile)
		rm(out.model)
 #need to clear out the prior looped output before running it again in case the mipp.seq returns nothing, then out could contain data from the prior run

} #Models with ACC29 as the test set 
}

for (k in 1:length(rulelist)){ #within this loop: NCI60 only used as training/test set.
	rule <- rulelist[k]
	print(c("NCI60 only",drugs,rule,probs))
	x <- HumMiPP_split
	y <- factor(c("S","S","S","S","S","S","S","S","S","S","S","S","R","R","R","R","R","R","R","R","R","R","R","R"))
	out<-mipp(x=x,y=y,n.fold=5,p.test=1/3,n.split=40,n.split.eval=100,percent.cut=1,rule=rule)
	out.model<-out$model
	out.model.eval<-out$model.eval
	colnames(out.model.eval)<-gsub(" ",".",colnames(out.model.eval))
	
	outfile1<-sprintf("%s_%s.%s.%s_%s_%s.%s","outmodel",drugs,probs,"all","NCI60only",rule,"csv")
	outfile2<-sprintf("%s_%s.%s.%s_%s_%s.%s","outmodel_eval",drugs,probs,"all","NCI60only", rule,"csv")
	
	write.csv(out.model,file=outfile1)
	write.csv(out.model.eval,file=outfile2)

	#sort table by mean sMiPP score
	out.model.eval<-out.model.eval[order(out.model.eval$mean.sMiPP, decreasing=TRUE),]
	genesonly<-out.model.eval[,2:(ncol(out.model.eval)-9)]


	for (t in 1:nrow(genesonly)){ #sort the gene sequences numerically
		genesonly[t,]<-genesonly[t,order(genesonly[t,])]
	}
	droprows<-c() #list of dupilcated rows that need to be dropped from the list..
	genesonly[is.na(genesonly)]<-0
	for (v in 1:(nrow(genesonly)-1)){
			for (u in 2:nrow(genesonly)){
				if (u!=v && all((genesonly[v,])==(genesonly[u,]))){
				droprows<-c(droprows,u)
				}
			}
	} 
	droprows<-unique(droprows)
	genesonly<-as.matrix(genesonly)
	genesonly<-t(genesonly)
	if (is.null(droprows==FALSE)){
	genesonly<-genesonly[-droprows,]
	}
	#pull out individual predictions on ACC29 from these gene models...

	for (q in 1:5){
		
	MiPP2 <- cbind(HumMiPP_split,K9MiPP)
	MiPP2<-as.matrix(MiPP2)
	x.train<-MiPP2[,1:24]
	y.train<-(c(rep("Sen",12),rep("Res",12)))
	x.test<-MiPP2[,25:ncol(MiPP2)]	
		if (any(drugs==c("lom","ptx", "vbl"))==TRUE){
			y.test<-(c(rep("Sen",13),rep("Res",14)))
		}
		else y.test<-(c(rep("Sen",14),rep("Res",15)))
		
		geneseq<-genesonly[,q]
		geneseq<-as.vector(geneseq)
		geneseq<-geneseq[geneseq!=0]
		
		if (rule=="lda"){
			colnames(x.train)<-c(1:ncol(x.train))
			colnames(x.test)<-c(1:ncol(x.test))
			x.train<-t(x.train)
			x.train<-as.matrix(x.train[,geneseq])
			x.test<-t(x.test)
			x.test<-as.matrix(x.test[,geneseq])
			fit<-lda(x.train,y.train)
			predict<-predict(fit,x.test)	
		}
		else if (rule=="qda"){
			colnames(x.train)<-c(1:ncol(x.train))
			colnames(x.test)<-c(1:ncol(x.test))
			x.train<-t(x.train)
			x.train<-as.matrix(x.train[,geneseq])
			x.test<-t(x.test)
			x.test<-as.matrix(x.test[,geneseq])
			fit<-qda(x.train,y.train)
			predict<-predict(fit,x.test)	
		}
		else if (rule=="svmlin"){
			x.train <- t(x.train)
			x.test <- t(x.test)
			if(is.data.frame(x.train)) x.train <- as.matrix(x.train)
			if(is.data.frame(x.test)) x.test <- as.matrix(x.test)
			y.train <- factor(y.train)
			y.test <- factor(y.test)
			fit <- svm(x.train[,geneseq],y.train, kernel="linear")
			True.class <- y.test
			Pred.class <- predict(fit, x.test[,geneseq])
			fofx <- numeric(length(y.test))
			for(h in 1:length(y.test)){
				xin <- x.test[h,]
				fofx[h] <- linearkernel.decision.function(xin, x.train, fit)
			}
			c <-1 #optimal parameter
			prob <- 1/(1+c*exp(-fofx))
			postdf <- data.frame(prob, True.class)
			post.prob <- ifelse(postdf$True.class==Pred.class, 1-postdf$prob,postdf$prob)
			N <- length(y.test)
			nMiss <- N - sum(True.class==Pred.class)
			Er <- nMiss/N
			MiPP <- sum(post.prob)-nMiss
			sMiPP <- MiPP/N
			predict <-cbind(Pred.class,postdf,post.prob)
		}
		else if (rule=="svmrbf"){
			x.train <- t(x.train)
			x.test <- t(x.test)
			if(is.data.frame(x.train)) x.train <- as.matrix(x.train)
			if(is.data.frame(x.test)) x.test <- as.matrix(x.test)
			y.train <- factor(y.train)
			y.test <- factor(y.test)
			gammap <- 1/length(ncol(x.train))
			fit <- svm(x.train[,c(geneseq)], y.train, kernel="radial", gamma=gammap)
			True.class <- y.test
			Pred.class <- predict(fit, x.test[,c(geneseq)])
			fofx <- numeric(length(y.test))
			for(h in 1:length(y.test)){
				xin <- x.test[h,]
				fofx[h] <- rbfkernel.decision.function(xin, x.train, fit)
			}
			c <-100 #optimal parameter
			prob <- 1/(1+c*exp(-fofx))
			postdf <- data.frame(prob, True.class)
			post.prob <- ifelse(postdf$True.class==Pred.class, 1-postdf$prob,postdf$prob)
			N <- length(y.test)
			nMiss <- N - sum(True.class==Pred.class)
			Er <- nMiss/N
			MiPP <- sum(post.prob)-nMiss
			sMiPP <- MiPP/N
			predict <-cbind(Pred.class,postdf,post.prob)
		}
		else { #rule = logistic
			x.train<-t(x.train)
			x.test<-t(x.test)
			y.train <- factor(y.train); levels(y.train) <- c("1","0")
			y.test <- factor(y.test); levels(y.test) <- c("1","0")
			if(is.data.frame(x.train)) x.train <- as.matrix(x.train)
			if(is.data.frame(x.test)) x.test <- as.matrix(x.test)
			fit <- glm(y.train ~ x.train[,geneseq], family="binomial")
			predx <- cbind(1,x.test[,geneseq])%*%t(matrix(fit$coef, nrow=1))
			prob <- 1/(1+exp(-predx))
			postdf <- data.frame(prob, y.test)
			post.prob <- ifelse(postdf$y.test=="1", 1-postdf$prob, postdf$prob)
			ind <- ifelse(post.prob > .5, 1, 0)
			predict <-cbind(postdf,post.prob,ind)
		}
		
		predictoutfile<-sprintf("%s%s_%s.%s.%s_%s_%s.%s",rule,"predict",drugs,probs,"all","NCI60only_ACC30_model",q,"csv")
		write.csv(predict,file=predictoutfile) 
			
	}	
} #Models only using NCI60 for training/testing
}
else {
	rm(PermRC_prob)
	rm(RCout_prob)
	rm(ModelGenes)}
}
}

