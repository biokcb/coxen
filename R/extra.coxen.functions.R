#Extra COXEN related functions 
# Kristen Brown

# Collapse repeated probesets/RNAs
#collapse.probes<-function(genematrix,method="max",)

# Convert probe IDs to gene Names for Human-Canine COXEN
#gene.convert<-function()

# Select max differentially expressed genes using SAMR
degs.select<-function(genematrix,fdr=0.1,method="hi.lo", num.sen=12, num.res=12, max.gene=300){
		Y <- c(rep(2,num.sen), rep(1,num.res))
		samdata<-list(x=genematrix, y=Y, geneid=as.character(1:nrow(genematrix)), genenames=row.names(genematrix), logged2=TRUE)
		samr.pac<-samr(samdata, resp.typ="Two class unpaired", nperms=500)

		#pick out delta value closest to 0.1
		delta.table <- samr.compute.delta.table(samr.pac)
		colnames(delta.table) <- gsub(" ",".",colnames(delta.table))
		delta.table<-delta.table[complete.cases(delta.table),]
    delta.table<-delta.table[delta.table[,5]!=0,]
		delta.FDR <- which(abs(delta.table[,5] - 0.1)==min(abs(0.1 - delta.table[,5])))	
		del.max <- delta.table[delta.FDR,1]
		del.min <- delta.table[(delta.FDR-1),1]
		cat(del.min,"to ", del.max)
		delta.table <- samr.compute.delta.table(samr.pac,dels=seq(del.min, del.max, 0.002))
		delta.table2 <- delta.table[delta.table[,5] <= 0.1,]
		delta.FDR <- which(abs(delta.table[,5] - 0.1)==min(abs(0.1 - delta.table[,5])))
    if (length(delta.FDR)>1){
      delta.FDR.pick <- which(abs(delta.table[delta.FDR[1],5])==max(abs(delta.table[delta.FDR[2],5])))
      delta.FDR <- delta.FDR[delta.FDR.pick]
    }
		delta.opt <- delta.table[delta.FDR,1]

		#DEGs (Differentially Expressed Genes) Summary.
		siggenes.table<-samr.compute.siggenes.table(samr.pac, delta.opt, samdata, delta.table, all.genes=FALSE)
		
		if (method == "hi.lo"){
			split <- max.gene/2
			up <- siggenes.table$genes.up
			down <- siggenes.table$genes.lo
			if (nrow(up)>split){
				up <- up[1:split,]
				up.tot <- split
			} else up.tot <- nrow(up)
			down.max <- max.gene - up.tot
			if (nrow(down)>down.max){
				split2 <- max.gene - down.max
				down <- down[1:split2,]
				down.tot <- split2
			} 
			DEGs<-rbind(up,down) #bind high/low table rows
		}
		else if (method == "lo"){
			DEGs <- siggenes.table$genes.lo 
		}
		else {
			DEGs<-siggenes.table$genes.up 
		}
		return(DEGs)

}


# Calculate Regression coefficients after removing extremes
reg.coef.calc <- function(drugdata,num.sen=12,num.res=12){
	drug <- as.matrix(drugdata)
	n.cell <- length(drug)
	if (ncol(drug)>1){
	drug <- t(drug)
	}
	drug <- as.matrix(drug[order(drug),])
	if (sd(drug)!=0){
	drug<-scale(drug)
	aa= as.matrix(drug[(num.sen+1):(n.cell-num.res),])
	bb=scale(1:(n.cell-num.sen-num.res))
	reg.coef = lm( aa ~ bb)$coefficients[[2]]
	} else reg.coef = 0
	return(reg.coef)
}

#coexp

coexpression<-function(trainmat,coxenmat,prob=0.98, gene.num=gene.num){

CMat1 <- cor(t(trainmat))
CMat2 <- cor(t(coxenmat))

##generate 2nd correlation matrix between datasets##

CMat1 <- as.matrix(CMat1)
CMat2 <- as.matrix(CMat2)
CMat1m <- CMat1
is.na(diag(CMat1m)) <- TRUE
CMat2m <- CMat2
is.na(diag(CMat2m)) <- TRUE
COR <- cor(CMat1m,CMat2m,use="pairwise.complete")
RCout <- diag(cor(CMat1m,CMat2m,use="pairwise.complete"))
PermRC <- c()
for(j in 1:500)
{
order <- sample(seq(1:gene.num), gene.num, replace = FALSE)
Perm1 <- CMat1m[order,order]
RCO <- diag(cor(Perm1,CMat2m,use="pairwise.complete"))
PermRC <- cbind(PermRC,RCO)
rm(order,Perm1,RCO)
cat(j, '...')
}
PermRC_prob<-quantile(PermRC,prob=prob)
RCout_prob<- RCout[RCout>PermRC_prob]
ModelGenes<-names(RCout_prob)

return(ModelGenes)
}

# LOOCV

# LOHCV

# Percentile Rank 

newscale <- function(scores, min = -1, max = 1){
	old.min <- min(scores)
	old.max <- max(scores)
	new.scores <- c()	
	for (i in 1:length(scores)){
		new.score <- (max - min)/(old.max - old.min)*(scores[i]-old.min) + min
		new.scores[i]<-new.score
	}
	return(new.scores)
}
