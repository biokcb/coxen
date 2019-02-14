############################################################################
#                                                                          #
# SWOG Trial COXEN Score calculation script using equations from the model #
#                                                                          #
############################################################################

# Replace the directory here with the directory where the files to be analyzed are stored, don't add a '/' to the end folder.
input.dir <- "/Users/Kristen/Documents/GustafsonLab"

# Replace the directory here with the directory where the files should be saved after analysis, don't add a '/' to the end folder.
save.dir <- "/Users/Kristen/Documents/GustafsonLab"

# Packages that need to be loaded prior to running the script
library(affy)

# RMA Processed file to be analyzed as a .csv file. MDAnderson, PilotStudy
input.rma.c <- "PilotRMA.csv"

# Data set name, for saving purposes
data <- "PilotStudy"

# Number of years for probability calculation
time <- 5

# gene model using biomarkers (TRUE) or gene names (FALSE)?
gen.bio <- TRUE

#or UnigeneID?
uniID <- FALSE

###### Single drug score calculations ########

# Read in the biomarker information file that Jared created
drug.data <- read.csv(file=paste(input.dir,'/',"MVACGdata.csv",sep=''),header=TRUE)

if (uniID==TRUE){
	gene.data <- read.csv(file=paste(input.dir,'/',"MVACGunigenedata.csv",sep=''),header=TRUE)
	M.gene <- as.vector(gene.data[2:21,1])
	V.gene <- as.vector(gene.data[2:61,3])
	A.gene <- as.vector(gene.data[2:41,5])
	C.gene <- as.vector(gene.data[2:36,7])
	G.gene <- as.vector(gene.data[2:36,9])
	} else {
	gene.data <- read.csv(file=paste(input.dir,'/',"MVACGgenedata.csv",sep=''),header=TRUE)
	M.gene <- as.vector(gene.data[2:18,1])
	V.gene <- as.vector(gene.data[2:59,3])
	A.gene <- as.vector(gene.data[2:36,5])
	C.gene <- as.vector(gene.data[2:35,7])
	G.gene <- as.vector(gene.data[2:33,9])
	}

# Methotrexate variables
M.biom <- as.vector(drug.data[2:21,1])
M.int <- as.vector(drug.data[1,2])
M.coef <- as.vector(drug.data[2:21,2])

# Vinblastine variables
V.biom <- as.vector(drug.data[2:61,3])
V.int <- as.vector(drug.data[1,4])
V.coef <- as.vector(drug.data[2:61,4])

# Adriamycin variables
A.biom <- as.vector(drug.data[2:41,5])
A.int <- as.vector(drug.data[1,6])
A.coef <- as.vector(drug.data[2:41,6])

# Cisplatin variables
C.biom <- as.vector(drug.data[2:36,7])
C.int <- as.vector(drug.data[1,8])
C.coef <- as.vector(drug.data[2:36,8])

# Gemcitabine variables
G.biom <- as.vector(drug.data[2:36,9])
G.int <- as.vector(drug.data[1,10])
G.coef <- as.vector(drug.data[2:36,10])

# Put these all together into a list
biom <- objects(pattern=".biom")
int <- objects(pattern=".int")
coef <- objects(pattern=".coef")
gene <- objects(pattern=".gene")
drugs <- c("Adriamycin","Cisplatin","Gemcitabine","Methotrexate","Vinblastine")

# Load in the RMA processed expression data
input.rma <- read.table(paste(input.dir,'/',input.rma.c,sep=''),header=TRUE,sep=",", row.names=1)

# Normalize/scale the expression data gene by gene
input.scale <- apply(input.rma,1,scale)
input.scale <- t(input.scale)

# Create empy score matrix to store scores for each drug.
scores.final <- c()

# Loop over all the drugs listed for single drug scores
for (i in 1:length(drugs)){
	# Establish variables to use that were defined above
	drug <- drugs[i]
	cat('Calculating single drug scores for', drug, '...')
	biom.drug <- get(biom[i])
	int.drug <- get(int[i])
	coef.drug <- get(coef[i])
	gene.drug <- get(gene[i])
	
	# Subset biomarkers from dataset
	if (gen.bio==TRUE){
		input.bm <- input.scale[biom.drug,]
		} else  if (uniID==FALSE){
			input.bm <- input.scale[(rownames(input.scale)%in%gene.drug),]
			input.biom <- c()
			for (j in 1:length(gene.drug)){
				rep.gene <- gene.drug[j]
				input.biom <- rbind(input.biom,input.bm[(rownames(input.bm)%in%rep.gene),])
			}
			rownames(input.biom) <- gene.drug
			colnames(input.biom) <- colnames(input.bm)
			input.bm <- input.biom
			coef.drug <- coef.drug[which(gene.drug%in%rownames(input.bm))]
			
		} else input.bm <- input.scale[(rownames(input.scale)%in%gene.drug),]

	# Multiple biomarker coefficients by biomarker expression data	

	mult.coef <- function(x){
		result <- coef.drug * x
	}
	scorematrix <- apply(input.bm, 2, mult.coef)
	
	# Add the results together by columns and add the intercept to get scores for each patient
	scores <- apply(scorematrix,2,sum)
	scores <- as.matrix(scores)
	scores <- int.drug + scores

	# Now convert these to percentile rank scores
	scores <- -1 * scores
	rank.scores <- rank(scores)/length(scores)
	rank.scores <- as.matrix(rank.scores)
	rownames(rank.scores) <- colnames(input.rma)
	colnames(rank.scores) <- drug
	scores.final <- cbind(scores.final,rank.scores)
	cat('done. \n')
}

# Needs to be a data frame for remaining calculations
scores.final <- as.data.frame(scores.final)

###### Combination Score Calculation ######

# Define variables needed for the equations, as outlined in the model

MVAC.int <- 0.761211
MVAC.scale <- 1.091486
MVAC.M <- 0.766369
MVAC.V <- 0.3657074
MVAC.A <- 0.7634502
MVAC.C <- 0.8308788

GC.int <- -0.2456275
GC.scale <- 0.8961938
GC.G <- 2.8852677
GC.C <- 1.6755004

# First calculate the coefficients * scores for the patients. Defining this as beta for the MVAC calculation.
M.beta <- as.matrix(MVAC.M*scores.final$Methotrexate)
V.beta <- as.matrix(MVAC.V*scores.final$Vinblastine)
A.beta <- as.matrix(MVAC.A*scores.final$Adriamycin)
C.beta <- as.matrix(MVAC.C*scores.final$Cisplatin)

rownames(M.beta) <- rownames(scores.final)
rownames(V.beta) <- rownames(scores.final)
rownames(A.beta) <- rownames(scores.final)
rownames(C.beta) <- rownames(scores.final)

# gamma for the GC calculation
G.gamma <- as.matrix(GC.G*scores.final$Gemcitabine)
C.gamma <- as.matrix(GC.C*scores.final$Cisplatin)

rownames(G.gamma) <- rownames(scores.final)
rownames(C.gamma) <- rownames(scores.final)

# Sum for each regimen, I'll name this alpha. Need to add V.beta back in if using.
MVAC.alpha <- as.matrix(rowSums(cbind(M.beta,V.beta,A.beta,C.beta)))
MVAC.alpha <- MVAC.alpha + MVAC.int
GC.alpha <- as.matrix(rowSums(cbind(G.gamma,C.gamma)))
GC.alpha <- GC.alpha + GC.int

# Lambda calculation, e^(-alpha)
MVAC.lambda <- exp(-1*MVAC.alpha)
GC.lambda <- exp(-1*GC.alpha)

# Final probability calculation, as per the model using P = exp(-(5*lambda^scale))
MVAC.mu <- (time*MVAC.lambda)^(MVAC.scale)
GC.mu <- (time*GC.lambda)^(GC.scale)

P.MVAC <- exp(-MVAC.mu)
P.GC <- exp(-GC.mu)

# Save the probability scores
write.csv(P.MVAC,file=paste(save.dir,'/',data,"MVAC.Prob_",time,".YearSurv.csv",sep=''))
write.csv(P.GC,file=paste(save.dir,'/',data,"GC.Prob_",time,".YearSurv.csv",sep=''))

# Save a plot for scores
pdf(file=paste(save.dir,'/',data,'MVAC.vs.GC.Prob_',time,".YearSurv.pdf",sep=''))
# plot MVAC probabilities vs GC probabilities
plot(P.MVAC,P.GC, pch=20, pty="s", ylab="GC Probability Score", xlab="MVAC Probability Score", xlim=c(0,1),ylim=c(0,1), main=paste("COXEN Scores for Combination Treatment"))

#In order to save a file that has sample names instead of points, use the following
#plot(P.MVAC,P.GC, pch=NA, pty="s", ylab="GC Probability Score", xlab="MVAC Probability Score", xlim=c(0,1),ylim=c(0,1), main=paste("COXEN Scores for Combination Treatment"))
#text(P.MVAC,P.GC,rownames(P.GC))

# with margins
abline(0,1)
abline(0.25,1,col=2,lty="dashed")
abline(-0.25,1,col=2,lty="dashed")
abline(0.15,1,col=3,lty="dashed")
abline(-0.15,1,col=3,lty="dashed")

dev.off()

### Repeat plot commands so that it shows up in R window ###

# plot MVAC probabilities vs GC probabilities
plot(P.MVAC,P.GC, pch=20, pty="s", ylab="GC Probability Score", xlab="MVAC Probability Score", xlim=c(0,1),ylim=c(0,1), main=paste("COXEN Scores for Combination Treatment"))
#text(P.MVAC,P.GC,rownames(P.GC))

# with margins
abline(0,1)
abline(0.25,1,col=2,lty="dashed")
abline(-0.25,1,col=2,lty="dashed")
abline(0.15,1,col=3,lty="dashed")
abline(-0.15,1,col=3,lty="dashed")

### Survival Curves ###





