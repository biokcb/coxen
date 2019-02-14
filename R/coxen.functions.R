# Originally written by Paul Williams
#match.by.rownames is a function that is designed to process two gene expression datasets which share some common rownames.
#The output of this function is a list of two matrices with matching rownames, subsets of the input matrices

match.by.rownames <- function( genematrixa, genematrixb ){
	#Do the actual matching
	psmatch <- match( rownames(genematrixa), rownames(genematrixb) )
	if( all(is.na(psmatch)) ) stop('In function match.by.rownames: gene expression matrices share no common row names')
	
	#select subsets
	submatrixa <- genematrixa[!is.na(psmatch),]
	submatrixb <- genematrixb[psmatch[!is.na(psmatch)],]

	if( nrow(submatrixa) != nrow(submatrixb) )
		stop('In function match.by.rownames: submatrices have differing numbers of rows')

	return( list(matcheda=submatrixa,matchedb=submatrixb) )
}




#match.by.file is a function that is designed to process two gene expression datasets which do not share common rownames.
#First the function attempts to load the indicated file, then match the appropriate dataset rownames to the appropriate 
#columns in the file
#File should be set up so that there the first two columns correspond to rownames from either gene expression matrix

match.by.file <- function( genematrixa, genematrixb, filename ){
	#load in the matchfile
	matchfile <- read.table(filename,sep='\t',header=TRUE,as.is=TRUE)

	#Check to see if we can actually match the rownames in the gene matrices to values in the file
	#also I am going to try to set this up so that the order of genematrices a and b do not have to depend on the order of the columns in the file
	tempmatch1 <- match(matchfile[,1], rownames(genematrixa))
	tempmatch2 <- match(matchfile[,2], rownames(genematrixa))
	tempmatch3 <- match(matchfile[,2], rownames(genematrixb))
	tempmatch4 <- match(matchfile[,1], rownames(genematrixb))

	#If either gene matrix has no matching rownames to either column something's wrong and we'll have to quit
	if( all(is.na(tempmatch1)) && all(is.na(tempmatch2)) )
		stop('In function match.by.file: found no matching names for genematrix a in first two columns of file')
	if( all(is.na(tempmatch3)) && all(is.na(tempmatch4)) )
		stop('In function match.by.file: found no matching names for genematrix b in first two columns of file')

	#Check to see which columns have more matches to the gene matrix rownames
	#This will help select which column to match to which gene matrix
	#Also, if one column matches both gene matrices something weird may be going on
	testval1 <- sum(!is.na(tempmatch1)) >= sum(!is.na(tempmatch2))
	testval2 <- sum(!is.na(tempmatch3)) >= sum(!is.na(tempmatch4))
	if( testval1 != testval2 )
		warning(paste('In function match.by.file: one column in file',filename,'seems to match both gene matrices'))

	#select which matches to actually use
	if( testval1 ) amatch <- tempmatch1
	else amatch <- tempmatch2
	if( testval2 ) bmatch <- tempmatch3
	else bmatch <- tempmatch4

	#Get the subsets of the matching files
	nas <- is.na(amatch) | is.na(bmatch)
	submatrixa <- genematrixa[amatch[!nas],]
	submatrixb <- genematrixb[bmatch[!nas],]
	
	if( nrow(submatrixa) != nrow(submatrixb) )
		stop('In function match.by.rownames: submatrices have differing numbers of rows')

	return( list(matcheda=submatrixa,matchedb=submatrixb) )
}




#normalize is a vectorizable function that standardizes a vector x to have mean=0, sd=1
#Can use apply to call this function on the rows or columns of a matrix, but requires transposition to get original dimensions

normalize <- function(x) { (x-mean(x,na.rm=T))/sqrt(var(x,na.rm=T)) }




#cor.test.p is a vectorizable function that returns the p-value from a cor.test call
#Can use apply to call this function on the rows or columns of a matrix

cor.test.p <- function(...) { return( cor.test(...)$p.value ) }




#cor.test.discovery is a function to calculate how correlated gene expression is with drug effectiveness
#genematrix is the matrix of gene expression
#effectiveness is the sensitivity data corresponding to gi50 or sf values
#cormethod is the method of correlation desired and passed on to cor.test.p
#corexact is passed on to cor.test.p as the "exact" variable
#This function returns a matrix of three columns corresponding to
#	[,1] the row number for the given gene
#	[,2] the cor-test p-value for that gene
#	[,3] the FDR q-value calculated from the p-values
#There are three different options for selecting the number of genes returned in the matrix
#	1) if the function is called as cor.test.discovery(x,y,numgenes = n), where n > 0, the function returns the n most significant genes
#	2) if the function is called as cor.test.discovery(x,y,cutoff = q), where 0 < q < 1, the function returns the genes with  FDR q-value <= q
#	3) if neither numgenes or cutoff is specified, all genes are returned
#	numgenes takes precedence over cutoff

cor.test.discovery <- function( genematrix, effectiveness,
	cormethod = 'spearman', corexact = FALSE, numgenes = 0, cutoff = 0.0,... ){
	require(qvalue)
	gene.pval <- apply( genematrix,1,cor.test.p, y=effectiveness, method=cormethod,
					exact=corexact, ...)
	gene.pval[is.na(gene.pval)] <- 1.0
	gene.qval = qvalue(gene.pval)$qvalue

	outmat <- matrix(nrow=length(gene.pval),ncol=3)
	colnames(outmat) <- c('Gene Index','Cor-Test P-Value','Cor-Test FDR Q-Value')
	rownames(outmat) <- rownames(genematrix)
	outmat[,1] <- seq(along=gene.pval)
	outmat[,2] <- gene.pval
	outmat[,3] <- gene.qval
	gene.order <- order(gene.pval)
	if( numgenes > 0 ){
		if( numgenes < nrow( outmat ) ){ return( outmat[gene.order[1:numgenes],,drop=FALSE] ) }
		else return( outmat[gene.order,,drop=FALSE] )
	} else if( cutoff > 0 ){
		return( outmat[gene.order[gene.qval<cutoff],,drop=FALSE] )
	} else {
		return( outmat[,,drop=FALSE] )
	}
} 




#t.test.approximate encapsulates jae's fast but approximate t-test code into a function
#x and y should be the trainingset[,drug.sensitive] and trainingset[,drug.resistant]
#This function returns a matrix of three columns corresponding to
#	[,1] the row number for the given gene
#	[,2] the approximate t-test p-value for that gene
#	[,3] the FDR q-value calculated from the p-values

t.test.approximate <- function( x, y ){
	require(qvalue)
	x.mean = apply(x, 1, mean, na.rm=T)
	y.mean = apply(y, 1, mean, na.rm=T)
	x.var = apply(x, 1, var, na.rm=T);
	y.var = apply(y, 1, var, na.rm=T)
	t.stat = (x.mean-y.mean)/sqrt(x.var/(ncol(x)-1)+y.var/(ncol(y)-1))
	t.pval = pnorm(t.stat)
	t.pval = apply(cbind(t.pval, 1-t.pval),1,min)*2
	t.pval[is.na(t.pval)] <- 1.0
	q.val = qvalue(t.pval)$qvalues
	outmat <- matrix(nrow=length(t.pval),ncol=3)
	outmat[,1] <- seq(along=t.pval)
	outmat[,2] <- t.pval
	outmat[,3] <- q.val
	colnames(outmat) <- c('Gene Index','T-Test P-Value','T-Test FDR Q-Value')
	rownames(outmat) <- rownames(x)
	return(outmat)
}




#t.test.discovery returns a list of significant genes from a t.test, ordered by the FDR Q-value (smallest to largest)
#This function is a wrapper for t.test.approximate
#x and y should be the trainingset[,drug.sensitive] and trainingset[,drug.resistant] and are passed on to t.test.approximate
#This function returns a matrix of three columns corresponding to
#	[,1] the row number for the given gene
#	[,2] the approximate t-test p-value for that gene
#	[,3] the FDR q-value calculated from the p-values
#There are three different options for selecting the number of genes returned in the matrix
#	1) if the function is called as t.test.discovery(x,y,numgenes = n), where n > 0, the function returns the n most significant genes
#	2) if the function is called as t.test.discovery(x,y,cutoff = q), where 0 < q < 1, the function returns the genes with  FDR q-value <= q
#	3) if neither numgenes or cutoff is specified, all genes are returned
#	numgenes takes precedence over cutoff

t.test.discovery <- function( x, y, numgenes=0, cutoff=0.0 ){
	sigmat <- t.test.approximate( x, y )
	gene.order <- order(sigmat[,2])
	if( numgenes > 0 ){
		if( numgenes < nrow( sigmat ) ){ return( sigmat[gene.order[1:numgenes],,drop=FALSE] )}
		else return( sigmat[gene.order,,drop=FALSE] )
	} else if( cutoff > 0 ){
		return( sigmat[gene.order[sigmat[,3]<cutoff],,drop=FALSE] )
	} else {
		return( sigmat )
	}
}




#select.extremes
#function to loop through various numbers of sensitive and resistant cell lines to determine the 'optimal'
#number of sensitive and resistant cell lines
#
#genematrix is the matrix of gene expression
#sensitivity.order is the order statistic representing how sensitive the corresponding cell line is to a drug
#cutoff is an FDR q-value passed on to the t.test.discovery function
#min. and max. num.sensitive and num.resistant define the limits of search-space

select.extremes <- function( genematrix, sensitivity.order, cutoff, 
	min.num.sensitive=8, max.num.sensitive=17, min.num.resistant=9, max.num.resistant=23 ){

	#table for the counts of significant genes with different numbers of sensitive or resistant cell lines
	ngene.matrix = matrix( 0, nrow=max.num.resistant-min.num.resistant+1, ncol=max.num.sensitive-min.num.sensitive+1)
	colnames(ngene.matrix) = min.num.sensitive:max.num.sensitive
	rownames(ngene.matrix) = min.num.resistant:max.num.resistant

	for( n.sensitive in min.num.sensitive:max.num.sensitive ){
		for( n.resistant in min.num.resistant:max.num.resistant ){

			drug.sensitive = sensitivity.order[1:n.sensitive]
			drug.resistant = rev(sensitivity.order)[1:n.resistant]

			#Calculate t-test p-values between sensitive and resistant cell lines
			#this form of the function call returns a matrix with nrow =  the number of genes with FDR < qval.cutoff
			sigmat <- t.test.discovery( x=genematrix[,drug.sensitive], y=genematrix[,drug.resistant],cutoff = cutoff)
			ngene.matrix[n.resistant-min.num.resistant+1,n.sensitive-min.num.sensitive+1] <- nrow(sigmat)
		}
	}

	#Select the optimal numbers of sensitive and resistant cell lines
	max.comb = which(ngene.matrix == max(ngene.matrix), arr.ind=T)
	n.sensitive = max.comb[1,2] + min.num.sensitive - 1
	n.resistant = max.comb[1,1] + min.num.resistant - 1
	drug.sensitive = sensitivity.order[1:n.sensitive]
	drug.resistant = rev(sensitivity.order)[1:n.resistant]
	return( list(drug.sensitive=drug.sensitive, drug.resistant=drug.resistant) )
	#return(ngene.matrix)
}




#COXEN STEP
#Implementation of the correlation of correlations coxen filtration set.  Determines which genes are concordantly regulated
#(relative to other genes in the dataset) between two different datasets.  Returns the indices of the genes which are
#found to be significantly co-regulated.
#
#trainmat is the trainingset gene expression matrix
#coxenmat is the gene expression matrix to compare with the training set matrix
#gene.final is a list of gene indices to be filtered---NOTE!!! the gene.indices MUST CORRESPOND to the rows of the input matrices!
#min.coxen.cor is the minimum allowable coxen correlation coefficient
#max.coxen.cor.pval is the maximum allowable cor-test p-value

coxen.step <- function( trainmat, coxenmat, gene.final,  min.coxen.cor=0.01, max.coxen.cor.pval = 0.05 ){
	cat('Running COXEN filtration step\n')

	n1 <- nrow( trainmat )
	if( nrow(coxenmat) != n1 ){
		stop( 'Error in coxen.step(), trainmat and coxenmat have different numbers of genes!')
	}

	trainmat.cor = cor(t(trainmat))
	coxenmat.cor <- cor(t(coxenmat))

	coxen.cor = matrix(1,nrow=n1,ncol=2)
	for (i in 1:n1) {
   		coxen.cor[i,1] = cor(trainmat.cor[i,-i], coxenmat.cor[i,-i], method="spearman")
		coxen.cor[i,2] = cor.test(trainmat.cor[i,-i], coxenmat.cor[i,-i], method="spearman",exact=FALSE)$p.value
	}
	surviving.genes <- (1:n1)[coxen.cor[,1] > min.coxen.cor & coxen.cor[,2] < max.coxen.cor.pval]
	cat( 'After coxen step,', length(surviving.genes), 'out of', nrow(trainmat), 'genes survive filter\n')
	gene.final = gene.final[surviving.genes ]
	return(gene.final)
}




#LDA-based prediction
#Implementation of singular-value-decomposition of training set, mapping of predictionset onto training set coordinates,
#then linear discriminant analysis of training set followed by prediction on prediction set
#
#trainingset is the gene-expression matrix for the training set
#drug.sensitive and drug.resistant are indices corresponding to the cell lines sensitive or resistant to treatment
#predictionset is the gene-expression matrix for the set for which the predictions will be generated
#best1 is the indices corresponding to the singular values to be used for the predictions

lda.prediction <- function( trainingset, drug.sensitive, drug.resistant, predictionset, best1=1:3){
	require(MASS)

	training.factor = factor(c(rep("Sensitive",length(drug.sensitive)), rep("Resistant",length(drug.resistant))),levels=c("Sensitive","Resistant"))
	#n genes, p chips, svd: M = u d v', n x p  = n x n, n x n, n x p
	svd.out = svd(trainingset)
	U = svd.out$u
	D = svd.out$d
	V = svd.out$v
	V = V[c(drug.sensitive,drug.resistant),]

	# prediction.v = t( 1/d * t(u) * predictionset )
	prediction.v  = t( diag(1/D) %*% t(U) %*% predictionset )

	best1 = best1
	if( is.null(ncol(V)) || ncol(V) < 3 ){
		lda.best1 = lda( V, grouping = training.factor )
	} else {
		lda.best1 = lda( V[,best1], grouping = training.factor )
	}
	lda.best1.out.prediction = predict( lda.best1, prediction.v[,best1] )

	#the lda.best1.out.prediction$x is like a z-score
	model.score = c(lda.best1.out.prediction$x)
	model.score = normalize(model.score)
	coxen.score = 1-pnorm(model.score)

	out <- cbind( lda.best1.out.prediction$posterior[,1], lda.best1.out.prediction$x, coxen.score )
	colnames(out)[1] <- paste('Posterior Probability',colnames(lda.best1.out.prediction$posterior)[1])
	return( out )
		
}
