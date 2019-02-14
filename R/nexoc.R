# function to use the coxen algorithm to predict how well drugs will perform per patient
# created by Kristen Brown
# based on work in collaboration with Jared Fowles, Paul Williams

nexoc <- function(genematrix, coxenmats, indepmat, scale=TRUE, pred.all=TRUE, nsc, drugdata, rev=TRUE, tissue, num=c(12,12), deg=c("t.test", "cor", "SAM", "randomForest"), min.gene=15, max.gene=1000, method=c("mipp", "mod.mipp", "superlearner", "pcr.simple","pcr.step", "randomForest")){
  ####### 1 - DEG selection
  num.sen <- as.numeric(num[1])
  num.res <- as.numeric(num[2])
  if (pred.all){
    final <- matrix(data=NA, nrow=1, ncol=(ncol(coxenmats)+ncol(indepmat)+8))
  } else final <- matrix(data=NA, nrow=1, ncol=(ncol(indepmat)+8))
  if (rev){
    drug.eff <- -1*drugdata[,,drop=FALSE]
  } else drug.eff <- drugdata[,,drop==FALSE]
  
  if (scale){
    # scale input rmas by row/gene
    names1 <- colnames(genematrix)
    genematrix <- apply(genematrix,1,scale)
    genematrix <- t(genematrix)
    colnames(genematrix) <- names1
    
    names2 <- colnames(coxenmats)
    coxenmats <- apply(coxenmats,1,scale)
    coxenmats <- t(coxenmats)
    colnames(coxenmats) <- names2
    
    names3 <- colnames(indepmat)
    indepmat <- apply(indepmat,1,scale)
    indepmat <- t(indepmat)
    colnames(indepmat) <- names3
    
    rm(names1,names2,names3)
    
    cat('calculating scores for ', nsc, '...\n')
    
  } else cat('calculating scores for ', nsc, '...\n')
  
  #remove cell lines missing drug values
  genemat.order <- genematrix[,!is.na(drug.eff)]
  drug.eff <- drug.eff[!is.na(drug.eff)]  
  tissue <- tissue[!is.na(drug.eff),]
  # make sure genematrix is ordered by sensitivity
  genemat.order <- genemat.order[,order(drug.eff)]
  tissue <- tissue[order(drug.eff),]
  drug.eff <- drug.eff[order(drug.eff)]
  # cut nci60 expression matrix into extreme groups
  genemat.col1 <- ncol(genemat.order) - num.res + 1
  genemat.col2 <- ncol(genemat.order)
  genemat.sen <- genemat.order[,c(1:num.sen)]
  genemat.res <- genemat.order[,c(genemat.col1:genemat.col2)]
  
  if (any(deg==c("t.test","cor","SAM"))){
    # t test calculation
    if (any(deg==c("t.test"))){
      gene.tt <- t.test.discovery(genemat.sen,genemat.res)
      gene.tt <- gene.tt[order(gene.tt[,2]),]
      if (gene.tt[1,3]<0.1){
        gene.final <- gene.tt[gene.tt[,3]<0.1,]
        gene.final <- rownames(gene.final)
        final[4]<-length(gene.final)
      } else gene.final <- NULL
      final[4] <- length(gene.final)
    } else if (any(deg==c("SAM"))){
      #sam gene selection 
      sam.final <- degs.select(genemat.order, fdr=0.1, method="hi.lo", num.sen=num.sen, num.res=num.res)
      gene.final <- sam.final[[1]]
      final[6] <- sam.final[[2]]
      if (is.na(gene.final)){
        final[5] <- NA
      } else {
        final[5] <- nrow(gene.final)
        gene.final <- gene.final[,2]
     }
    }
    if (is.null(gene.final)){
      gene.final <- NULL
      cat('method returned 0 differentially expressed biomarkers\n')
    } else {
      cat('method returned',length(gene.final), 'differentially expressed biomarkers\n' )
    }
    if (any(deg==c("cor"))){ 
      # subset if prior testing happened
      if (is.null(gene.final)){
        genemat.cut <- genemat.order
      } else {genemat.cut <- genemat.order[rownames(genemat.order) %in% gene.final,]}
      # correlation test
      gene.ct <- cor.test.discovery(genemat.cut, as.numeric(drug.eff[order(drug.eff)]))
      gene.final <- cbind(gene.tt,gene.ct[gene.tt[,1],])
      gene.final <- gene.final[order(gene.final[,6]),]
      gene.final <- gene.final[gene.final[,5]<0.01,]
      cat('correlation test returned',nrow(gene.final),'genes also highly correlated to drug response\n')   
      if (length(gene.ct)>min.gene){
        gene.final <- rownames(gene.final)[rownames(gene.final) %in% gene.ct]
        final[7] <- length(gene.ct)
        cat('correlation test returned', length(gene.ct), 'biomarkers highly correlated with drug response\n')
      } else {
        gene.final <- rownames(gene.final)
        cat('correlation test returned too few biomarkers before COXEN step, using other biomarkers only if available.... \n')
      }
    }
    if (any(deg==("cor")) && length(gene.final)<min.gene && !is.null(gene.final)){
      cat('Not enough biomarkers for correlation testing, using initial biomarkers only...\n')
    }
  } else if (any(deg==c("randomForest")) && method==c("randomForest")){
    if (is.na(gene.final)){
    names(drug.eff)<- colnames(genemat.order)
    sig <- makesig(drug.eff,genemat.order,tissue,fdr=0.1)  
    gene.final <- sig$signature
    } else {
      names(drug.eff)<- colnames(genemat.order[rownames(genemat.order) %in% gene.final,])
      sig <- makesig(drug.eff,genemat.order[rownames(genemat.order) %in% gene.final,],tissue,fdr=0.1)  
      gene.final <- sig$signature
    }
    cat('randomForest selection leaves ', length(gene.final), 'biomarkers for COXEN step..\n')
  } else {
    cat('deg method not recognized with this model building method...\n')
    gene.final <- NULL
  }
  
  genemat.sr <- cbind(genemat.sen,genemat.res)
  genemat.sr.degs <- genemat.sr[rownames(genemat.sr) %in% gene.final,]
  coxenmat.degs <- coxenmats[rownames(coxenmats) %in% gene.final,]
  indep.degs <- indepmat[rownames(indepmat) %in% gene.final,]
  
  ####### 2 - COXEN Step
  if (length(gene.final)>min.gene){
    gene.coxen <- coxen.step(genemat.sr.degs, coxenmat.degs, c(1:length(gene.final)), min.coxen.cor=0.1, max.coxen.cor.pval = 0.05)
    if (length(gene.coxen)>max.gene && length(gene.coxen)>min.gene){
    gene.coxen <-  gene.coxen[1:max.gene]
      } else gene.coxen <- gene.coxen
    if (length(gene.coxen)>min.gene && length(gene.coxen)<=max.gene){
      gene.final <- gene.final[gene.coxen]
      cat('Left with', length(gene.final), 'genes for model building after coxen filter step\n')
      final[8]<- length(gene.final)
      coxenmat.degs <- coxenmat.degs[rownames(coxenmat.degs) %in% gene.final,]
      indep.degs <- indep.degs[rownames(indep.degs) %in% gene.final,]
      
      if (method=="mipp"){
        cat('Evaluating COXEN biomarkers using MiPP algorithm (LDA)...\n')
        x <- genemat.sr.degs[rownames(genemat.sr.degs) %in% gene.final,]
        y <- factor(c(rep("S",num.sen),rep("R",num.res)))
        
        if (length(gene.final)<100){
          pcut <- 1
        } else if (length(gene.final)<600){
          pcut <- 0.1
        } else pcut <- 0.01
        
        out<-mipp(x=x,y=y,n.fold=5,p.test=1/3,n.split=20,n.split.eval=100,percent.cut=pcut,rule="lda")
        out.model<-out$model
        out.model.eval<-out$model.eval
        colnames(out.model.eval)<-gsub(" ",".",colnames(out.model.eval))
        
        rownames(out.model.eval) <- out.model.eval[,1]
        out.model.eval <- out.model.eval[,2:ncol(out.model.eval)]
        
        #sort table by mean sMiPP score
        out.model.eval<-out.model.eval[order(out.model.eval$mean.sMiPP, decreasing=TRUE),]
        genesonly<-out.model.eval[,1:(ncol(out.model.eval)-9)]
        if (is.matrix(genesonly)){
          for (u in 1:nrow(genesonly)){ #sort the gene sequences numerically
            genesonly[u,]<-genesonly[u,order(genesonly[u,])]
          }
        } else genesonly <- as.matrix(genesonly)
        
        genesonly[is.na(genesonly)]<-0
        genesonly <- unique(genesonly)
        genesonly<-as.matrix(genesonly)
        genesonly <- t(genesonly)  
        
        if (is.matrix(genesonly) && ncol(genesonly)>5){
          qend <- 5
        } else qend <- ncol(genesonly)
        
        models.p<-c()
        probs.sen <- c()
        geness <- c()
        for (q in 1:qend){
          x.train<-as.matrix(x)
          y.train<-(c(rep("Sen",num.sen),rep("Res",num.res)))
          if (pred.all){
            x.test <- cbind(coxenmat.degs, indep.degs)
          } else x.test<- indep.degs
          
          geneseq<-genesonly[,q]
          geneseq<-as.vector(geneseq)
          geneseq<-geneseq[geneseq!=0]
          
          colnames(x.train)<-c(1:ncol(x.train))
          colnames(x.test)<-c(1:ncol(x.test))
          x.train<-t(x.train)
          x.train<-as.matrix(x.train[,geneseq])
          x.test<-t(x.test)
          x.test<-as.matrix(x.test[,geneseq])
          fit<-lda(x.train,y.train)
          predict<-predict(fit,x.test)
          
          probs <- predict$posterior
          prob.sen <- probs[,2]
          probs.sen <- rbind(probs.sen,prob.sen)
          
          geness <- c(geness,length(geneseq))      
        }         
        models.p <- apply(probs.sen,2,mean) 
        final[3] <-  mean(geness)
        final[9:length(final)] <- models.p
        final[1] <- mean(models.p)
        final[2] <- qend
        
      } else if (method=="superlearner"){
        cat('Evaluating COXEN biomarkers using SuperLearner method...\n')
        x <- genemat.order[rownames(genemat.order) %in% gene.final,]
        x <- as.data.frame(t(x))
        y <- drug.eff
        SL.library <- c("SL.bayesglm", "SL.randomForest", "SL.stepAIC", "SL.svm")
        out <- CV.SuperLearner(Y=y,X=x,V=5, SL.library=SL.library, verbose=TRUE)
        models <- predict(out, coxen.mat, X=x, Y=y, onlySL=TRUE)
        
        qend <- 1   
      } else if (method=="pcr.simple"){
        cat('method simple pcr not supported at this time...\n')
      } else if (method=="randomForest"){
        cat('Evaluating COXEN biomarkers using randomForests method...\n')
        #generate biomarker list based on COXEN gene list, based on ivDrug package
        x <- genemat.order[rownames(genemat.order) %in% gene.final,]
       
        common.features <- intersect(row.names(x),row.names(coxenmat.degs))
        trainset <- x[row.names(x) %in% common.features,]
        trainset <- trainset[sort(row.names(trainset)),]
        
        if (pred.all){
          testset <- cbind(coxenmat.degs, indep.degs)
          testset <- testset[row.names(testset) %in% common.features,]
          testset <- testset[sort(row.names(testset)),]
        } else {
          testset <- indep.degs[row.names(indep.degs) %in% common.features,]
          testset <- testset[sort(row.names(testset)),]
        }
        
        if (deg!=c("randomForest")){
          bottomLines <- sig$bottomlines
          topLines <- sig$toplines
        } else {
          bottomLines <- colnames(genemat.res)
          topLines <- colnames(genemat.sen)
        }
        
        bottomgenes <- trainset[row.names(trainset) %in% gene.final,bottomLines]
        topgenes <- trainset[row.names(trainset) %in% gene.final,topLines]
        drugresponse <- c(rep(0,length(bottomLines)),rep(1,length(topLines)))

        train <- cbind(bottomgenes,topgenes)
        test <- testset[row.names(testset) %in% row.names(train),]
        
        rf <- randomForest(t(train),as.factor(drugresponse),ntree=10000)
        pred <- predict(rf,t(test), type="prob")
        
        final[9:length(final)] <- pred[,1]
        final[1] <- mean(pred[,1])
        final[2] <- 1      
      } else if (method=="pcr.mod"){
        cat('method modified pcr not supported at this time...\n')
      } else if (method=="mod.mipp"){
        cat('method modified mipp not supported at this time...\n')
      } else {cat('method not recognized...\n')}
    } else if (length(gene.coxen)<min.gene){
      cat('Too few genes for model building...\n')    
    }
  } else {
    cat('Too few genes for COXEN step....\n')
  }   
  return(final)
}

########## extra required functions 

degs.select<-function(genematrix,fdr=0.1,method="hi.lo", num.sen=12, num.res=12){
  cat('Finding DEGs by SAM...\n')
  Y <- c(rep(2,num.sen), rep(1,num.res))
  samdata<-list(x=genematrix, y=Y, geneid=as.character(1:nrow(genematrix)), genenames=row.names(genematrix), logged2=TRUE)
  log<-capture.output(samr.pac<-samr(samdata, resp.typ="Two class unpaired", nperms=500))
  
  #pick out delta value closest to 0.1
  log2<-capture.output(delta.table <- samr.compute.delta.table(samr.pac))
  colnames(delta.table) <- gsub(" ",".",colnames(delta.table))
  delta.table<-delta.table[complete.cases(delta.table),]
  delta.table<-delta.table[delta.table[,5]!=0,]
  delta.FDR <- which(abs(delta.table[,5] - 0.1)==min(abs(0.1 - delta.table[,5])))	
  if (length(delta.FDR)>1){
    delta.FDR.pick <- which(abs(delta.table[delta.FDR[1],5])==max(abs(delta.table[delta.FDR[2],5])))
    delta.FDR <- delta.FDR[delta.FDR.pick]
  }
  del.max <- delta.table[delta.FDR,1]
  del.min <- delta.table[(delta.FDR-1),1]
  if (delta.table[delta.FDR,5]<0.5){
  cat(del.min,"to ", del.max,'with FDR', delta.table[delta.FDR,5], 'to', delta.table[delta.FDR-1,5],'..\n')
  log2<-capture.output(delta.table <- samr.compute.delta.table(samr.pac,dels=seq(del.min, del.max, 0.002)))
  delta.FDR <- which(abs(delta.table[,5] - 0.1)==min(abs(0.1 - delta.table[,5])))
  if (length(delta.FDR)>1){
    delta.FDR.pick <- which(abs(delta.table[delta.FDR[1],5])==max(abs(delta.table[delta.FDR[2],5])))
    delta.FDR <- delta.FDR[delta.FDR.pick]
  }
  delta.opt <- delta.table[delta.FDR,1]
  
  #DEGs (Differentially Expressed Genes) Summary.
  siggenes.table<-samr.compute.siggenes.table(samr.pac, delta.opt, samdata, delta.table, all.genes=FALSE)
  if (method=="hi.lo"){
    up <- siggenes.table$genes.up
    down <- siggenes.table$genes.lo
    DEGs<-rbind(up,down) #bind high/low table rows
  }else if (method == "lo"){
    DEGs <- siggenes.table$genes.lo 
  }else {
    DEGs<-siggenes.table$genes.up 
  }
  FDR<-delta.table[delta.FDR,5]
  done <- list(DEGs,FDR)
  } else done <- list(NA,NA)
  return(done)
}

#functions from coxen.function file (PDW)

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


gene.coxen <- coxen.step(genemat.sr.degs, coxenmat.degs, c(1:length(gene.final)), min.coxen.cor=0.1, max.coxen.cor.pval = 0.05)


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

