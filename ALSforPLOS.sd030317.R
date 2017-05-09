#'#'#'#'#'#'#'#'#'#'#'#'#'#'#'#'#'#'#'#'#'#'#'#'#'#'#'#
# Select Biomarkers for the ALS motor neuron data set #
#   See the included .pdf for project details         #
#																											#
#                                                     #
# @author Sam Danziger, PhD                           #
# @institution ISB / Seattle Biomed                   #
#'#'#'#'#'#'#'#'#'#'#'#'#'#'#'#'#'#'#'#'#'#'#'#'#'#'#'#
if (require(foreach) && !getDoParRegistered() && require(doParallel) ) { registerDoParallel(cores = parallel::detectCores()-1) }
#if (require(doMC, quietly = T) && !getDoParRegistered() && require(multicore) ) { registerDoMC(cores = multicore:::detectCores()-1) }
library(preprocessCore)
library(verification)
library(ROCR)


#' Call xvalidate.120214 to test the ALS.59 data set
#'
#' @param consensusPreds  Set to TRUE to make multiple prediction for ordered subsets of the features (DEFAULT: FALSE)
#' @export
#' @return  a list containing the mnModel and the predicted HEMI expression
#' @usage resultsList <- xvalidate.032515(consensusPreds=FALSE)
xvalidate.032515 <- function(consensusPreds=FALSE){
	resultsList <- list()

	resultsList[['p1']] <- xvalidate.120214(refName='ALS.59', gsaCut=1, consensusPreds=consensusPreds, useABScut = FALSE, outStr='032515.p1')
	resultsList[['n1']] <- xvalidate.120214(refName='ALS.59', gsaCut=-1, consensusPreds=consensusPreds, useABScut = FALSE, outStr='032515.n1')
	resultsList[['a1']] <- xvalidate.120214(refName='ALS.59', gsaCut=1, consensusPreds=consensusPreds, useABScut = TRUE, outStr='032515.a1')

	#Tally the number of features for each gsa cut
	featureCount <- list()
	for (cutName in c('1', '-1', '1.useABScut')) {
		pattern <- paste( 'classifyDataSets.ALS.59', '[0-9]+', cutName, 'RData', sep=".")
		fNames <- list.files(path='./intermediate', pattern=pattern, full.names=TRUE)
		
		m1 <- c()
		m1.mProbes <- c()
		m2 <- c()
		m2.mProbes <- c()
		for (fName in fNames) {
			classifyDataSets <- get(load(fName)[1])
			m1 <- c(m1, ncol(classifyDataSets[grepl('^method1.train.ALS.59', names(classifyDataSets))][[1]]))
			m1.mProbes <- c(m1.mProbes, ncol(classifyDataSets[grepl('^method1.mProbes.train.ALS.59', names(classifyDataSets))][[1]]))
			m2 <- c(m2, ncol(classifyDataSets[grepl('^method2.train.ALS.59', names(classifyDataSets))][[1]]))
			m2.mProbes <- c(m2.mProbes, ncol(classifyDataSets[grepl('^method2.mProbes.train.ALS.59', names(classifyDataSets))][[1]]))
			featureCount[[cutName]][[fName]]
		}
		featureCount[[cutName]] <- list(m1=m1,m1.mProbes=m1.mProbes,m2=m2,m2=m2.mProbes)
	}

	invisible(resultsList)
}

#' Make predictions for 'ALS.new.122012' and 'ALS.new.090513' from the 'ALS' data set
#'
#' @param gsaCut  The GSA pathway cutoff score (DEFAULT: 1.0)
#' @param gsaPval  Set to cut GSA based on the pvalue.lo and the pvalue.hi.  NULL means no cutoff (DEFAULT: 0.05)
#' @param firmaCut  The FIRMA cutoff score (DEFAULT: 1.5)
#' @param cutLowExp  The low expression cutoff (DEFAULT: 200)
#' @param runXval  Set to TRUE to xvalidate training set (DEFAULT: FALSE)
#' @param numMprobes  Set to NULL for no mProbes, inf for all (DEFAULT: 600)
#' @param standardize  Set to TRUE to standardize data (Z-score) before classification (DEFAULT: TRUE)
#' @param rmSDpaths  Set to TRUE to remove the sd pathways. sd pathways were made before test sets (DEFAULT: FALSE)
#' @param getClusterStack  Set to TRUE to use cMonkey results instead of gene sets (DEFAULT: FALSE)
#' @param consensusPreds  Set to TRUE to make multiple prediction for ordered subsets of the features (DEFAULT: FALSE)
#' @param useABSgsa  Set to TRUE to use abs(GSA) > gsaCut rather than GSA > gsaCut (DEFAULT: FALSE)
#' @param outStr  The string to use for the output files (DEFAULT: '120214')
#' @export
#' @return  a list containing the mnModel and the predicted HEMI expression
#' @usage resultsList <- xvalidate.120214(refName='ALS', gsaCut=1, gsaPval=NULL, firmaCut=1.5, cutLowExp=200, runXval = FALSE, numMprobes=600, standardize=FALSE, rmSDpaths=FALSE, getClusterStack=FALSE, consensusPreds=FALSE, useABScut = FALSE, outStr='120214')
xvalidate.120214 <- function(refName='ALS', gsaCut=1, gsaPval=NULL, firmaCut=1.5, cutLowExp=200, runXval = FALSE, numMprobes=600, standardize=FALSE, rmSDpaths=FALSE, getClusterStack=FALSE, consensusPreds=FALSE, useABScut = FALSE, outStr='120214'){
	#Load the gene and pathway expression scores for the three data sets
	topGSA <- findTopGSA(runNames=c(refName), forceRebuild=FALSE, pVal.cut=gsaPval, gsaCut=1, absGSAcut=TRUE, stripCR = TRUE, getClusterStack=FALSE, xval=TRUE)
	
	#makeClassifierDataSets should be modified so that if it detects a cross-validated dataset
	#  it will make cross-validated datasets that include a train and test subset.
	nFolds <- sum(grepl('ALS', colnames(topGSA)))
	resultLists <- list() #The NB.SVM.EV.RF prediction
	resultLists.svm <- list() #The SVM only prediction
		
	datasetLists <- list()
	for (foldNum in 1:nFolds) {
		fName <- paste('./intermediate/classifyDataSets', refName, foldNum, gsaCut, 'RData', sep=".")
		if(useABScut == TRUE) { fName <- sub('\\.RData$', '.useABScut.RData', fName) }
		if (!file.exists(fName)) {
			classifyDataSets <- makeClassifierDataSets(refName=refName, topGSA = topGSA, cutScore=gsaCut, cutLowExp=cutLowExp, rmSDpaths=rmSDpaths, writeOut=FALSE, includeMprobes = TRUE, xValIdx=foldNum, useABScut = useABScut)
			names(classifyDataSets) <- sub( '\\.[-]+1$', '', names(classifyDataSets))
			save(classifyDataSets, file=fName)
		} else {
			classifyDataSets <- get(load(fName)[1])
		}
		datasetLists[[foldNum]] <- classifyDataSets
		for (meth in c('method1', 'method2')) {
		#for (meth in c('method2')) {
        		for (useMprobes in c(FALSE, TRUE) ) {
        		#for (useMprobes in c(TRUE) ) {
        			if (useMprobes == TRUE) {meth <- paste(meth, 'mProbes', sep=".")}
				trainName <- paste(meth, 'train', refName, foldNum+7, sep='.') #+7.  Hackish
				testName <- paste(meth, 'test', refName, foldNum+7, sep='.') #+7.  Hackish
				#trainLists <- classifyDataSets[trainName]
				trainLists <- classifyDataSets[grepl(paste('^', trainName, sep=""), names(classifyDataSets))]
				#testLists <- classifyDataSets[testName]
				testLists <- classifyDataSets[grepl(paste('^', testName, sep=""), names(classifyDataSets))]

				if (consensusPreds == TRUE) {
					numsFeatures <- 1:(ncol(trainLists[[1]])-1) #Better preds, slower
				} else {
					numsFeatures <- ncol(trainLists[[1]])-1
				}
				predResults <- makePreds(trainData=trainLists[[1]], testData=testLists[[1]], 
					numsFeatures=numsFeatures, onlySVM=FALSE)
				if (is.null(resultLists[[meth]])) { 
					resultLists[[meth]] <- predResults$predTally$All
					resultLists.svm[[meth]] <- predResults$predTally$SVM
				} else {
					resultLists[[meth]] <- rbind(resultLists[[meth]], predResults$predTally$All)
					resultLists.svm[[meth]] <- rbind(resultLists.svm[[meth]], predResults$predTally$SVM)
				}
        		} #for (useMprobes in c(FALSE, TRUE) ) {
        	} #for (meth in c('method1', 'method2')) {
	}
	fName <- paste("./intermediate/xvalResults.120214", refName, Sys.Date(), 'RData', sep='.')
	if (consensusPreds == TRUE) { fName <- sub('RData$', 'consensus.RData', fName) }
	if (useABScut == TRUE) { fName <- sub('\\.RData$', '.useABScut.RData', fName) }
	save(resultLists, resultLists.svm, file=fName)
	
	aggXval <- NULL
	for (meth in c('method1', 'method2')) {
        	for (useMprobes in c(FALSE, TRUE) ) {
        		if (useMprobes == TRUE) {meth <- paste(meth, 'mProbes', sep=".")}
        		accRow <- calcAccs(resultLists[[meth]]$prediction, resultLists[[meth]]$actual)
        		rownames(accRow) <- meth
        		if (is.null(aggXval)) { aggXval <- accRow } else { aggXval <- rbind(aggXval, accRow)}
        		
        		accRow <- calcAccs(resultLists.svm[[meth]]$prediction, resultLists.svm[[meth]]$actual)
			rownames(accRow) <- paste(meth, 'svm', sep=".")
			if (is.null(aggXval)) { aggXval <- accRow } else { aggXval <- rbind(aggXval, accRow)}
		}
	}
	fName2 <- paste("./intermediate/xvalResults", outStr, refName, Sys.Date(), 'csv', sep='.')
	if (consensusPreds == TRUE) { fName2 <- sub('csv$', 'consensus.csv', fName2) }
	write.csv(aggXval, file=fName2)
	
	rv <- list(resultLists=resultLists, resultLists.svm=resultLists.svm, fName=fName,
		aggXval=aggXval, fName2=fName2)
	return(rv)
}


#' Rank features using mutual information
#'   Cross validate to pick consistently best features
#'
#' @param dataMat  The attribute matrix
#' @param classIdx  The index for the class variable (DEAFULT: nrow(dataMat))
#' @param folds  The number of folds to break up the data (DEFAULT: 10)
#' @param verbose  Set to TRUE to show output (DEFAULT: FALSE)
#' @export
#' @usage featureList <- xvalMIranker(dataMat, classIdx=nrow(dataMat), folds=10, verbose=FALSE)
xvalMIranker <- function(dataMat, classIdx=ncol(dataMat), folds=10, verbose=FALSE) {
	cuts <- cut(1:nrow(dataMat), folds)
	randIdxs <- order(runif(nrow(dataMat)))
	
	featureMatrix <- NULL
	for (curCut in unique(cuts)) {
		if (verbose == TRUE) { cat(curCut, ", ") }
		idxs <- randIdxs[cuts!=curCut]
		curTest <- dataMat[idxs,]
		#featureList <- order(useMutualInformation(curTest), decreasing=FALSE)
		featureList <- useMutualInformation(curTest)
		names(featureList) <- colnames(curTest)
		featureMatrix <- cbind(featureMatrix, featureList)
		#rownames(featureMatrix) <- colnames(curTest)
	}
	colnames(featureMatrix) <- unique(cuts)

	rankMatrix <- apply(featureMatrix, 2, function(x) {rank(x, na.last = FALSE, ties.method='min')})
	rankMatrix <- nrow(featureMatrix) - rankMatrix

	featureMeans <- rowMeans(rankMatrix)
	featureSds <- apply(rankMatrix, 1, sd)
	rv <- data.frame(rankMean = featureMeans, rankSds = featureSds)
	rv <- rv[order(rv$rankMean), ]
	
	return(rv)
}

#' Rank features using mutual information
#'   This only seems to be the method in [R] that's reasonably quick
#'   The maximum score has the highest importance
#'
#' @param dataMat  The attribute matrix
#' @param classIdx  The index for the class variable (DEAFULT: nrow(dataMat))
#' @export
#' @usage featureList <- useMutualInformation(dataMat, classIdx=ncol(dataMat))
useMutualInformation <- function(dataMat, classIdx=ncol(dataMat)) {
	library(entropy)
	
	classVector <- as.factor(dataMat[,classIdx])
	miVect <- rep(NA, ncol(dataMat))
	names(miVect) <- colnames(dataMat)
	for ( i in 1:ncol(dataMat)) {
		featureName <- colnames(dataMat)[i]
		if ( i != classIdx ) {
			curVector <- as.factor(dataMat[,i])
			miVect[[featureName]] <- mi.plugin(rbind(curVector, classVector))
		}
	}
	#miVect <- sort(miVect)
	return(miVect)
}

#' Our mouse motor neurons are derived from mouse stem cells, and it's not clear
#'   what percentage of cells are transform.  This function uses mnGenes as proxy
#'   markers for the fraction of motor neuron cells and scales accordingly.
#'   mnGenes <- c(Olig2, Mnx1, Isl1, Lhx3)
#'
#' This function will use NCAR data and mnGenes to make a linear model of the 
#'   expected expression of the target gene.
#'
#' @param targetGene  The target gene to make the model of (DEFAULT: 'Slc7a1')
#'                    NOTE: gene names with dashes will be changed to periods by make.names
#' @param NCARdata  The gene expression for the controls to build the model.  Make sure HEMI and NCAR in the column name (DEFAULT: NULL, i.e. load ALS data)
#' @param mnGenes  The motor neuron marker genes (DEFAULT: c('Olig2', 'Mnx1', 'Isl1', 'Lhx3'))
#' @param plotIt  Set to TRUE to plot the genes (DEFAULT: FALSE)
#' @export
#' @return  a list containing the mnModel and the predicted HEMI expression
#' @usage mnModel <- makeMNmodel(targetGene='Slc7a1', NCARdata=NULL, mnGenes=c('Olig2', 'Mnx1', 'Isl1', 'Lhx3'), plotIt=FALSE)
makeMNmodel <- function(targetGene='Slc7a1', expData=NULL, mnGenes=c('Olig2', 'Mnx1', 'Isl1', 'Lhx3'), plotIt=FALSE) {
	targetGene <- make.names(targetGene)
	mnGenes <- mnGenes[!(mnGenes %in% targetGene)]
	
	if(is.null(expData)) { 
		expData <- get(load('./intermediate/geneEx.ALS.RData')[1])
	}
	
	if (plotIt == TRUE) { #Used for debugging
		ncarExp <- expData[mnGenes , grepl('_NCAR', colnames(expData))]
		
		hemiExp <- expData[mnGenes , grepl('_HEMI', colnames(expData))]
		
		cols <- rainbow(length(mnGenes))
		names(cols) <- mnGenes
		
		ylim <- c(.9*min(hemiExp, ncarExp), 1.1*max(hemiExp, ncarExp))	
		
		for (curGene in mnGenes) {
			if (curGene == mnGenes[1]) {
				plot(hemiExp[curGene,], ylim=ylim, ylab="Gene Expression", xlab="Experiment #", col=cols[curGene], main="Motor Neuron Gene Expression", type='o')
			} else {
				lines(hemiExp[curGene,], ylim=ylim, col=cols[curGene], type='o')
			}
			
			lines(ncarExp[curGene,], ylim=ylim, col=cols[curGene], type='o', lty=3)
			
		}
		legend('topright', legend=mnGenes, col=cols, lty=1)
		legend('topleft', legend=c('HEMI', 'NCAR'), lty=c(1,3))
	}
	
	ncarExp <- expData[, grepl('_NCAR', colnames(expData))]
	hemiExp <- expData[, grepl('_HEMI', colnames(expData))]		

	#Note: Nkx2-5 causes trouble in the formula because the - is interpreted as a minus.
	#  Either I can change the annotation upstream, or I can just modify the equation.
	#  I will chose to modify the equation

	expForm <- formula(paste(targetGene, ' ~ ', paste(rownames(ncarExp[mnGenes,]), collapse=' + '), sep=""))
	curData <- as.data.frame(t(ncarExp[c(targetGene, mnGenes),]))
	mnModel <- lm(expForm, curData)
	curData <- as.data.frame(t(hemiExp[c(targetGene, mnGenes),]))
	predHemi <- predict(mnModel, newdata=curData)
	
	return(list(mnModel=mnModel, predHemi=predHemi))
}




#' Modified from print4Weka, but restructured to base filters on first dataset and then match other datasets to that.
#'   From the gene expression data
#'   Method 1: Use genes in initial ALS data set that have a GSA score > a cutoff
#'   Method 2: Use aggregate annotation expression in initial ALS data set that have a GSA score > a cutoff
#'   mnGenes <- c(Olig2, Mnx1, Isl1, Lhx3)
#'   12-03-14: Modified to make it possible to use this for n-fold X-validation.  
#'             This also required modifications to 'findTopGSA', to return cross-validated data with names such as ALS.8.22exp.
#'             ALS.8.22exp refers to the first experiment being excluded (i.e. 8-7 = 1).  ALS.9.22exp 
#'
#' @param refName  The reference dataset name.   (DEFAULT: 'ALS')
#' @param topGSA  The top GSA scores (DEFAULT: findTopGSA())
#' @param qNorm  Set to TRUE to se mean to 0 and SD to 1 (DEFAULT: TRUE)
#' @param cutScore  The ALS score to cut at.  If negative, then will use <= unless useABScut=TRUE (DEFAULT: 1)
#' @param fakeRatios  Set to TRUE to pick a NCAR reference (DEFAULT: FALSE)
#' @param useMNmodel  Set to TRUE to build a linear model for gene expression based on motor neuron marker genes.  Don't use with fakeRatios (DEFAULT: TRUE)
#' @param cutLowExp  Any expression level must be above this threshold to be included.  100 is probably a good number (DEFAULT: NULL)
#'                   NOTE: Potential x-val bug - it uses *ALL* data, not just ALS
#' @param forceRebuild  TRUE to rebuild intermediate files.  Not implemented 1/3/14 (DEFAULT: FALSE)
#' @param rmSDpaths  Set to TRUE to remove the SD pathways (DEFAULT: FALSE)
#' @param writeOut  Set to TRUE to write output files (DEFAULT: TRUE)
#' @param includeMprobes Include datasets that filters out features with mProbes FDR < 1  (DEFAULT: TRUE)
#' @param xValIdx  An index into topGSA and the data set to indicate which topGSA column to use and which experiment to exclude (DEFAULT: NULL)
#' @param useABScut  Set to true to use |GSA| > cutScore rather than GSA > cutScore (DEFAULT: FALSE)
#' @param sortByRefGSA  Set to true to force features to be sorted by GSA scores in reference, not overall score.  Should be TRUE for blind predictions. 04-10-15 (DEFAULT: FALSE)
#' @param ntree  The number of tree in the randomForest mProbes selection. 04-24-15 (DEFAULT: 1000)
#' @export
#' @usage classifyDataSets <- makeClassifierDataSets(refName='ALS', topGSA=NULL, qNorm=FALSE, fakeRatios=FALSE, cutScore=1, useMNmodel=TRUE, cutLowExp=NULL, rmSDpaths=FALSE, writeOut=TRUE, includeMprobes=TRUE, xValIdx=NULL, useABScut=FALSE, sortByRefGSA=FALSE, ntree=1000)
makeClassifierDataSets <- function(refName='ALS', topGSA=NULL, qNorm=FALSE, fakeRatios=FALSE, cutScore=1, useMNmodel=TRUE, cutLowExp=NULL, rmSDpaths=FALSE, writeOut=TRUE, includeMprobes=TRUE, xValIdx=NULL, useABScut=FALSE, sortByRefGSA=FALSE, ntree=1000) {
	#Method 1: Use genes in initial ALS data set that have a GSA score > cutScore
	#Method 2: Use aggregate annotation expression in initial ALS data set that have a GSA score > cutScore
	
	#TODO: Return features for other data sets selected on the first data set.  THIS IS NECESSARY!!!
	if (is.null(topGSA)) {
		if (is.null(xValIdx)) {
			topGSA <- findTopGSA()  #This will need to be updated later to return new data set
			#Note 01-12-14:  When I ran the real AL experiment, did it do this, or did I run findTopGSA separately...
		} else {
			topGSA <- findTopGSA(runNames = refName, xval = TRUE)  #This will need to be updated later to return new data set
		} #if (is.null(xValIdx)) {
	}
	
	if (useMNmodel == TRUE) { fakeRatios <- FALSE }
	mnModels <- NULL #The MN model will be calculated on the first data set, but not on any other
	refCol <- NULL #A reference column will be calculated only on the first data set.
	mpNames.m1 <- mpNames.m2 <- NULL #Use mp names only from reference
	
	if (rmSDpaths==TRUE) {
		sdIdx <- grepl('^sd\\.', rownames(topGSA))
		topGSA <- topGSA[!sdIdx,]
	}
	
	expNames <- sub('\\.[0-9]+exp$', '', colnames(topGSA)[grepl('ALS' , colnames(topGSA))])
	colnames(topGSA)[grepl('ALS' , colnames(topGSA))] <- expNames
	
	#If the is a cross-validated run, select the appropriate topGSA column and data set
	if (!is.null(xValIdx)) {
		expName <- expNames[xValIdx]
		expNames <- strsplit(expName, '\\.')[[1]][1]
		idx <- as.numeric(strsplit(expName, '\\.')[[1]][2]) - 7
		refName.GSA <- expName
		#Decision time.  Should this use refName and the xVal flag to break this up into a single data set.
		#  If that's the case, then I will have to call this multiple times. 
		#  I think that is the best solution.  Strip the .8 or whatever at the end of the *first* data column
		#    then parse it to create the relevant training and test data sets.
	} else {
		refName.GSA <- refName
	}
	#01-08-14 - for blind predictions I used GSA > 1 rather than |GSA| > 1.  
	#  Modify the code accordingly to redo the xvalidation and active learning.
	topGSA <- topGSA[order(abs(topGSA[,which(colnames(topGSA)==refName.GSA)[1]]),decreasing=TRUE),]  #Added on 04-10-15.  Should have done this for blind predictions for consensus to work properly
	if (useABScut == FALSE) {
		if (cutScore >= 0) {
			cutTopGSA <- topGSA[which(topGSA[,which(colnames(topGSA)==refName.GSA)[1]] >= cutScore),] 
		} else { #negative
			cutTopGSA <- topGSA[which(topGSA[,which(colnames(topGSA)==refName.GSA)[1]] <= cutScore),] 
		}
	} else {
		cutTopGSA <- topGSA[which(abs(topGSA[,which(colnames(topGSA)==refName.GSA)[1]]) >= abs(cutScore)),] 	
	}
	topGenes <- unique(sub(' ', '', unlist(strsplit(as.character(cutTopGSA$genes), ';'))))
	topGenes <- make.names(topGenes[!is.na(topGenes)])

	#calcuate which genes never have high expression
	if (!is.null(cutLowExp)) {
		cutLowBool <- NULL
		for (expName in expNames) {
			#cat(expName, '\n')
			geneExFile <- paste('./intermediate/geneEx.', refName, '.', expName, '.RData', sep="")
			if (file.exists(geneExFile)) {
				geneEx <- get(load(geneExFile)[1])
			} else {
				if(is.null(xValIdx)) {
					expData <- getFIRMAandGENE(expName)
				} else {
					expData <- getFIRMAandGENE(refName)
				}
				geneEx <- exonToGene(origData=expData$geneExpr.vals)
				save(geneEx, file=geneExFile)
			}
			if (is.null(cutLowBool)) {
				cutLowBool <- rep(TRUE, nrow(geneEx))
				names(cutLowBool) <- rownames(geneEx)
			}
			geneEx <- geneEx[rownames(geneEx) %in% names(cutLowBool),]
			cutLowBool <- cutLowBool & apply(geneEx, 1, function(x) {all(x < cutLowExp)})
		} #for (expName in expNames) {
	} #if (!is.null(cutLowExp)) {

	#Break up the training set if this is supposed to be a cross-validated run
	if (!is.null(xValIdx)) {
		inGeneExFile <- paste('./intermediate/geneEx.', refName, '.', expNames[1], '.RData', sep="") #expNames should only have 1 element
		#This block of code should probably be its own function. Mind the data set loaded
		if (file.exists(inGeneExFile)) {
			geneEx <- get(load(inGeneExFile)[1])
		} else {
			expData <- getFIRMAandGENE(refName)
			geneEx <- exonToGene(origData=expData$geneExpr.vals)
			save(geneEx, file=inGeneExFile)
		}
		
		testName <- colnames(geneEx)[xValIdx]
		geneEx.train <- geneEx[,colnames(geneEx) != testName, drop=FALSE]
		geneEx.test <- geneEx[,testName, drop=FALSE]
		
		trainNamePat <- paste('train', refName.GSA, sep='.')
		testNamePat <- paste('test', refName.GSA, sep='.')
		expNames <- c(trainNamePat, testNamePat)

		trainNameFile <- paste('./intermediate/geneEx.', refName, '.', trainNamePat, '.RData', sep="")
		save(geneEx.train, file=trainNameFile)
		testNameFile <- paste('./intermediate/geneEx.', refName, '.', testNamePat, '.RData', sep="")
		save(geneEx.test, file=testNameFile)
	} else {
		#Reorder expNames so refName comes first
		expNames <- c(expNames[expNames==refName], expNames[expNames!=refName])
	} #if (!is.null(xValIdx)) {
	
	outDFs <- list()
	for (expName in expNames) {
		cat(expName, '\n')
		geneExFile <- paste('./intermediate/geneEx.', refName, '.', expName, '.RData', sep="")
		if (file.exists(geneExFile)) {
			geneEx <- get(load(geneExFile)[1])
		} else {
			if(is.null(xValIdx)) {
				expData <- getFIRMAandGENE(expName)
			} else {
				expData <- getFIRMAandGENE(refName)
			}
			geneEx <- exonToGene(origData=expData$geneExpr.vals)
			save(geneEx, file=geneExFile)
		}
		
		if (is.null(cutLowExp)) {
			cutLowBool <- rep(FALSE, nrow(geneEx))
			names(cutLowBool) <- rownames(geneEx)
		}
		
		#Remove the low expression genes from consideration
		geneEx <- geneEx[!cutLowBool,,drop=FALSE]
		
		#Select set first NCAR as a reference column
		#browser()
		if ( fakeRatios == TRUE ) {
			#Only use the first dataset to calculate the reference column (refCol)
			if (is.null(refCol)) {
				if (any(colnames(geneEx) == "20110127_motorneuron_NCAR_24hr_rep4")) {
					refColName <- "20110127_motorneuron_NCAR_24hr_rep4"
				} else {
					refColName <- colnames(geneEx)[ grep('NCAR', colnames(geneEx))[1] ]
				}
				refCol <- geneEx[, refColName, drop=FALSE]
			} #if (is.null(refCol)) {
			geneEx <- geneEx[, colnames(geneEx) != refColName, drop=FALSE]
			curRowNames <- rownames(geneEx)
			geneEx <- log2(apply(geneEx, 2, function(x) {x/refCol} ))
			rownames(geneEx) <- curRowNames
		}
		
		#Quantile normalize?  Convert to Z scores
		if (qNorm==TRUE) {
			geneEx <- apply(geneEx, 2, function(x) { (x-mean(x)) / sd(x)})
		}
		
		#NOTE! Only  35/82 match.  Perhaps this is due to MOUSE / HUMAN mismatch
		geneEx.all <- geneEx
		geneEx <- geneEx.all[toupper(rownames(geneEx)) %in% topGenes,,drop=FALSE]
		
		if (useMNmodel == TRUE) {
			#First generate models
			if (is.null(mnModels)) {
				#Note, it's very easy to run out of memory here
				#mnModels <- foreach (curGene = rownames(geneEx)) %dopar% {
				mnModels <- foreach (curGene = rownames(geneEx)) %do% {
					makeMNmodel(targetGene = curGene, expData = geneEx.all)
				}
				names(mnModels) <- rownames(geneEx)
			} #if (is.null(mnModels)) {
			
			#for (curGene in rownames(geneEx)) {
			geneEx.new <- foreach (curGene = rownames(geneEx), .combine=rbind) %dopar% {
				
				#Only run on the first dataset, for all others, reuse the same model
				mnModel <- mnModels[[curGene]]
				
				#Use mnModel to make some ratios
				refGenes <- names(mnModel$mnModel$coefficients)
				refGenes <- refGenes[refGenes != "(Intercept)"]
				predExp <- predict(mnModel$mnModel, newdata=as.data.frame(t(geneEx.all[refGenes,])))
				#Correct for values falling below the noise level
				#predExp[predExp < min(geneEx[curGene,])] <- min(geneEx[curGene,])
				if (any(predExp <= 0)) {  #Correct for negative predictions when gene expression is very low
					bias <- 1.1 * abs(min(predExp[predExp != 0]))
				} else {
					bias <- 0
				}
				outrow <- log2((geneEx.all[curGene,]+bias) / (predExp+bias)) #log2(geneEx[curGene,] / predExp)
				outrow #geneEx[curGene,]
			}
			rownames(geneEx.new) <- rownames(geneEx)
			colnames(geneEx.new) <- colnames(geneEx)
			geneEx <- geneEx.new
		}
		
		classNames <- rep('?', ncol(geneEx)) #makes vector for class labels
		names(classNames) <- colnames(geneEx) 
		classNames[grepl('_NCAR', colnames(geneEx))] <- 'NCAR'
		classNames[grepl('_HEMI', colnames(geneEx))] <- 'HEMI'
		outDF <- cbind(t(geneEx), data.frame(class=classNames))
		outFile <- paste('./intermediate/weka.Method1', refName, expName, 'ref', expNames[1], 'cut', cutScore, 'qNorm', qNorm, 'fakeRatios', fakeRatios, 'useMNmodel', useMNmodel, 'cutLowExp', cutLowExp , 'csv', sep=".")
		if (useABScut == TRUE) {
			outFile <- sub('\\.csv$', '.useABScut.csv', outFile)
		}
		if(writeOut == TRUE) { write.csv(outDF, outFile, row.names=FALSE) }
		outDFs[[paste('method1', expName, cutScore, sep=".")]] <- outDF
		
		#mProbes
		if (includeMprobes == TRUE) {
			if (is.null(mpNames.m1)) {  #Use the first reference set
				#mp <- mProbes(outDF)
				mp <- mProbes(outDF, type='randomForest', ntree=ntree) #4/15/15.  MI runs out of memory
				mpNames.m1 <- c(names(mp)[mp < 1], colnames(outDF)[ncol(outDF)])
			}
			outDFs[[paste('method1.mProbes', expName, cutScore, sep=".")]] <- outDF[, mpNames.m1, drop=FALSE]
		}
		
		outDF <- matrix(0, nrow=nrow(cutTopGSA), ncol=ncol(geneEx), dimnames=list(rownames(cutTopGSA), colnames(geneEx)))
		for (curAnnot in rownames(outDF)) {
			genes <- sub( ' ', '', unlist(strsplit(as.character(cutTopGSA[curAnnot, 'genes']), ';')))
			curCols <- colMeans(geneEx[toupper(rownames(geneEx)) %in% genes, , drop=FALSE])
			outDF[curAnnot, names(curCols)] <- curCols
		}
		
		#Filter out repeated annotations
		outDF <- outDF[!duplicated(cutTopGSA$'genes'),,drop=FALSE]
		rep('HEMI', ncol(outDF))
		classNames[grepl('_NCAR', colnames(outDF))] <- 'NCAR'
		outDF <- cbind(t(outDF), data.frame(class=classNames))

		#Remove any column that is all NAs, caused by too mant genes being removed to low signal?
		outDF <- outDF[,!apply(outDF, 2, function(x) { all(is.na(x)) }), drop=FALSE]
		
		outFile <- paste('./intermediate/weka.Method2', refName, expName, 'ref', expNames[1], 'cut', cutScore, 'qNorm', qNorm, 'fakeRatios', fakeRatios, 'useMNmodel', useMNmodel, 'cutLowExp', cutLowExp, 'csv', sep=".")
		if (useABScut == TRUE) {
			outFile <- sub('\\.csv$', '.useABScut.csv', outFile)
		}
		if(writeOut == TRUE) { write.csv(outDF, outFile, row.names=FALSE) }
		outDFs[[paste('method2', expName, cutScore, sep=".")]] <- outDF

		#mProbes
		if (includeMprobes == TRUE) {
			if (is.null(mpNames.m2)) {
				mp <- mProbes(outDF)
				mpNames.m2 <- c(names(mp)[mp <1], colnames(outDF)[ncol(outDF)])
			}
			outDFs[[paste('method2.mProbes', expName, cutScore, sep=".")]] <- outDF[, mpNames.m2, drop=FALSE]
		}
	} #for (expName in expNames) {
	
	return(outDFs)	
}

#' Compare the gene enrichments to find annotations that are enriched across all runs
#'   By default, use everything with a GSA score > 1?
#'   To Do 1/3/14 - Include flag to force rebuild of intermediate files.
#'   Use annotationEnrich to fix this
#'
#' @param runNames  The data set names (DEFAULT: c('ALS', 'ALS.new.122012', 'ALS.new.090513'))
#' @param forceRebuild  Set to TRUE to rebuild (DEFAULT: FALSE)
#' @param pVal.cut  Set to cut based on the pvalue.lo and the pvalue.hi (DEFAULT: NULL)
#' @param gsaCut  DEPRICATED AND NON-FUNCTIONAL Remove pathways with gsa scores < gsaCut (DEFAULT: 1)
#' @param absGSAcut  Set to TRUE to use the absolute value - was FALSE pre070814 (DEFAULT: TRUE)
#' @param stripCR  Set to TRUE to strip Circadian and Riluzole (DEFAULT: TRUE)
#' @param getClusterStack  Set to TRUE to use the biclusters instead (DEFAULT: FALSE)
#' @param xval  Set to true to build leave one out cross-validated sets  (DEFAULT: FALSE)
#' @param nperms  The number of permutations  (DEFAULT: 300)
#' @export
#' @usage topGSA <- findTopGSA(runNames=c('ALS', 'ALS.new.122012', 'ALS.new.090513'), forceRebuild=FALSE, pVal.cut=NULL, gsaCut=1, absGSAcut=TRUE, stripCR = TRUE, getClusterStack=FALSE, xval=FALSE, nperms=300)
findTopGSA <- function(runNames=c('ALS', 'ALS.new.122012', 'ALS.new.090513'), forceRebuild=FALSE, pVal.cut=NULL, gsaCut=1, absGSAcut=TRUE, stripCR = TRUE, getClusterStack=FALSE, xval=FALSE, nperms=300) {
	cat('NOTE: absGSAcut only effects scoring function.\n Use makeClassifierDataSets() for true filtering.\n')
	#Load the different runs
	
	enrichList <- list()
	numExp <- NULL
	numExpNames <- NULL
	for (runName in runNames) {
		outFile.geneEnrich <- paste('./intermediate/', paste(runName, "geneExpr.Enrich.csv",sep="."), sep="")
		if (getClusterStack==TRUE) {outFile.geneEnrich <- sub('csv', 'clusterStack.csv', outFile.geneEnrich)}
		
		if (xval==FALSE) { numExpNames <- c(numExpNames, runName) }
		if (file.exists(outFile.geneEnrich) & !forceRebuild & xval==FALSE) {
			enrichList[[runName]] <- read.csv(outFile.geneEnrich, row.names=1)
		
			#Count the number of experiments in each dataset
			numExp[runName] <- length(list.files(paste('./rawData/', runName, sep=""), pattern='*.CEL', recursive=TRUE))
		} else {
			expData <- getFIRMAandGENE(runName, recreate=forceRebuild)

			geneExpr.vals <- expData$geneExpr.vals
			
			if (xval==FALSE) {
				testCols <- grep('HEMI', colnames(geneExpr.vals))
				controlCols <- grep('NCAR', colnames(geneExpr.vals))

				#If there are no test / controls, just pick some at random.
				#  This will happen if this is unlabeled.
				#But this doesn't seem to actually do anything?
				if (length(testCols) == 0 | length(controlCols) == 0) {
					cIdx <- grep('cell', colnames(geneExpr.vals))
					allCols <- (cIdx + 1):ncol(geneExpr.vals)
					half <- floor(length(allCols)/2)
					testCols <- allCols[1:half]
					controlCols <- allCols[(half+1):length(allCols)]
					nperms <- 10 #No reason since it's random
				} 

				x <- exonToGene(origData=geneExpr.vals)
				y <- as.numeric(grepl('NCAR', colnames(x))) + 1
				if(all(y == 1)) { y[floor(length(y)/2):length(y)] <- 2 }
				annotEnrich <- AnnotateGeneExpr(x, y, getClusterStack=getClusterStack, nperms=nperms)
				write.csv(annotEnrich, file=outFile.geneEnrich)

				enrichList[[runName]] <- annotEnrich

				#Count the number of experiments in each dataset
				numExp <- c(numExp, length(list.files(paste('./rawData/', runName, sep=""), pattern='*.CEL', recursive=TRUE)))
			}  else  {#if (xval==TRUE) {
				exclNames <- colnames(geneExpr.vals)[grepl('_HEMI', colnames(geneExpr.vals)) | grepl('_NCAR', colnames(geneExpr.vals))]
				for (exclName in exclNames) {
					curGeneExpr <- geneExpr.vals[,colnames(geneExpr.vals) != exclName]
					exclIdx <- which(colnames(geneExpr.vals) == exclName)
					curRunName <- paste(runName, exclIdx, sep='.')
				
					testCols <- grep('HEMI', colnames(curGeneExpr))
					controlCols <- grep('NCAR', colnames(curGeneExpr))

					curOutFile <- sub('Enrich', paste('Enrich.xval', exclIdx, sep='.'), outFile.geneEnrich)
					if (!file.exists(curOutFile))  {
						x <- exonToGene(origData=curGeneExpr)
						y <- as.numeric(grepl('NCAR', colnames(x))) + 1
						annotEnrich <- AnnotateGeneExpr(x, y, getClusterStack=getClusterStack, nperms = 300)
						write.csv(annotEnrich, file=curOutFile)
					} else {
						annotEnrich <- read.csv(curOutFile, row.names=1)
					}

					enrichList[[curRunName]] <- annotEnrich

					#Count the number of experiments in each dataset
					numExp <- c(numExp, length(exclNames)-1)
					numExpNames <- c(numExpNames, curRunName)
					#numExp[runName] <- length(list.files(paste('./rawData/', runName, sep=""), pattern='*.CEL', recursive=TRUE))
					#construct numExp sequentially, along with numExpNames
				} #for (exclName in colnames(geneExpr.vals)) {
			} #if (xval==FALSE) {
		} #if (file.exists()) {
	} #for (runName in runNames) {
	names(numExp) <- numExpNames

	#Optional: strip circadian and Riluzole pathways.  They were hand picked from 47 experiment data set
	if (stripCR == TRUE) {
		for (runName in names(enrichList)){
			pn <- rownames(enrichList[[runName]])
			enrichList[[runName]][!grepl('^js\\.', pn),]
		} #for (runName in names(enrichList)){
	}
	#01-09-13 USE Either ABS(GSA) or FDR.lo | FDR.hi to rank pathways
	#   New modifications will return FDR

	#Record the GSA scores Update 07-08-14 -- NOT IMPLEMENT, NOT APPROPRIATE HERE
	#There are three different methods here:
	#  A) GSAscore > 1 //The original way
	#  B) abs(GSAscore) > GSA.cutoff //The best way?
	#  C) pvalues.lo | pvalues.hi <= pVal.cut //The best way?
	#  D) A combination of B and C //The best way?
	# pVal.cut  Set to cut based on the pvalue.lo and the pvalue.hi (DEFAULT: NULL)
	# gsaCut  Remove pathways with gsa scores < gsaCut (DEFAULT: 1)
	# absGSAcut  Set to TRUE to use the absolute value (DEFAULT: FALSE)
	
	#outDF records the GSA scores.
	rNames <- unique(unlist(lapply(enrichList, function(x) {rownames(x)})))
	gsaScores <- matrix(0, nrow=length(rNames), ncol=length(enrichList))
	rownames(gsaScores) <- rNames
	colnames(gsaScores) <- names(enrichList)
	for(runName in colnames(gsaScores) ) {
		#Set anything that doesn't pass the p-Value cutoff to 0
		curGSA <- enrichList[[runName]]$gsaScore
		if (!is.null(pVal.cut)) {
			curGSA[(enrichList[[runName]]$pvalues.lo > pVal.cut) & (enrichList[[runName]]$pvalues.hi > pVal.cut)] <- 0
		}
		gsaScores[rownames(enrichList[[runName]]), runName] <- curGSA
		colnames(gsaScores)[colnames(gsaScores) == runName] <- paste(runName, '.', numExp[runName], 'exp', sep="")
	}
	
	outDF <- gsaScores
	curWeights <- numExp[sub('\\.[0-9]+exp', '', colnames(outDF))]
	curWeights <- curWeights / sum(curWeights)
	meanGSA <- apply(outDF, 1, function(x) { sum(x * curWeights) })

	if (is.null((pVal.cut))) {
		outDF <- cbind(data.frame(meanGSA=meanGSA, numGT1=apply(abs(outDF) > 1, 1, sum, na.rm=TRUE)), outDF, data.frame(genes=enrichList[[1]][rownames(outDF), 'genes']))
	} else {
		outDF <- cbind(data.frame(meanGSA=meanGSA, numGT1=apply(abs(outDF) > 0, 1, sum, na.rm=TRUE)), outDF, data.frame(genes=enrichList[[1]][rownames(outDF), 'genes']))
	}

	#Note there is a misonomer.  If there is a pVal.cut, then numGT1 should be numGT0
	if (absGSAcut == TRUE) {
		outDF <- cbind(data.frame(score=abs(outDF$meanGSA)+outDF$numGT1), outDF)
	} else {
		outDF <- cbind(data.frame(score=outDF$meanGSA+outDF$numGT1), outDF)		
	}

	outDF <- outDF[order(outDF$score, decreasing=TRUE), ]
	outFile.tally <- paste('./intermediate/geneExpr.Enrich.', paste(runNames, sep=".", collapse="."), ".gen.", Sys.Date(), ".csv", sep="")
	if (!is.null(pVal.cut)){outFile.tally <- sub('.csv', paste('.pCut', pVal.cut, 'csv', sep='.'), outFile.tally)}
	if (xval == TRUE) { outFile.tally <- sub('.csv', 'xval.csv', outFile.tally) }
	write.csv(outDF, file=outFile.tally)
	
	return(outDF)
}

#' Load the list of genesets
#'
#' @param geneNames  A vector of gene names
#' @param loadExtra  Set to FALSE to only load the GO list (DEFAULT: TRUE)
#' @param addCR  Set to TRUE to load circadian in riluzole (DEFAULT: FALSE)
#' @param getClusterStack  Set to TRUE to load only cMonkey clusterStack data (DEFAULT: FALSE)
#' @export
#' @usage genesets <- getGenesets(geneNames, loadExtra=TRUE, addCR=FALSE, getClusterStack=FALSE)
getGenesets <- function(geneNames, loadExtra=TRUE, addCR=FALSE, getClusterStack=FALSE) {
	if(getClusterStack == TRUE) {
		print('cMonkey intergration not implemented')
	} else {
	
		geneNames <- toupper(geneNames)
		goList <- GetHumanGOlists(geneNames)
		goList <- c(goList, GetBlists())

		if (loadExtra==FALSE) {
			genesets <- annotList <- goList
		} else {
			keggList <- GetKEGGlists(orgCode='hsa')
			reactomeList <- GetReactomeLists(orgCode='Homo sapiens')
			annotList <- c(goList, keggList, reactomeList)

			#Load the Immunology data set.  http://www.broadinstitute.org/gsea/msigdb/cards/GSE10325_CD4_TCELL_VS_LUPUS_CD4_TCELL_UP
			#'c7.all.v4.0.symbols.gmt'
			immunList <- GetImmunelists()
			annotList <- c(annotList, immunList)

			slc7 <- geneNames[grepl('SLC7', geneNames)]
			slc.top <- toupper(c("Slc7a1", "Slc7a3", "Slc7a5", "Slc7a11"))
			slc <- geneNames[grepl('SLC', geneNames)]
			olf <- geneNames[grepl('OLF', geneNames)]
			tas <- geneNames[grepl('TAS', geneNames)]
			tRNAsyn <- toupper(c("Cars", "Yars", "Aars", "Iars", "Gars", "Mars", "Eprs", "Lars", "Nars", "Sars", "Tars", "Wars2", "Vars", "Tars2", "Hars2", "Kars", "Sars2", "Lars2", "Ears2", "Rars", "Tarsl2", "Dars", "Iars2", "Cars2", "Aars2", "Aimp1", "Rars2", "Yars2", "Farsa", "Aimp2", "Mars2", "Prorsd1", "Nars2", "Pars2", "Vars2", "Dars2", "Qars", "Farsb", "Hars"))

			genesets <- c(list(sd.slc7=slc7, sd.slc.top=slc.top, sd.slc=slc, sd.olf=olf, sd.tas=tas, sd.tRNAsyn=tRNAsyn), annotList)
			#ALS specific pathways noticed by J Smith and added on 11-20-14, did not improve classifier
			if ( addCR == TRUE ) {
				circadian <- c('ARNTL', 'BHLHE40', 'BHLHE41', 'CLOCK', 'CRY1', 'CRY2', 'NR1D1', 'NR1D2', 'PER1', 'PER2', 'PER3', 'GM129', 'RORA', 'RORB', 'RORC', 'DBP', 'NFIL3')
				riluzole.sodium <- c('SCN5A', 'SCN2A1', 'SCN4A', 'SCN1B', 'SCN2B', 'SCN1A', 'SCN8A', 'SCN9A', 'SCN3B', 'SCN4B', 'SCN3A', 'SCN7A', 'SCN11A', 'SCN10A')
				riluzole.glutamate <- c('SLC1A6', 'SLC1A2', 'SLC1A3', 'SLC1A1', 'SLC1A4', 'SLC1A7')
				genesets <- c(list(js.circadian=circadian, js.riluzole.sodium=riluzole.sodium, js.riluzole.glutamate=riluzole.glutamate), genesets)
			} #if ( addCR == TRUE ) {
		} #if (loadExtra==FALSE) {
		names(genesets) <- make.names(names(genesets), unique = TRUE)
	} #if(getClusterStack == TRUE) {
	return(genesets)
}

#' This function should generate the training and test set for ML
#'  NOTE:  This function has not been finalized, rather it is intended to give a clue to the user
#' 
#' Note, Method 3 lets different days be swapped in and out for training and testing sets
#'   Will return NAs for pathways that do not have enough genes (i.e. atleast 3 or 4)
#'
#' @param x  The raw gene expression scores
#' @param y  The ranked index to the exons to use as features (DEFAULT: NULL, i.e. all exons)
#' @param cutoff  A False discovery rate cutoff.  Results not sane 010913 (DEFAULT: 1)
#' @param nperms  nperms parameter for GSA function (DEFAULT: 300)
#' @param loadExtra  Set to FALSE to only load the GO list (DEFAULT: TRUE)
#' @param genesets  A pre-loaded geneset from getGenesets (DEFAULT: NULL)
#' @param getClusterStack  Set to TRUE to load only cMonkey clusterStack data (DEFAULT: FALSE)
#' @export
#' @usage annotEnrich <- AnnotateGeneExpr(x, y, cutoff=1, nperms=300, loadExtra=TRUE, genesets=NULL, getClusterStack=FALSE)
AnnotateGeneExpr <- function(x, y, cutoff=1, nperms=300, loadExtra=TRUE, genesets=NULL, getClusterStack=FALSE) {
	library(GSA)
	
	rownames(x) <- toupper(rownames(x))
	if (is.null(genesets)) {
		genesets <- getGenesets(geneNames=rownames(x), loadExtra=loadExtra, getClusterStack=getClusterStack)
	}
	
	#Filter out gene sets to only include those genes that appear in x
	for (i in 1:length(genesets)) {
		genesets[[i]] <- genesets[[i]][genesets[[i]] %in% rownames(x)]
	}
	#sapply(genesets, length)
	GSAres <- GSA(as.matrix(x), y, genesets=genesets, genenames=rownames(x), minsize=4, maxsize=max(sapply(genesets, length))+1, nperms=nperms)

	#Make the gene lists
	geneLists <- sapply(genesets, function(x) {paste(unique(x), collapse="; ")})
	
	#  These can be used to determine a p-value cutoff that corresponds to FDR
	#fdrData <- GSA.listsets(GSAres, geneset.names=1:length(genesets), FDRcut=cutoff)
	#fdrSetIdx <- c(fdrData$negative[,'Gene_set_name'], fdrData$positive[,'Gene_set_name'])
	#fdrs <- as.numeric(c(fdrData$negative[,'FDR'], fdrData$positive[,'FDR']))
	#fdrSetNames <- names(genesets)[as.numeric(fdrSetIdx)]
	#GSAfdrs <- data.frame(setNames=fdrSetNames, fdrs=fdrs)
	
	GSAscores <- data.frame(gsaScore = GSAres$GSA.scores, pvalues.lo=GSAres$pvalues.lo, pvalues.hi=GSAres$pvalues.hi)
	rownames(GSAscores) <- names(genesets)
	#GSAscores <- cbind(GSAscores[GSAfdrs$setNames,], data.frame(fdr=GSAfdrs$fdrs))
	GSAscores <- GSAscores[order(abs(GSAscores$gsaScore), decreasing=TRUE),]
	GSAscores <- cbind(GSAscores, data.frame(genes=geneLists[rownames(GSAscores)]))
	return(GSAscores)
}

#' Load the B-cell specific pathways
#'
#' @param inFile  Input file type (DEFAULT: './intermediate/B.Pathways.csv')
#' @param recalc  Set to TRUE to force the recalculation of the list (DEFAULT: FALSE)
#' @export
#' @usage annotList <- GetBlists(inFile='./intermediate/B.Pathways.csv', recalc=FALSE)
GetBlists <- function(inFile='./intermediate/B.Pathways.csv', recalc=FALSE) {
	outFile <- paste('./intermediate/', paste("bList.RData",sep="."), sep="")

	if (file.exists(outFile) & (recalc==FALSE)) {
		annotList <- get(load(outFile)[1])
	} else {
		immunDF <- read.table(inFile, sep=',', header=TRUE, fill=TRUE)
		annotList <- list()

		for (i in 1:nrow(immunDF)) {
			annotName <- paste('BCELL.', as.character(immunDF[i,1, drop=TRUE]), sep="")
			genes <- toupper(strsplit(as.character(immunDF[i,4]), '; *')[[1]])
			genes <- genes[!grepl('^N/A', genes)]
			annotList[[annotName]] <- genes
		}
		save(annotList, file=outFile)
	} #if (file.exists(outFile) & (recalc==FALSE)) {
	
	return(annotList)
}

#' Load the Immunology data set.  http://www.broadinstitute.org/gsea/msigdb/index.jsp
#'
#' @param iFile  The gmt file name. (DEFAULT: './genesetAnnots/c7.all.v4.0.symbols.gmt')
#' @param recalc  Set to TRUE to force recalculation. (DEFAULT: FALSE)
#' @export
#' @usage annotList <- GetImmunelists(iFile='./genesetAnnots/c7.all.v4.0.symbols.gmt')
GetImmunelists <- function(iFile='./genesetAnnots/c7.all.v4.0.symbols.gmt', recalc=FALSE) {
	outFile <- paste('./intermediate/', paste("immuneList.RData",sep="."), sep="")

	if (file.exists(outFile) & (recalc==FALSE)) {
		annotList <- get(load(outFile)[1])
	} else {
		immunDF <- read.table(iFile, sep='\t', header=FALSE, fill=TRUE)
		annotList <- list()

		for (i in 1:nrow(immunDF)) {
			annotName <- as.character(immunDF[i,1, drop=TRUE])
			genes <- as.character(unlist(immunDF[i,3:ncol(immunDF),drop=TRUE]))
			annotList[[annotName]] <- genes
		}
		save(annotList, file=outFile)
	} #if (file.exists(outFile) & (recalc==FALSE)) {
	
	return(annotList)
}

#' This function should generate the training and test set for ML
#'
#' @param geneNames  The names of the genes (DEFAULT: NULL)
#' @param recalc  Set to true to recalculate the list (DEFAULT: FALSE)
#' @export
#' @usage annotEnrich <- GetHumanGOlists(geneNames=NULL, recalc=FALSE)
GetHumanGOlists <- function(geneNames=NULL, recalc=FALSE) {
	outFile <- paste('./intermediate/', paste("goList.RData",sep="."), sep="")

	if (file.exists(outFile) & (recalc==FALSE)) {
		goList <- get(load(outFile)[1])
	} else {
		library("hgu95av2.db")
		library('GO.db')
		#keytypes(hgu95av2.db)
	
		#Get GO Annotations for genes
		geneNames.all <- keys(hgu95av2.db, "ALIAS")
		geneNames.all <- geneNames.all[!grepl('\'', geneNames.all)] #Fix a bug
		#tst <- select(hgu95av2.db, keys, "GO", "ALIAS")
		goIDs <- suppressWarnings(select(hgu95av2.db, geneNames.all, "GO", "ALIAS"))[,c('ALIAS', 'GO')]
		goLUT <- suppressWarnings(select(GO.db, unique(goIDs$GO), "TERM", "GOID"))
		goList <- list()
		for (goID in unique(goIDs$GO)) {
			if (is.na(goID)) {next;}
			goName <- goLUT[goLUT$GOID==goID, 'TERM']
			goList[[goName]] <- unlist(goIDs$ALIAS[which(goIDs$GO==goID)])
		}
		save(goList, file=outFile)
	}	

	if (!is.null(geneNames)) {
		for (goName in names(goList)) {
			goList[[goName]] <- goList[[goName]][goList[[goName]] %in% geneNames]
		} #for (goName in names()) {
		goList <- goList[sapply(goList, function(x) {length(x) != 0})]
	}
		
	return(goList)
}

#' This function should generate the training and test set for ML
#'   Note: Using HUMAN gene names even if you use a non-human pathway
#'
#' @param orgCode  The organism name for the KEGG pathway (DEFAULT: 'hsa', i.e. human)
#' @param recalc  Set to true to recalculate the list (DEFAULT: FALSE)
#' @export
#' @usage keggList <- GetKEGGlists(orgCode='hsa', recalc=FALSE)
GetKEGGlists <- function(orgCode='hsa', recalc=FALSE) {
	outFile <- paste('./intermediate/', paste("keggList.RData",sep="."), sep="")

	if (file.exists(outFile) & (recalc==FALSE)) {
		keggGenes <- get(load(outFile)[1])
	} else {
		library("hgu95av2.db")
		library('KEGG.db')
		#library('reactome.db')
		#keytypes(hgu95av2.db)
	
		keggGenes <- as.list(KEGGPATHID2EXTID)
		keggGenes <- keggGenes[grepl(orgCode, names(keggGenes))]
	
		#Get GO Annotations for genes
		for (keggName in names(keggGenes)) {
			geneNames <- try(suppressWarnings(select(hgu95av2.db, keggGenes[[keggName]], "ALIAS", "ENTREZID"))[,c('ALIAS')], silent=TRUE)
			keggGenes[[keggName]] <- geneNames[!is.na(geneNames)]
		}
	
		keggNames <- unlist(as.list(KEGGPATHID2NAME))
		curNames <- keggNames[sub( orgCode, '', names(keggGenes))]
		names(keggGenes) <- sub('^', 'KEGG.', make.names(curNames))
		keggGenes <- keggGenes[sapply(keggGenes, function(x) {!inherits(x, "try-error")})]
		save(keggGenes, file=outFile)
	} # if (file.exists(outFile) & (recalc==FALSE)) {	
	return(keggGenes)
}

#' This function should generate the training and test set for ML
#'   Note: Using HUMAN gene names even if you use a non-human pathway
#'   DOESN'T WORK AS OF 04-24-14
#'
#' @param orgCode  The organism name for the KEGG pathway (DEFAULT: 'hsa', i.e. human)
#' @param recalc  Set to true to recalculate the list (DEFAULT: FALSE)
#' @export
#' @usage keggList <- GetKEGGRESTlists(orgCode='hsa', recalc=FALSE)
GetKEGGRESTlists <- function(orgCode='hsa', recalc=FALSE) {
	outFile <- paste('./intermediate/', paste("keggrestList.RData",sep="."), sep="")

	if (file.exists(outFile) & (recalc==FALSE)) {
		keggGenes <- get(load(outFile)[1])
	} else {
		library("hgu95av2.db")
		library('KEGGREST')
		#library('reactome.db')
		#keytypes(hgu95av2.db)

		keggLists <- keggList(orgCode)
		keggLUT <- sapply(keggLists, function(x) {sub(';.*$', '', sub(',.*$', '', x))})
	
		#pathGenes <- keggLink("pathway", "hsa")
		genePath <- keggLink("hsa", "pathway")
		
		keggGenes <- list()
		for(curPath in unique(names(genePath))) {
			keggGenes[[curPath]] <- keggLUT[genePath[[curPath]]]
		}
		
		save(keggGenes, file=outFile)
	} # if (file.exists(outFile) & (recalc==FALSE)) {	
	return(keggGenes)
}

#' This function should generate the training and test set for ML
#'   Note: Using HUMAN gene names even if you use a non-human pathway
#'
#' @param orgCode  The organism name for the KEGG pathway (DEFAULT: 'Homo sapiens', i.e. human)
#' @param recalc  Set to true to recalculate the list (DEFAULT: FALSE)
#' @export
#' @usage reactomeList <- GetReactomeLists(orgCode='Homo sapiens', recalc=FALSE)
GetReactomeLists <- function(orgCode='Homo sapiens', recalc=FALSE) {
outFile <- paste('./intermediate/', paste("reactomeList.RData",sep="."), sep="")

	if (file.exists(outFile) & (recalc==FALSE)) {
		reactomeGenes <- get(load(outFile)[1])
	} else {
		library("hgu95av2.db")
		#library('KEGG.db')
		library('reactome.db')
		#keytypes(hgu95av2.db)
	
		reactomeGenes <- as.list(reactomePATHID2EXTID)
		reactomeNames <- unlist(as.list(reactomePATHID2NAME))
		names(reactomeGenes) <- make.names(reactomeNames[names(reactomeGenes)])
		reactomeGenes <- reactomeGenes[grepl( make.names(orgCode), names(reactomeGenes))]
	
		#Get GO Annotations for genes
		for (reactomeName in names(reactomeGenes)) {
			geneNames <- try(suppressWarnings(select(hgu95av2.db, reactomeGenes[[reactomeName]], "ALIAS", "ENTREZID"))[,c('ALIAS')], silent=TRUE)
			reactomeGenes[[reactomeName]] <- geneNames[!is.na(geneNames)]
		}
	
		names(reactomeGenes) <- sub('^', 'reactome.', names(reactomeGenes))
		reactomeGenes <- reactomeGenes[sapply(reactomeGenes, function(x) {!inherits(x, "try-error")})]
		save(reactomeGenes, file=outFile)
	} #if (file.exists(outFile) & (recalc==FALSE)) {
	return(reactomeGenes)
}


#' Run FIRMA and RMA on a data set using the doFIRMA and doRMA functions
#'   NOTE: This uses the latest doFIRMA and doRMA, other functions are obsolete
#'
#' @param dataSet  The exon expression array (DEFAULT: 'ALS' )
#'                 Alternative: ALS.all.122613
#' @param subSet Set to a range of numbers to only process a subset of the data (DEFAULT: NULL, i.e. 'all')
#'          NOTE: Not implemented 12/6/13
#' @param recreate Set to TRUE to prevent loading saved file(DEFAULT: FALSE)
#' @param chipType Usually "MoEx-1_0-st-v1" for mouse or 'HuEx-1_0-st-v2' for human (DEFAULT: "MoEx-1_0-st-v1")
#' @param tags Used for file nameing, (DEFAULT: NULL, i.e. mouse: 'U-Ensembl50,G-Affy,EP', human: 'coreR2,A20070914,EP')
#' @param loadFIRMA Set to FALSE to not load the FIRMA scores (DEFAULT: TRUE)
#' @export
#' @usage fsScores <- getFIRMAandGENE(dataSetName='ALS.new.122012', subSet=NULL, recreate=FALSE, chipType="MoEx-1_0-st-v1", tags=NULL, loadFIRMA=TRUE)
getFIRMAandGENE <- function(dataSetName='ALS.new.122012', subSet=NULL, recreate=FALSE, chipType="MoEx-1_0-st-v1", tags=NULL, loadFIRMA=TRUE) {
	library(affy)     # bioconductor
	library(Biobase)  # bioconductor
	library(gplots)   # CRAN
	library(gdata)    # CRAN
	library(limma)    # bioconductor
	library(GenomeGraphs)   # bioconductor
	library(aroma.affymetrix)

	if (is.null(tags)) {
		if (chipType == "MoEx-1_0-st-v1") {
			tags <- 'U-Ensembl50,G-Affy,EP'
		} else if (chipType == "HuEx-1_0-st-v2") {
			tags <- 'coreR2,A20070914,EP'
		} else {
			tags <- ""
		}
		cat('tags set to \'', tags, "\'\n", sep="")
	}

	if(!is.null(subSet)) { 
		sTag <- paste('*', paste(subSet, collapse='.'), sep=",") #This assumes a continuous tag
	} else {
		sTag <- ""
	}

	outFile.gene <- paste('./intermediate/', paste(dataSetName, chipType, tags, sTag, "exonVals.RData",sep="."), sep="")
	outFile.exon <- paste('./intermediate/', paste(dataSetName, chipType, tags, sTag, "firmaVals.RData",sep="."), sep="")
	
	if (!file.exists(outFile.gene) | !file.exists(outFile.exon) | recreate == TRUE) {
		verbose <- Arguments$getVerbose(-8, timestamp=TRUE)  # specify full messages to the R console
		#From http://aroma-project.org/chipTypes/MoEx-1_0-st-v1
		cdf <- AffymetrixCdfFile$byChipType(chipType, tags=tags)
		#print(cdf)
		cs <- AffymetrixCelSet$byName(dataSetName, cdf=cdf)

		#print(cs)
		setCdf(cs, cdf)

		#Load the standard gene names
		geneExonLUT <- loadGeneExonLUT(cdf)
	}
	
	if (!file.exists(outFile.gene) | recreate == TRUE) {
		#   03-24-14: Added, recreate == TRUE
		geneExpr <- doRMA(cs) #, addNames=TRUE)
		geneExpr.vals <- extractDataFrame(geneExpr, addNames=TRUE)
		rm(geneExpr)
		gc(reset=TRUE)
		LUTidxs <- match(as.character(geneExpr.vals$groupName), as.character(geneExonLUT$stable_id))
		if(all(is.na(LUTidxs))) {
			LUTidxs <- match(as.character(geneExpr.vals$unitName), as.character(geneExonLUT$stable_id))
		}
		
		geneExpr.vals <- cbind(geneExonLUT[LUTidxs, c("symbol","symbol_description")], geneExpr.vals)
		save(geneExpr.vals, file=outFile.gene)
	} else {
		cat('Loading:', outFile.gene, '\n')
		geneExpr.vals <- get(load(outFile.gene)[1])
	}
	
	if (loadFIRMA == TRUE) {
		if (!file.exists(outFile.exon) | (recreate == TRUE)) {
			cdf <<- cdf  #Work around for FIRMA bug
			exonExpr <- doFIRMA(cs) #, addNames=TRUE)
			exonExpr.vals <- extractDataFrame(exonExpr, addNames=TRUE)
			rm(exonExpr)
			gc(reset=TRUE)
			LUTidxs <- match(as.character(geneExpr.vals$groupName), as.character(geneExonLUT$stable_id))
			if(all(is.na(LUTidxs))) {
				LUTidxs <- match(as.character(geneExpr.vals$unitName), as.character(geneExonLUT$stable_id))
			}

			exonExpr.vals <- cbind(geneExonLUT[LUTidxs, c("symbol","symbol_description")], exonExpr.vals)
			save(exonExpr.vals, file=outFile.exon)
		} else {
			cat('Loading:', outFile.exon, '\n')
			exonExpr.vals <- get(load(outFile.exon)[1])
		} #if (!file.exists(outFile.exon)  | recreate == TRUE) {
	} else {
		exonExpr.vals <- NULL
	} #if (loadFIRMA == TRUE) {
	
	outVals <- list(geneExpr.vals=geneExpr.vals, exonExpr.vals=exonExpr.vals)
	return(outVals)
}


#' Load the gene exon LUT.
#'
#' @param cdf  The cdf file with the gene probe name
#' @export
#' @usage geneExonLUT <- loadGeneExonLUT(cdf=NULL)
loadGeneExonLUT <- function(cdf=NULL) {
	
	if (is.null(cdf)) {
		chipType <- "MoEx-1_0-st-v1"
		tags <- 'U-Ensembl50,G-Affy,EP'
		#From http://aroma-project.org/chipTypes/MoEx-1_0-st-v1
		cdf <- AffymetrixCdfFile$byChipType(chipType, tags=tags)
		#print(cdf)
		cs <- AffymetrixCelSet$byName('ALS', cdf=cdf)
	
		#print(cs)
		setCdf(cs, cdf)
	}
	
	cdfName <- capture.output(cdf)
	cdfName <- sub('^.*: ', '', cdfName[grepl('^Chip', cdfName)])
	saveFile <- paste('./intermediate/LUT', cdfName, 'RData', sep=".")

	if(file.exists(saveFile)) {
		naData <- get(load(saveFile)[1])
	} else {
		#To look up the names, I need a csv file from affy.
		#http://www.aroma-project.org/vignettes/UsingGenomeGraphsWithFIRMA
		chipType <- getChipType(cdf, fullname=FALSE);
		na <- AffymetrixNetAffxCsvFile$byChipType(chipType, pattern=".*.probeset.csv$");
		colPatterns <- c("(probesetId|seqname|strand|transcriptClusterId)"="character", "(start|stop)"=NA);
		#naData <- readDataFrame(na, colClassPatterns=colPatterns);
		naData <- readDataFrame(na, colClasses=colPatterns); #update 12/6/13
		dim(naData);

		na.trans <- AffymetrixNetAffxCsvFile$byChipType(chipType, pattern=".*.transcript.csv$");
		#"transcript_cluster_id","probeset_id","seqname","strand","start","stop","total_probes","gene_assignment","mrna_assignment","swissprot","unigene","GO_biological_process","GO_cellular_component","GO_molecular_function","pathway","protein_domains","category"
		#colPatterns <- c("(transcriptClusterId|seqname|strand|geneAssignment|mrnaAssignment|swissprot|unigene)"="character", "(start|stop)"=NA);
		naData.trans <- readDataFrame(na.trans);
		dim(naData.trans);

		LUTidxs <- match(as.character(naData$transcriptClusterId), as.character(naData.trans$transcriptClusterId))
		naData <- cbind(naData, naData.trans[LUTidxs,])

		#Duplicate this behavior to get symbol / gene names
		#geneExonLUT <- loadGeneExonLUT()
		#LUTidxs <- match(as.character(fsScores$unitName), as.character(geneExonLUT$stable_id))
		#geneAnnots <- geneExonLUT[LUTidxs, c("symbol","symbol_description")]
		# Set up BioMart
		#mart <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl");
		# Transcript annotations
		#tst <- getBM( filters= "ensembl_gene_id", attributes= c("ensembl_gene_id", "external_gene_id", 
		#	"entrezgene", "description"), values= ensembl_genes, mart= mart)
		colnames(naData)[colnames(naData) == 'probesetId'] <- 'stable_id' #042914 -- why?
		symbol <- sapply(naData$geneAssignment, function(x) { gsub(' ', '', as.character(strsplit(x, '//')[[1]][2])) })
		symbol_description <- sapply(naData$geneAssignment, function(x) { gsub(' ', '', as.character(strsplit(x, '//')[[1]][3])) })
		names(symbol) <- names(symbol_description) <- NULL
		
		naData <- naData[,1:(which(colnames(naData) == 'geneAssignment')-1)]
		#Remove the extraneous information
		
		#Duplicate this behavior to get symbol / gene names
		#LUTidxs <- match(as.character(fsScores$unitName), as.character(geneExonLUT$stable_id))
		#geneAnnots <- geneExonLUT[LUTidxs, c("symbol","symbol_description")]
		naData$symbol <- symbol
		naData$symbol_description <- symbol_description

		save(naData, file=saveFile)
	} #if(file.exists(saveFile)) {
	
	return(naData)
}

#' Generate estimate FWER for ranked features using the mProbes method outlined in 
#'  Van Anh Huynh-Thu, Yvan Saeys, Louis Wehenkel, and Pierre Geurts
#'  Statistical interpretation of machine learning-based feature importance scores for biomarker discovery
#'  Bioinformatics Advance Access published April 25, 2012
#'  See notebook #3, page 164 for algorithm
#'  NOTE:  mProbes is pretty strict, so it wouldn't be aweful to select ANY that have
#'    a score less than 1.
#'  Can I use CMIM?
#'  rFerns should be two orders of magnitude faster than RandomForest but similarly ordered. 
#'
#' @param dataMat  The attribute matrix
#' @param classIdx  The index for the class variable (DEFAULT: ncol(dataMat))
#' @param type  The feature selection:  Currently 'MI', 'randomForest', 'rFerns', or 'Dynamic'.  Dynamic sets MI if there's more than 700 features, or randomForest otherwise. (DEFAULT: "Dynamic")
#' @param dynamicCutoff  The cutoff if 'type' is set to "Dynamic". (DEFAULT: 800)
#' @param ntree  The number of trees to include in a random forest or ferns in rFerns if used. (DEFAULT: 1000)
#' @param mcWorkers  Set to override the number of threads. (DEFAULT: getDoParWorkers())
#' @return a ranked list of features with their estimated FWER
#' @export
#' @usage features.FWER <- mProbes(dataMat, classIdx=ncol(dataMat), type="Dynamic", dynamicCutoff=800, ntree=1000, mcWorkers=getDoParWorkers())
mProbes <- function (dataMat, classIdx=ncol(dataMat), type="Dynamic", dynamicCutoff=800, ntree=1000, mcWorkers=getDoParWorkers()) {
	gc(reset=TRUE)
	library(randomForest)
	library(rFerns)

	classDF <- dataMat[, classIdx, drop=FALSE]
	className <- colnames(dataMat)[classIdx]
	dataPart <- dataMat[,1:ncol(dataMat) != classIdx]
	features.FWER <- rep(1, ncol(dataPart))
	names(features.FWER) <- colnames(dataPart)

	#Dynamic select MI method.  Used MI only if there's too many features from randomForest
	#Change 05-13-15 to use rFerns if there are too many features
	if (type == "Dynamic") {
		if (ncol(dataPart) > dynamicCutoff) {
			#type <- 'MI'
			type <- 'rFerns'
		} else {
			type <- 'randomForest'
		}
	}

	if (type == 'randomForest' | type == 'rFerns') {
		numFeat <- ncol(dataPart)
		if (ntree <= 2 * numFeat) {
			cat(numFeat, ' features but only ', ntree, ' trees/ferns for ', type , '. Consider using more trees/ferns.', '\n', sep="")
		}
	}

	cat('\nRunning mProbes with feature selector: ', type, ' on ', ncol(dataPart), ' candidates\n', sep="")

	# I. Use randomForest to rank features 
	# For Xi in all X features
	#   1) Generate |X| permutations of Xi
	#   2) Rank the 2*|X| features
	#rankedFeat <- foreach (i = 1:ncol(dataPart)) %dopar% {
	####
	oldWorkers <- getDoParWorkers()
	if (mcWorkers != oldWorkers) { registerDoMC(cores = mcWorkers) }
	rankedFeat <- foreach (i = 1:ncol(dataPart)) %dopar% {
		if(i %% ceiling(ncol(dataPart)/100) == 0) { cat(round(i/ncol(dataPart),2)*100,'%,',sep="") }

		Xi <- dataPart[,i]
		randData <- as.matrix(sapply(1:ncol(dataPart), function(x) {Xi[order(runif(length(Xi)))]}))
		rownames(randData) <- rownames(dataPart)
		colnames(randData) <- paste('random.', 1:ncol(randData), sep="")
		curExp <- cbind(dataPart, randData)
		curExp <- apply(curExp, 2, as.numeric)
		rownames(curExp) <- rownames(dataPart)
		curExp <- cbind(as.data.frame(curExp), classDF)
		curForm <- as.formula(paste(colnames(curExp)[ncol(curExp)], "~", "."))
		if (type == "MI") {
			featureList <- xvalMIranker(curExp, classIdx=nrow(curExp), folds=2, verbose=FALSE)
			features.ranked <- featureList[,1]
			names(features.ranked) <- rownames(featureList)
			features.ranked <- sort(features.ranked, decreasing=TRUE)
		} else if (type == 'rFerns') {
			#if(i %% ceiling(ncol(dataPart)/100) == 0) {cat(' rferns ')}
			rf <- rFerns(x=curExp[,1:(ncol(curExp)-1),drop=FALSE], y=curExp[,ncol(curExp)], ferns=ntree, importance=TRUE, saveForest=FALSE)
			inf <- rf$importance[,'MeanScoreLoss']
			names(inf) <- rownames(rf$importance)
			features.ranked <-  sort(inf, decreasing=FALSE)  
		} else if (type == 'randomForest') {
			#rf <- randomForest(curForm, data=curExp)
			rf <- randomForest(x=curExp[,1:(ncol(curExp)-1),drop=FALSE], y=curExp[,ncol(curExp)], ntree=ntree, importance=TRUE,  keep.forest=FALSE)
			#rf <- randomForest(curForm, data=curExp, ntree=20)
			features.ranked <-  importance(rf)[order(importance(rf)[,'MeanDecreaseAccuracy']),'MeanDecreaseAccuracy']
		} else {
			cat(type, 'is not a valid type of ranker\n')
			return(NULL)
		}

		#Enumerate the features that are better than random
		#Most important features should be on the bottom
		bestRandom <- grep('^random\\.[0-9]*$', names(features.ranked))
		bestRandom <- bestRandom[length(bestRandom)]
		if (bestRandom < length(features.ranked)) {
			features <- names(features.ranked)[(bestRandom+1):length(features.ranked)]
		} else {
			features <- NULL
		}

		return(features)
	} #rankedFeat <- foreach (i = 1:ncol(dataPart)) %dopar% {
	if (mcWorkers != oldWorkers) { registerDoMC(cores = oldWorkers) }
	cat('\n')

	# II. Calculate the FWERs using the previous list
	# For Xi in all X features
	#   1) FWER = fraction of the |X| runs where a rand feature ranks above Xi
	cur.FWERs <- 1 - (table(unlist(rankedFeat)) / length(rankedFeat))
	features.FWER[match(names(cur.FWERs), names(features.FWER))] <- cur.FWERs

	return(features.FWER)
}


#' Make predictions using multiple methods
#'
#' @param trainData  The training data matrix
#' @param testData  The test data matrix
#' @param numsFeatures  The number of top features to test (DEFAULT: 1:100)
#' @param onlySVM  SEt to TRUE to only run the support vector machine (DEFAULT: FALSE)
#' @param posClass  The class of positive examples (DEFAULT: 'HEMI')
#' @param negClass  The class of negative examples (DEFAULT: 'NCAR')
#' @param singleCore  Set to TRUE to only use a single core (DEFAULT: FALSE)
#' @return the predictions
#' @export
#' @usage predResults <- makePreds(trainData, testData, numsFeatures=1:100, onlySVM=FALSE, posClass='HEMI', negClass='NCAR', singleCore=FALSE)
makePreds <- function (trainData, testData, numsFeatures=1:100, onlySVM=FALSE, posClass='HEMI', negClass='NCAR', singleCore=FALSE) {
	library(e1071)
	library(evtree)
	library(rpart)
	library(kernlab)
	library(BayesTree, quietly = TRUE)
	
	print('To Do: Add BayesTree to the ensemple.\nKueffner et al. Nature Biotechnology 33, 51-57 (2015)\n')
	
	
	labelIdx <- ncol(trainData)
	colnames(trainData)[labelIdx] <- 'Class'
	colnames(testData)[labelIdx] <- 'Class'
	trainData$Class <- as.factor(as.character(trainData$Class))
	if(max(numsFeatures) >= ncol(trainData)) { numsFeatures <- min(numsFeatures):(ncol(trainData)-1) }
	
	if (singleCore == FALSE) {
		#for (num in numsFeatures) {
		allRes <- foreach (num = numsFeatures) %dopar% {
			results <- list()

			cat(num, ", ", sep="")

			svmModel <-  runSVM(trainData[,c(1:num,labelIdx),drop=FALSE])
			results[['SVM']] <-  predict(svmModel, newdata = as.matrix(testData[,1:num, drop=FALSE]))
			names(results[['SVM']]) <- rownames(testData)

			nbModel <- naiveBayes(trainData[,1:num], trainData[,labelIdx]) 
			results[['NB']] <- predict(nbModel, newdata = as.matrix(testData[,1:num, drop=FALSE]))
			names(results[['NB']]) <- rownames(testData)

			#NOTE: For evtree and rf, it may be best to tune speed based on data size.
			minsplit <- min(nrow(trainData)-1, 20)
			ev.control=evtree.control(niterations = 90000, ntrees = 300, alpha=0.1, minsplit=minsplit)
			evModel <- try(evtree(formula=as.formula('Class ~ .'), data=trainData[,c(1:num,labelIdx),drop=FALSE], control=ev.control), silent = TRUE)
			if (class(evModel)[1] == "try-error") {
				cat('EVTREE failed, using random labels\n')
				#Assign labels at random, could be done better
				newRes <- rep(levels(trainData[,labelIdx])[1], nrow(testData))
				frac2 <- mean(trainData[,labelIdx] == levels(trainData[,labelIdx])[2])
				newRes[runif(length(newRes)) < frac2] <- levels(trainData[,labelIdx])[2]
				results[['EV']] <- factor(newRes)
			} else {
				results[['EV']] <- predict(evModel, newdata = testData[,1:num, drop=FALSE])
			}
			names(results[['EV']]) <- rownames(testData)

			rfModel <- runRF(trainData[,c(1:num,labelIdx),drop=FALSE], numTrees=10000, mtry=max(ceiling(num/3), 1))
			results[['RF']] <- predict(rfModel, newdata = testData[,1:(ncol(testData)-1), drop=FALSE])
			names(results[['RF']]) <- rownames(testData)

			gc(reset=TRUE)
			return(results)
		} #Foreach
	} else {	
		allRes <- foreach (num = numsFeatures) %do% {
			results <- list()

			cat(num, ", ", sep="")

			svmModel <-  runSVM(trainData[,c(1:num,labelIdx),drop=FALSE])
			results[['SVM']] <-  predict(svmModel, newdata = as.matrix(testData[,1:num, drop=FALSE]))
			names(results[['SVM']]) <- rownames(testData)

			nbModel <- naiveBayes(trainData[,1:num], trainData[,labelIdx]) 
			results[['NB']] <- predict(nbModel, newdata = as.matrix(testData[,1:num, drop=FALSE]))
			names(results[['NB']]) <- rownames(testData)

			#NOTE: For evtree and rf, it may be best to tune speed based on data size.
			minsplit <- min(nrow(trainData)-1, 20)
			ev.control=evtree.control(niterations = 90000, ntrees = 300, alpha=0.1, minsplit=minsplit)
			evModel <- try(evtree(formula=as.formula('Class ~ .'), data=trainData[,c(1:num,labelIdx),drop=FALSE], control=ev.control), silent = TRUE)
			if (class(evModel)[1] == "try-error") {
				cat('EVTREE failed, using random labels\n')
				#Assign labels at random, could be done better
				newRes <- rep(levels(trainData[,labelIdx])[1], nrow(testData))
				frac2 <- mean(trainData[,labelIdx] == levels(trainData[,labelIdx])[2])
				newRes[runif(length(newRes)) < frac2] <- levels(trainData[,labelIdx])[2]
				results[['EV']] <- factor(newRes)
			} else {
				results[['EV']] <- predict(evModel, newdata = testData[,1:num, drop=FALSE])
			}
			names(results[['EV']]) <- rownames(testData)

			rfModel <- runRF(trainData[,c(1:num,labelIdx),drop=FALSE], numTrees=10000, mtry=max(ceiling(num/3), 1))
			results[['RF']] <- predict(rfModel, newdata = testData[,1:(ncol(testData)-1), drop=FALSE])
			names(results[['RF']]) <- rownames(testData)

			gc(reset=TRUE)
			return(results)
		} #Foreach
	} #if (singleCore == FALSE) {
	names(allRes) <- numsFeatures

	results <- list()
	results$NB <- list()
	results$SVM <- list()
	results$EV <- list()
	results$RF <- list()
	for(n in names(allRes)) {
		results$SVM[[n]] <- allRes[[n]]$SVM
		results$NB[[n]] <- allRes[[n]]$NB
		results$EV[[n]] <- allRes[[n]]$EV
		results$RF[[n]] <- allRes[[n]]$RF
	}


	cat('\n')
	
	#Count the predictions
	predTally <- list()
	for (className in names(results)) {
		curRes <- unlist(results[[className]], use.names=TRUE)
		if (!all(is.na(curRes))) {
			names(curRes) <- sub('^.*\\.', '', names(curRes))
			tRes <- tapply(curRes, names(curRes), table)
			resDF <- foreach(res = tRes, .combine=rbind) %do% { res }
			if (length(tRes) == 1) {  #Work around to force a DF
				resDF <- rbind(resDF, resDF)[1,,drop=FALSE]
			}
			rownames(resDF) <- names(tRes)
			resDF <- as.data.frame(resDF)
			pC1 <- resDF[, posClass] / rowSums(resDF)
			resDF$pHEMI <- pC1
			resDF$prediction <- negClass
			resDF$prediction[resDF$pHEMI > 0.5] <- posClass
			predTally[[className]] <- resDF
		}
	}
	
	#Count aggregate predictions
	for(resName in names(results)) {
		res <- results[[resName]]
		if( all(sapply(res, function(x) { all(is.na(x)) })) ) {
			results[[resName]] <- NULL
		}
	}
	
	curRes <- unlist(foreach(res = results, .combine=c) %do% { if(!all(is.na(res))) {res} })
	names(curRes) <- sub('^.*\\.', '', names(curRes))
	tRes <- tapply(curRes, names(curRes), table)
	resDF <- foreach(res = tRes, .combine=rbind) %do% { res }
	if (length(tRes) == 1) {  #Work around to force a DF
		resDF <- rbind(resDF, resDF)[1,,drop=FALSE]
	}
	rownames(resDF) <- names(tRes)
	resDF <- as.data.frame(resDF)
	resDF$pHEMI <- resDF[, posClass] / rowSums(resDF)
	resDF$prediction <- negClass
	resDF$prediction[resDF$pHEMI > 0.5] <- posClass
	predTally[['All']] <- resDF
	
	actual <- rep('?', nrow(predTally[[1]]))
	actual[grepl(paste('_', posClass, sep=""), rownames(predTally[[1]]))] <- posClass
	actual[grepl(paste('_', negClass, sep=""), rownames(predTally[[1]]))] <- negClass
	
	#TODO: Replace with dedicated function
	#  The function is: calcAccs.  Use it!!!
	for (i in 1:length(predTally)) {
		predTally[[i]]$actual <- actual
		tp <- sum(predTally[[i]]$actual==posClass & predTally[[i]]$prediction==posClass)
		tn <- sum(predTally[[i]]$actual==negClass & predTally[[i]]$prediction==negClass)
		fp <- sum(predTally[[i]]$actual==negClass & predTally[[i]]$prediction==posClass)
		fn <- sum(predTally[[i]]$actual==posClass & predTally[[i]]$prediction==negClass)
		acc <- (tp+tn)/(tp+tn+fp+fn)
		sens <- (tp)/(tp+fn)
		spec <- (tn)/(tn+fp)
		f1 <- (2*tp)/(2*tp+fp+fn)
		rho <- suppressWarnings(cor(predTally[[i]]$actual==posClass, predTally[[i]]$prediction==posClass))
		pVal <- pbinom(round(acc*length(actual)), length(actual), 0.5, lower.tail = FALSE, log.p = FALSE)
		
		if (length(actual) >= 12) {
			scores <- rep("", length(actual))
			scores[1:12] <- c('acc', round(acc,3), 'rho', round(rho,3), 'pVal', format(pVal, digits=3), 'sens', round(sens,3), 'spec', round(spec,3), 'f1', round(f1,3))
			predTally[[i]]$scores <- scores
		}
		
		if ((length(actual) >= 6) && (length(actual) < 12)) {
			scores <- rep("", length(actual))
			scores[1:6] <- c('acc', round(acc,3), 'rho', round(rho,3), 'pVal', format(pVal, digits=3))
			predTally[[i]]$scores <- scores
		}
	}
	
	rv <- list(predTally=predTally, results=results)
	save(rv, file=paste("./intermediate/results.", Sys.Date(), ".pred.RData", sep=""))
	return(rv)
}

#' Create a support vector machine for the data
#'   TODO:  Used best.tune to optimize parameters
#'
#' @param dataMat  The attribute matrix
#' @export
#' @usage svmModel <- runSVM(dataMat, prob.model = FALSE)
runSVM <- function(dataMat, prob.model = FALSE) {
	library(kernlab)
	
	dataPart <- data.frame(apply(dataMat[,1:(ncol(dataMat)-1), drop=FALSE], 2, as.numeric))
	curMat <- cbind(dataPart, data.frame(Class=dataMat[,ncol(dataMat)]))
	
	svmModel <- ksvm(x=as.formula('Class ~ .'), data=curMat, method='class', kernel='polydot', kpar=list(degree=1), shrinking=FALSE, prob.model = prob.model)

	return(svmModel)
}

#' Create a random forest classifier on the data
#'
#' @param dataMat  The attribute matrix
#' @param numTrees  The number of trees (DEAFULT: 1000)
#' @param mtry  # Feature to Try Each Tree.  Try max(floor(ncol(dataMat)/3), 1) (DEAFULT: NULL)
#' @export
#' @usage rfModel <- runRF(dataMat, numTrees=1000, mtry=NULL)
runRF <- function(dataMat, numTrees=1000, mtry=NULL) {
	library(randomForest)
	
	if (is.null(mtry)) {
		rfModel <- randomForest(x = apply(dataMat[,1:(ncol(dataMat)-1), drop=FALSE], 2, as.numeric), y = as.factor(dataMat[,ncol(dataMat)]), ntree=numTrees)
	} else {
		rfModel <- randomForest(x = apply(dataMat[,1:(ncol(dataMat)-1), drop=FALSE], 2, as.numeric), y = as.factor(dataMat[,ncol(dataMat)]), ntree=numTrees, mtry=mtry)
	}

	return(rfModel)
}

#' Cross validate predictions
#'
#' @param dataMatrix  Each row is an example, each col a feature, in order of importance, last col is label
#' @param numsFeatures  The number of top features to select
#' @param folds  The number of folds to cross-validate (DEAFULT: 10)
#' @param verbose  Set to TRUE to display the current cut (DEAFULT: FALSE)
#' @param classType  Classifier Subset to use.  Usually 'SVM' or 'ALL' (DEAFULT: 'ALL')
#' @export
#' @usage xvAcc <- xvalPreds(dataMatrix, numsFeatures=1:100, folds=10, verbose=FALSE, classType='ALL')
xvalPreds <- function(dataMatrix, numsFeatures=1:100, folds=10, verbose=FALSE, classType='ALL') {
	cuts <- cut(1:nrow(dataMatrix), folds)
	randIdxs <- order(runif(nrow(dataMatrix)))
	
	allPreds <- NULL
	for (curCut in unique(cuts)) {
		if (verbose == TRUE) { cat(curCut, ", ") }
		trainIdxs <- randIdxs[cuts!=curCut]
		testIdxs <- randIdxs[cuts==curCut]
		
		curPreds <- makePreds(trainData=dataMatrix[trainIdxs,,drop=FALSE], testData=dataMatrix[testIdxs,,drop=FALSE], numsFeatures=numsFeatures)
		allPreds <- rbind(allPreds, curPreds$predTally[[classType]])
	}
	
	accRow <- calcAccs(allPreds$prediction, allPreds$actual)
	
	return(accRow)
}

#' Estimate the TP, TN, FP, FN, Agreement, F1, MCC, Sensitivity, Specificity, FPR
#'  Note: Assume every relationship not reported in "refRelations" does not exist
#'  Note: This uses an undirect graph
#'
#' @param predRelations  Three Columns: Regulator, Target, Relationship
#' @param refRelations  The reference relationships (usually take from literature)
#' @export
#' @usage accRow <- calcAccs(predRelations, refRelations)
calcAccs <- function(prediction, actual) {
	prediction <- as.character(prediction)
	actual <- as.character(actual)
	
	pos <- prediction[[1]]
	tp <- sum(prediction == pos & actual == pos)
	tn <- sum(prediction != pos & actual != pos)
	fp <- sum(prediction == pos & actual != pos)
	fn <- sum(prediction != pos & actual == pos)
	
	stats <- getF1mcc(tp, fp, tn, fn)
	if (length(as.numeric(prediction==pos)) > 0 & length(as.numeric(actual==pos)) > 0) {
		mcc.p <- cor.test(as.numeric(prediction==pos), as.numeric(actual==pos))$p.value
	} else {
		mcc.p <- 1
	}
	#Sloppy to mix cor.test with the mcc calculation below.  
	
	rv <- cbind(data.frame(tp=tp, fp=fp, tn=tn, fn=fn), as.data.frame(stats), data.frame(mcc.p=mcc.p))
	return(rv)
}

#' Get f1 / mcc and other accuracy measurements
#'  p-value is calculated using binomial distribution (p=0.5) on agreement
#'
#' @param tps  
#' @param fps  
#' @param tns  
#' @param fns  
#'
#' @export
#' @usage pVals <- getF1mcc(tps, fps, tns, fns)
getF1mcc <- function(tps, fps, tns, fns) {
	sensitivity <- tps/(tps+fns)
	specificity <- tns/(fps+tns)
	fpr <- fps / (fps+tns)
	fdr <- fps / (fps+tps)
	f1 <- (2*tps) / (2*tps + fps + fns)
	N <- tns+tps+fns+fps
	S <- (tps+fns) / N
	P <- (tps+fps) / N
	mcc <- (tps/N - S * P) / sqrt(P*S*(1-S)*(1-P))
	if (is.na(mcc) | mcc==Inf) { mcc <- 0 }

	agreement <- (tps+tns)/(tps+fps+tns+fns)
	p.value <- pbinom(q=tps+tns, size=tps+fps+tns+fns, prob=0.5, lower.tail = FALSE)
	
	#TODO: Add correlation based p-value
	
	return(list(sensitivity=sensitivity, specificity=specificity, fpr=fpr, fdr=fdr, f1=f1, mcc=mcc, agreement=agreement, p.value=p.value))
}

#' Convert the data so that each gene has a single expression valu in each timepoint
#'   This is designed for gene level analyse
#'
#' @param origData  The expression data including 'symbol' and 'symbol_description'
#' @export
#' @usage geneData <- exonToGene(origData)
exonToGene <- function (origData) {
	if(is.null(origData$SystematicName)) {
		origData <- cbind(data.frame(SystematicName=origData$symbol), origData)
	}
	origData <- origData[origData$SystematicName!='NA' & !is.na(origData$SystematicName),]
	
	dataIdx <- grep('cell', colnames(origData)) + 1
	outData <- matrix(0, nrow=length(unique(origData$symbol)), ncol=(ncol(origData) - dataIdx + 1))
	rownames(outData) <- unique(origData$symbol)
	colnames(outData) <- colnames(origData)[dataIdx:ncol(origData)]
	
	for (curExp in colnames(outData)) {
		curCol <- tapply(origData[,curExp], origData$symbol, mean, na.rm=TRUE)
		outData[names(curCol), curExp] <- curCol
	}
	rownames(outData) <- make.names(rownames(outData))
	colnames(outData) <- make.names(colnames(outData))
	
	return(outData)
}


#' Calculate the entropy for each experiment in the testSet
#'   H(D) = -p(active|D,M)*log(p(active|D,M)) - p(inactive|D,M)*log(p(inactive|D,M))
#' This calculates the p-Value directly from the hyperplane.
#' 
#' @param trainSet  The training data set to build the classifier
#' @param testSet  The test data to calculate the margins 
#' @param useConsensus  Set to TRUE to use p-Values from consensus preds rather than margin (DEFAULT: FALSE)
#' @return a data fram including the entropy for each example
#' @export
#' @usage dataTable <- getEntropy(trainSet, testSet, useConsensus=FALSE)
getEntropy <- function(trainSet, testSet, useConsensus=FALSE) {
	#classifyDataSets <- makeClassifierDataSets(refName='ALS') #Reference
	# trainSet = classifyDataSets$method2.mProbes.ALS.1
	# testSet = classifyDataSets$method2.mProbes.ALS.new.122012.1
	
	if (useConsensus == FALSE) {
		#Build the model (svm)
		svmModel <- runSVM(trainSet, prob.model = TRUE)
		#predHemi <- predict(mnModel, newdata=curData)
		preds <- predict(svmModel, newdata=testSet, type = "probabilities")
		rownames(preds) <- rownames(testSet)
	
		
	} else { #Use the num predictions to calculate the probabilities
		preds.mp <- makePreds(trainSet, testSet, numsFeatures= 1:(ncol(testSet)-1))
		preds <- preds.mp$predTally$SVM[, 1:2]
		preds <- preds / sum(preds[1,])
	} #if (useConsensus == FALSE) {
	
	#Entropy: 
	entropys <- apply(preds, 1, function(x) {-1*x[1]*log(x[1]) - x[2]*log(x[2])})
	preds <- as.data.frame(preds)
	preds$entropy <- entropys
	preds$entropy[is.nan(preds$entropy)] <- 0
	
	return(preds)
}

#' In depth analysis of the predictions on Huntingtons.
#'
#'   @param alPattern  A pattern to match to fund the correct results form al.59onHunt
#' 		(DEFAULT: '^ActiveLearning.forReal.Huntingtons.NoHEMI.*.2015-04-21\\.RData$')
#'              Preds with fake HEMIs: '^ActiveLearning.forReal.Huntingtons.*2015-04-17\\.RData$'
#'              No HEMIs: '^ActiveLearning.forReal.Huntingtons.NoHEMI.*.2015-04-21\\.RData$'
#'              No HEMIs, 10000 tree feature selection: '^ActiveLearning.forReal.Huntingtons.nTree.*2015-04-24\\.RData$'
#' 
#' @return 
#' @export
#' @usage alResults <- analyze.al.59onHunt(alPattern = '^ActiveLearning.forReal.Huntingtons.NoHEMI.*.2015-04-21\\.RData$')
analyze.al.59onHunt <- function(alPattern = '^ActiveLearning.forReal.Huntingtons.NoHEMI.*.2015-04-21\\.RData$') {
	library(ROCR)

	#Calculate background p(pred = HEMI| HEMI), p(pred = NCAR | NCAR)
	#p.HEMI.HEMI <- .8 
	#p.NCAR.NCAR <- .2
	xvFiles <- list.files(path='./intermediate/', pattern='^xvalResults.032515.*\\.csv$')
	xvStats <- list()
	for (gsaType in c('a1','p1')) {
		xvStats[[gsaType]] <- list()
		xvFile <- xvFiles[grepl(gsaType, xvFiles)]
		xvTable <- read.csv(paste('./intermediate/',xvFile,sep=""))
		for (meth in as.character(xvTable$X)) {
			stats <- xvTable[xvTable$X==meth, c('tp','fp','tn','fn')]
			moreStats <- c(stats['tp']/(stats['tp']+stats['fn']), 
					stats['tn']/(stats['tn']+stats['fp']), 
					stats['fp']/(stats['tn']+stats['fp']), 
					stats['fn']/(stats['tp']+stats['fn']))
			names(moreStats) <- c('p.HEMI.HEMI', 'p.NCAR.NCAR', 'p.HEMI.NCAR', 'p.NCAR.HEMI') #Sensitivity, specificity
			stats <- unlist(c(stats, moreStats))
			xvStats[[gsaType]][[meth]] <- stats
		}
	} #for (gsaType in c('a1','p1')) {
	names(xvStats) <- c('abs1','pos1')
	
	#Analyze the prediction results from al.59onHunt
	alRES <- list()
	alFiles <- list.files(path='.', pattern=alPattern)
	huntPredHEMIs <- list()
	totalHunts <- list()
	hemiPredHEMIs <- list()
	totHemis <- list()
	
	fParts <- strsplit(alFiles[1], '\\.')[[1]]
	fileID <- fParts[length(fParts)-1]
	if(fParts[4] == 'NoHEMI') {
		#TO DO:  Update to record #tree
		fileID <- paste('NoHEMI', fileID, sep='.')
	} else if (fParts[4] == 'nTree') {
		fileID <- paste(fParts[4], fParts[5], fParts[6], fileID, sep='.')
	}
	
	for (alFile in alFiles) {
		gsaType <- 'abs1'
		if (grepl('\\.pos1\\.', alFile)){ gsaType <- 'pos1' }
		meth <- 'method1'
		if (grepl('\\.Method2\\.', alFile)){ meth <- 'method2' }
		if (grepl('\\.mProbes\\.', alFile)){ meth <- paste(meth, 'mProbes', sep='.') }
		meth <- paste(meth, 'svm', sep=".")
		res <- load(alFile)
		allPredsAtOnce <- svmPreds
		alPreds <- alResults.cons$predResults
		
		#This is a bad way to estimate probability from entropy
		pHEMIs <- 1 - alPreds$entropy
		pHEMIs[alPreds$prediction == 'NCAR'] <- 1 - pHEMIs[alPreds$prediction == 'NCAR']
		
		alPreds$pHEMI <- pHEMIs
		alPreds$actual <- alPreds$experiment
		
		#Active learning predictions are probably invalid since we don't know that HEMIs are HEMI
		for (predType in c('ActiveLearning', 'AllAtOnce')) {
			titleStr <- paste(gsaType,meth, predType,sep=".")
			cat(titleStr, '\n')
			if (predType == 'ActiveLearning') {
				svmPreds <- alPreds
			} else { #predType == 'AllAtOnce'
				svmPreds <- allPredsAtOnce
			}
			
		
			#All predictions at once
			#Probability that these NCARs come from the same pool as the HEMIs
			p.HEMI.HEMI <- xvStats[[gsaType]][[meth]][['p.HEMI.HEMI']]
			huntPredHEMI <- sum(svmPreds$actual=='NCAR' & svmPreds$prediction=='HEMI')
			totalHunt <- sum(svmPreds$actual=='NCAR')
			#p.NCAR.HEMI <- xvStats[[gsaType]][[meth]][['p.NCAR.HEMI']]
			#pVal <- pbinom(q=huntPredHEMI, size=totalHunt, prob=p.HEMI.HEMI)
			pVal <- binom.test(huntPredHEMI, totalHunt, p.HEMI.HEMI, 'two.sided')$p.value
			
			#Expected #NCAR predicted Huntingtons
			expectedPredHunt <- totalHunt * xvStats[[gsaType]][[meth]][['p.HEMI.NCAR']]

			hemiPredHEMI <- sum(svmPreds$actual=='HEMI' & svmPreds$prediction=='HEMI')
			totHemi <- sum(svmPreds$actual=='HEMI')
			if (totHemi > 0) {
				pVal.HEMI <- binom.test(hemiPredHEMI, totHemi, p.HEMI.HEMI, 'two.sided')$p.value
				pred <- prediction(svmPreds$pHEMI, svmPreds$actual=='HEMI')
				auc <- performance(pred,"auc")
				AUC <- round(auc@y.values[[1]],3)
			} else {
				pVal.HEMI <- NA
				AUC <- NA
			}
			
			alRES[[titleStr]] <- data.frame(Xval.Sensitivity=p.HEMI.HEMI, huntPredHEMI=huntPredHEMI, totalHunt=totalHunt, pVal=pVal, expectedPredHunt=expectedPredHunt, hemiPredHemi=hemiPredHEMI, totHemi=totHemi, pVal.HEMI=pVal.HEMI, AUC=AUC)
			
			huntPredHEMIs[[titleStr]] <- huntPredHEMI
			totalHunts[[titleStr]] <- totalHunt
			hemiPredHEMIs[[titleStr]] <- hemiPredHEMI
			totHemis[[titleStr]] <- totHemi
			
		} #for (predType in c('ActiveLearning', 'AllAtOnce')) {
		
	} #for (alFile in alFiles) {
	
	#Calculate an aggregate statistic. 
	xvEval <- xvStats
	for (cn in names(xvEval)) {
		xvEval[[cn]] <- xvEval[[cn]][grepl('svm$', names(xvEval[[cn]]))]
	}
	
	tp <- sum(sapply(xvEval, function(x) {sapply(x, function(y) {y['tp']})}))
	fp <- sum(sapply(xvEval, function(x) {sapply(x, function(y) {y['fp']})}))
	tn <- sum(sapply(xvEval, function(x) {sapply(x, function(y) {y['tn']})}))
	fn <- sum(sapply(xvEval, function(x) {sapply(x, function(y) {y['fn']})}))
	p.HEMI.HEMI <- tp / (tp + fn)
	
	huntPredHEMI <- sum(unlist(huntPredHEMIs))
	totalHunt <- sum(unlist(totalHunts))
	hemiPredHEMI <- sum(unlist(hemiPredHEMIs))
	totHemi <- sum(unlist(totHemis))
	
	p.HEMI.NCAR <- fp / (tn + fp)
	expectedPredHunt <- totalHunt * p.HEMI.NCAR
	
	pVal <- binom.test(huntPredHEMI, totalHunt, p.HEMI.HEMI, 'two.sided')$p.value
	if (totHemi > 0) {
		pVal.HEMI <- binom.test(hemiPredHEMI, totHemi, p.HEMI.HEMI, 'two.sided')$p.value
	} else {
		pVal.HEMI <- NA
	}
	alRES[['Aggregate']] <- data.frame(Xval.Sensitivity=p.HEMI.HEMI, huntPredHEMI=huntPredHEMI, totalHunt=totalHunt, pVal=pVal, expectedPredHunt=expectedPredHunt, hemiPredHemi=hemiPredHEMI, totHemi=totHemi, pVal.HEMI=pVal.HEMI, AUC=NA)

	alRES.df <- do.call(rbind.data.frame, alRES)	
	
	write.csv(alRES.df, file=paste("analyze.al.59onHunt", fileID, Sys.Date(), "csv", sep="."))
	return(alRES.df)
}

#' Use 'runActiveLearning' to predict 12 using the 47 example training set
#'   Note 'margin' or 'consensus' should come from the results of testActiveLearning
#'  
#' @return 
#' @export
#' @usage alResults <- testActiveLearning()
al.47on12 <- function() {
	topGSA <- findTopGSA(runNames = c("ALS.47", "ALS.new.053113"))
	classifyDataSets <- makeClassifierDataSets(refName='ALS.47', topGSA=topGSA)
	trainSet <- classifyDataSets$method2.mProbes.ALS.47.1
	testSet <- classifyDataSets$method2.mProbes.ALS.new.053113.1
	
	###!!!###
	alResults.cons <- runActiveLearning(trainSet, testSet, useConsensus=TRUE)
	
	save(alResults.cons, file=paste('ActiveLearning.forReal', Sys.Date(), 'RData', sep="."))
	
	#Now calculate the summary statistics.
	#calcAccs(prediction, actual)
	
	return(list(alResults.cons=alResults.cons))
}

#' Test 'runActiveLearning' using 59 predicts Huntingtons samples to determine if we have specificity
#'   Note: Cross-validation studies imply that method1, ABS(gsa) > 1, SVM is best
#'  
#' @param auto Set auto to 'TRUE' to skip the real active learning run and just use the simulated ones (DEFAULT:FALSE)
#'             Note: with auto on, this will take about a day to run due to the 900 features in method 1
#' @param stripHEMI Set to true to remove the two HEMI examples from Active Learning (DEFAULT:FALSE)
#' 		This is desirable because these samples do not appear to be actually HEMI
#' @param gsaTypes GSA types to include (DEFAULT: c('pos1', 'abs1'))
#' @param methNames Method types to include (DEFAULT: c('Method2.mProbes','Method2','Method1.mProbes','Method1'))
#' @param ntree Number of trees in random forest feature selection (DEFAULT: 1000)
#'
#' @return 
#' @export
#' @usage alResults <- al.59onHunt(auto=FALSE, stripHEMI=FALSE, gsaTypes=c('pos1', 'abs1'), methNames=c('Method2.mProbes','Method2','Method1.mProbes','Method1'), ntree=1000)
al.59onHunt <- function(auto=FALSE, stripHEMI=FALSE, gsaTypes=c('pos1', 'abs1'), methNames=c('Method2.mProbes','Method2','Method1.mProbes','Method1'), ntree=1000) {
	#A common (inner) function to run the classification for all examples
	#Prefix should be 'method1', 'method1.mProbes', 'method2', 'method2.mProbes'
	runClassifications <- function(classifyDataSets, prefix, hnStrs) {
		dsNames <- names(classifyDataSets)[grepl(paste('^', prefix, '.', sep=""), names(classifyDataSets))]
		if (!grepl('mProbes', prefix)) {
			dsNames <- dsNames[!grepl('mProbes',dsNames)]
		}

		trainSet <- classifyDataSets[[dsNames[!grepl('Huntingtons', dsNames)]]]
		testSet <- classifyDataSets[[dsNames[grepl('Huntingtons', dsNames)]]]
		rownames(testSet) <- sub('\\.', '_', rownames(testSet)) #Otherwise naming won't work
		rownames(testSet) <- paste(rownames(testSet), hnStrs, sep='_')

		alResults.cons <- runActiveLearning(trainSet, testSet, useConsensus=TRUE, auto=TRUE)
		allAtOnce <- makePreds(trainSet, testSet, numsFeatures=1:(ncol(trainSet)-1), onlySVM=TRUE)

		return(list(alResults.cons=alResults.cons, allAtOnce=allAtOnce))
	} #runClassifications <- function(classifyDataSets, prefix, hnStrs) {

	autoRes <- list()
	
	if (auto == FALSE) {
		#> tst <- load('preTrial.huntingtons.041015.RData')
		#> tst
		#[1] "trainSet"         "testSet"          "classifyDataSets" "topGSA"
	
		topGSA <- findTopGSA(runNames = c("ALS.59", "ALS.Huntingtons.100114"), absGSAcut = TRUE)
		classifyDataSets <- makeClassifierDataSets(refName='ALS.59', topGSA=topGSA, includeMprobes = FALSE, ntree=ntree) #sortByRefGSA = TRUE) Included during actual hand run
		trainSet <- classifyDataSets$method1.ALS.59.1 
		testSet <- classifyDataSets$method1.ALS.Huntingtons.100114.1

		rownames(testSet) <- sub('\\.', '_', rownames(testSet)) #Otherwise naming won't work

		alResults.cons <- runActiveLearning(trainSet, testSet, useConsensus=TRUE)
		save(alResults.cons, file=paste('ActiveLearning.forReal.Huntingtons', Sys.Date(), 'RData', sep="."))
		autoRes[['Not.Auto']] <- alResults.cons
	} #if (auto == FALSE) {
	
	if (auto == TRUE) {
		filePreamble <- 'ActiveLearning.forReal.Huntingtons'
		
		if (ntree != 1000) {
			filePreamble <- paste(filePreamble, 'nTree', ntree, sep='.')
		}
		
		if (stripHEMI == TRUE) {
			filePreamble  <- paste(filePreamble, 'NoHEMI', sep='.')
			hnStrs <- rep('NCAR', 6)
			#|GSA| > 1
			topGSA <- findTopGSA(runNames = c("ALS.59", "ALS.Huntingtons.100114.noHEMI"))
		} else {
			#04-14-15.  After blind predictions on 4/10/15, reveal names for scripted active learning
			#1 - HEMI, 2 - NCAR, 3 - HEMI, 4 - NCAR, 5 - NCAR, 6 - NCAR, 7 - NCAR, 8 - NCAR
			hnStrs <- rep('NCAR', 8)
			hnStrs[c(1,3)] <- 'HEMI'

			#|GSA| > 1
			topGSA <- findTopGSA(runNames = c("ALS.59", "ALS.Huntingtons.100114"))
		} #if (stripHEMI == TRUE) {
		
		
		for (gsaType in gsaTypes) {
			if (gsaType == 'pos1') { #GSA > 1
				classifyDataSets <- makeClassifierDataSets(refName='ALS.59', topGSA=topGSA, includeMprobes = TRUE, useABScut = FALSE, sortByRefGSA = TRUE, ntree=ntree)
			} else { #|GSA| > 1
				classifyDataSets <- makeClassifierDataSets(refName='ALS.59', topGSA=topGSA, includeMprobes = TRUE, useABScut = TRUE, sortByRefGSA = TRUE, ntree=ntree)
			}
			names(classifyDataSets) <- sub('\\.noHEMI','',names(classifyDataSets))

			for (methName in methNames) {
				methStr <- paste(methName, gsaType, sep='.')
				rcList <- runClassifications(classifyDataSets, chartr('^M', '^m', methName), hnStrs)
				autoRes[[methStr]] <- alResults.cons <- rcList$alResults.cons
				autoRes[[paste(methStr,'allAtOnce',sep='.')]] <- svmPreds <- rcList$allAtOnce$predTally$SVM
				
				fName <- fName.1 <- paste(filePreamble, methStr, Sys.Date(), 'RData', sep=".")
				idx <- 1
				while (file.exists(fName.1)) {
					fName.1 <- sub('.RData', '', fName)
					fName.1 <- paste(fName.1, idx, 'RData', sep=".")
					idx <- idx + 1
				}
				
				save(alResults.cons, svmPreds, file=fName.1)
			} #for (methName c('Method2.mProbes','Method2','Method1.mProbes','Method1')) {
		} #for (gsaType in c('pos1', 'abs1')) {
	} #if (auto ==TRUE) {
	return(autoRes)
}

#' Generate non-active learning predictions for ALS and Huntingtons separately for manuscript 
#'   Note: Cross-validation studies imply that method1, ABS(gsa) > 1, SVM is best
#'   Note: this will take about a day to run due to the 900 features in method 1
#' 
#' @param stripHEMI Set to true to remove the two HEMI examples from Active Learning (DEFAULT:FALSE)
#' @param gsaTypes GSA types to include (DEFAULT: c('pos1', 'abs1'))
#' @param methNames Method types to include (DEFAULT: c('Method2.mProbes','Method2','Method1.mProbes','Method1'))
#' @param ntree Number of trees in random forest feature selection (DEFAULT: 1000)
#'
#' @return 
#' @export
#' @usage fixedResults <- allAtOnce.47on12wHuntv2(tripHEMI=FALSE, gsaTypes=c('pos1', 'abs1'), methNames=c('Method2.mProbes','Method2','Method1.mProbes','Method1'), ntree=1000)
allAtOnce.47on12wHuntv2 <- function(stripHEMI=FALSE, gsaTypes=c('pos1', 'neg1', 'abs1'), methNames=c('Method2.mProbes','Method2','Method1.mProbes','Method1'), ntree=1000) {
	#A common (inner) function to run the classification for all examples
	#Prefix should be 'method1', 'method1.mProbes', 'method2', 'method2.mProbes'
	runClassifications <- function(classifyDataSets, prefix, hnStrs) {
		dsNames <- names(classifyDataSets)[grepl(paste('^', prefix, '.', sep=""), names(classifyDataSets))]
		if (!grepl('mProbes', prefix)) {
			dsNames <- dsNames[!grepl('mProbes',dsNames)]
		}
		
		trainSet <- classifyDataSets[[dsNames[grepl('ALS\\.47', dsNames)]]]
		testSets <- classifyDataSets[dsNames[!grepl('ALS\\.47', dsNames)]]
		testSet <- foreach(ds=names(testSets), .combine = rbind)%do%{testSets[[ds]]}
		cat(paste('trainSet dim:',dim(trainSet), 'testSet dim:', dim(testSet)))
		rownames(testSet) <- sub('\\.', '_', rownames(testSet)) #Otherwise naming won't work
		hnStrs <- rep('NULL', nrow(testSet))
		for (i in 1:nrow(testSet)){
			hnStrs[i] <- LUT[grep(substr(rownames(testSet), 11, nchar(rownames(testSet)))[i], LUT$file.name),'mouse.status']
		}
		
		rownames(testSet) <- paste(rownames(testSet), hnStrs, sep='_')
		
# 		alResults.cons <- runActiveLearning(trainSet, testSet, useConsensus=TRUE, auto=TRUE)
		allAtOnce <- makePreds(trainSet, testSet, numsFeatures=1:(ncol(trainSet)-1), onlySVM=TRUE)
		
		return(list(allAtOnce.cons=allAtOnce))
# 		return(list(alResults.cons=alResults.cons, allAtOnce=allAtOnce))
		
	} #runClassifications <- function(classifyDataSets, prefix, hnStrs) {
	
	LUT <- read.csv(file='~/gitScripts/iCAP/LUTs/ALS/ALS_allAtOnceKey.csv', stringsAsFactors = FALSE)
	autoRes <- list()

		filePreamble <- 'allAtOnce.47on12wHuntv2'
		
		if (ntree != 1000) {
			filePreamble <- paste(filePreamble, 'nTree', ntree, sep='.')
		}
		
		if (stripHEMI == TRUE) {
			filePreamble  <- paste(filePreamble, 'NoHEMI', sep='.')
			hnStrs <- rep('NCAR', 6)
			#|GSA| > 1
# 			topGSA <- findTopGSA(runNames = c("ALS.59", "ALS.Huntingtons.100114"))
			topGSA <- findTopGSA(runNames = c("ALS.47", "ALS.new.053113", "ALS.Huntingtons.100114.noHEMI"))
			
		} else {
			#04-14-15.  After blind predictions on 4/10/15, reveal names for scripted active learning
			#1 - HEMI, 2 - NCAR, 3 - HEMI, 4 - NCAR, 5 - NCAR, 6 - NCAR, 7 - NCAR, 8 - NCAR
			hnStrs <- rep('NCAR', 8)
			hnStrs[c(1,3)] <- 'HEMI'
			
			#|GSA| > 1
# 			topGSA <- findTopGSA(runNames = c("ALS.59",  "ALS.Huntingtons.100114"))
			topGSA <- findTopGSA(runNames = c("ALS.47",  "ALS.new.053113", "ALS.Huntingtons.100114"))
			
		} #if (stripHEMI == TRUE) {
		
		# Make Data Sets
		for (gsaType in gsaTypes) {
			if (gsaType == 'pos1') { #GSA > 1
				classifyDataSets <- makeClassifierDataSets(refName='ALS.47', topGSA=topGSA, includeMprobes = TRUE, useABScut = FALSE, sortByRefGSA = TRUE, ntree=ntree)
			} else if (gsaType == 'abs1') { #|GSA| > 1
				classifyDataSets <- makeClassifierDataSets(refName='ALS.47', topGSA=topGSA, includeMprobes = TRUE, useABScut = TRUE, sortByRefGSA = TRUE, ntree=ntree)
			} else if (gsaType == 'neg1') { #|GSA| > 1
				classifyDataSets <- makeClassifierDataSets(refName='ALS.47', topGSA=topGSA, includeMprobes = TRUE, cutScore = -1, useABScut = FALSE, sortByRefGSA = TRUE, ntree=ntree)
			}
			names(classifyDataSets) <- sub('\\.noHEMI','',names(classifyDataSets))
			
			for (methName in methNames) {
				methStr <- paste(methName, gsaType, sep='.')
				rcList <- runClassifications(classifyDataSets, chartr('^M', '^m', methName), hnStrs)
				autoRes[[methStr]] <- allAtOnce.cons <- rcList$allAtOnce.cons
				autoRes[[paste(methStr,'allAtOnce',sep='.')]] <- svmPreds <- rcList$allAtOnce.cons$predTally$SVM
				
				fName <- fName.1 <- paste(filePreamble, methStr, Sys.Date(), 'RData', sep=".")
				idx <- 1
				while (file.exists(fName.1)) {
					fName.1 <- sub('.RData', '', fName)
					fName.1 <- paste(fName.1, idx, 'RData', sep=".")
					idx <- idx + 1
				}
				
				save(allAtOnce.cons, svmPreds, file=fName.1)
				
# 				save(alResults.cons, svmPreds, file=fName.1)
			} #for (methName c('Method2.mProbes','Method2','Method1.mProbes','Method1')) {
		} #for (gsaType in c('pos1', 'abs1')) {

	return(autoRes)
}

#' predict the 12 samples used for blind predictions all at once, include Huntingtons samples
#' 	
#'   Based on al.47on12
#'  
#' 
#' @return 
#' @export
#' @usage allAtOnceResults <- allAtOnce.47on12wHunt()
allAtOnce.47on12wHunt <- function() {
	topGSA <- findTopGSA(runNames = c("ALS.47", "ALS.new.053113", "ALS.Huntingtons.100114"))
	classifyDataSets <- makeClassifierDataSets(refName='ALS.47', topGSA=topGSA)
	trainSet <- classifyDataSets$method2.mProbes.ALS.47.1 #samples in rows and features in columns
	testSet <- rbind(classifyDataSets$method2.mProbes.ALS.new.053113.1, classifyDataSets$method2.mProbes.ALS.Huntingtons.100114.1)
	cat(paste('dim trainSet:', dim(trainSet), 'dim testSet:', dim(testSet)),sep='\n')
  browser()
	allAtOnce.cons <- makePreds(trainSet, testSet, numsFeatures=1:(ncol(trainSet)-1), onlySVM=TRUE)

	save(allAtOnce.cons, file=paste('allAtOnce', Sys.Date(), 'RData', sep="."))
	
	#Now run again using an absolute value GSA cutoff
	classifyDataSets <- makeClassifierDataSets(refName='ALS.47', topGSA=topGSA, useABScut = TRUE)
	trainSet <- classifyDataSets$method2.mProbes.ALS.47.1
	testSet <- rbind(classifyDataSets$method2.mProbes.ALS.new.053113.1, classifyDataSets$method2.mProbes.ALS.Huntingtons.100114.1)
	
	allAtOnce.cons.absGSA <- makePreds(trainSet, testSet, numsFeatures=1:(ncol(trainSet)-1), onlySVM=TRUE)

	save(allAtOnce.cons.absGSA, file=paste('allAtOnce.absGSA', Sys.Date(), 'RData', sep="."))
	
	#Now run again using a GSA score <= -1
	classifyDataSets <- makeClassifierDataSets(refName='ALS.47', topGSA=topGSA, cutScore=-1, useABScut = FALSE)
	trainSet <- classifyDataSets$'method2.mProbes.ALS.47.-1'
	testSet <- rbind(classifyDataSets$'method2.mProbes.ALS.new.053113.-1', classifyDataSets$'method2.mProbes.ALS.Huntingtons.100114.-1')
	
	allAtOnce.cons.negGSA <- makePreds(trainSet, testSet, numsFeatures=1:(ncol(trainSet)-1), onlySVM=TRUE)

	save(allAtOnce.cons.negGSA, file=paste('allAtOnce.negGSA', Sys.Date(), 'RData', sep="."))
	
	#Now run again using GSA score > 1 and method 1
	classifyDataSets <- makeClassifierDataSets(refName='ALS.47', topGSA=topGSA, cutScore=1, useABScut = FALSE)
	##!!##
	trainSet <- classifyDataSets$'method1.mProbes.ALS.47.1'
	testSet <- rbind(classifyDataSets$'method1.mProbes.ALS.new.053113.1', classifyDataSets$'method1.mProbes.ALS.Huntingtons.100114.1')
	
	allAtOnce.cons.posGSA.method1 <- makePreds(trainSet, testSet, numsFeatures=1:(ncol(trainSet)-1), onlySVM=TRUE)
	save(trainSet, testSet, allAtOnce.cons.posGSA.method1, file=paste('allAtOnce.posGSA.method1', Sys.Date(), 'RData', sep="."))
	
	
	#Now run again using abs(GSA) score > 1 and method 1
	classifyDataSets <- makeClassifierDataSets(refName='ALS.47', topGSA=topGSA, cutScore=1, useABScut = TRUE)
	trainSet <- classifyDataSets$'method1.mProbes.ALS.47.1'
	testSet <- rbind(classifyDataSets$'method1.mProbes.ALS.new.053113.1', classifyDataSets$'method1.mProbes.ALS.Huntingtons.100114.1')
	
	allAtOnce.cons.absGSA.method1 <- makePreds(trainSet, testSet, numsFeatures=1:(ncol(trainSet)-1), onlySVM=TRUE)
# 	alResults.cons.absGSA.method1 <- runActiveLearning(trainSet, testSet, useConsensus=TRUE)
	save(trainSet, testSet, allAtOnce.cons.absGSA.method1, file=paste('allAtOnce.absGSA.method1', Sys.Date(), 'RData', sep="."))
	
	return(list(allAtOnce.cons=allAtOnce.cons, allAtOnce.cons.absGSA=allAtOnce.cons.absGSA, 
							allAtOnce.cons.negGSA=allAtOnce.cons.negGSA, allAtOnce.cons.posGSA.method1=allAtOnce.cons.posGSA.method1,
							allAtOnce.cons.absGSA.method1=allAtOnce.cons.absGSA.method1))
}

#' Given a training set and a test set,
#'   Iteratively chose the minimum entropy example and have the user
#'   interactively reveal the actual class
#'   Then retrain the classifier, and pick the next one.
#' 
#' @param trainSet  The training data set to build the classifier
#' @param testSet  The test data to calculate the margins 
#' @param useConsensus  Either 'margin' or 'consensus' (default: TRUE, i.e. consensus)
#'                    Note, predictions should always be consensus, they're better
#' @param auto  Set to TRUE to read actual class from experiment names (DEFAULT: FALSE)
#' @return 
#' @export
#' @usage alResults <- runActiveLearning(trainSet, testSet, useConsensus=FALSE, auto=FALSE)
runActiveLearning <- function(trainSet, testSet, useConsensus=FALSE, auto=FALSE) {

	trainSet.orig <- trainSet
	testSet.orig <- testSet

	numIters <- nrow(testSet)
	predResults <- NULL #a data frame with expNames, predicted, experiment
	
	for (i in 1:numIters) {
		preds.con <- preds.entropy <- getEntropy(trainSet, testSet, useConsensus=FALSE)
		preds.con <- getEntropy(trainSet, testSet, useConsensus=TRUE)
		 # comment out second like as a hack to speed up debugging
		cors <- try(cor.test(preds.con[,3], preds.entropy[,3]), silent=TRUE)
		if (class(cors) == "try-error") { cors <- list(estimate=NaN) }

		#Note, use other entropy method to break ties
		if (useConsensus == TRUE) {
			preds.use <- preds.con
			preds.xcheck <- preds.entropy
		} else {
			preds.use <- preds.entropy
			preds.xcheck <- preds.con
		}
		minVal <- min(preds.use$entropy)
		minIdxs <- which(preds.use$entropy == minVal) 
		if (length(minIdxs) > 1) {
			#Cross-check
			preds.sub <- preds.xcheck[minIdxs,]
			chosenExample <- rownames(preds.sub)[which.min(preds.sub$entropy)]
		} else {
			chosenExample <- rownames(preds.use)[minIdxs]
		}
		chosenEntropy <- preds.use[chosenExample, 'entropy']
		pHEMI <- preds.con[chosenExample, 'HEMI']
		prediction <- 'NCAR'
		if (pHEMI >= 0.5) { prediction <- 'HEMI' } #1 less HEMI example than NCAR, so bias towards HEMI
		
		cat('\nMinimum Entropy Example: ', chosenExample, '\n')
		cat('Entropy:', chosenEntropy, '\n', 'Prediction:', prediction, '\n')
		
		contin <- TRUE
		if (auto == TRUE) {
			class <- '?'
			if( grepl('_HEMI', chosenExample) ) { class <- 'HEMI' }
			if( grepl('_NCAR', chosenExample) ) { class <- 'NCAR' }
			if ( class != '?' ) {
				contin <- FALSE
			}
		} #if (auto == TRUE) {
		
		while(contin == TRUE) {
			answer <- readline(paste("Is ", chosenExample, " HEMI or NCAR? ", sep=""))
			if(grepl('^h', tolower(answer))) {class <- 'HEMI'} else {class <- 'NCAR'}

			answer <- readline(paste("Is ", chosenExample, ' ', class, ' (YES / NO)? ', sep=""))
			if(grepl('^y', tolower(answer))) {contin <- FALSE}
		} #while(contin == TRUE) {
		newRow <- data.frame(example=chosenExample, entropy=chosenEntropy, entropy.rho=cors$estimate, prediction=prediction, experiment=class)
		if (is.null(predResults)) {
			predResults <- newRow 
		} else {
			predResults <- rbind(predResults, newRow)
		}
		
		#Move the now classified example from the test the training set
		newRow <- testSet[chosenExample,, drop=FALSE]
		newRow$class <- class
		
		trainSet <- rbind(trainSet, newRow)
		ntrain <- nrow(trainSet)
		testSet <- testSet[rownames(testSet) != chosenExample, , drop=FALSE]
		ntest <- nrow(testSet)
		fName <- paste('./intermediate/ActiveLearning.train', ntrain, 'test', ntest, Sys.Date(), 'RData', sep=".")
		save(trainSet, testSet, file=fName)
		
		cat('Next training run will have:', ntrain, 'training examples and', ntest, 'testing examples.\n\n')
	} #for (i in 1:numIters) {
	
	#Add an accuracy row
	accRow <- calcAccs(prediction=predResults$prediction, actual=predResults$experiment)
	
	fName <- paste('./intermediate/ActiveLearning.results', Sys.Date(), 'RData', sep=".")
	save(predResults, accRow, file=fName)	
	return(list(predResults = predResults, accRow = accRow))
}

#' Convert an entropy value into possible p-Values
#'    entropys <- apply(preds, 1, function(x) {-1*x[1]*log(x[1]) - x[2]*log(x[2])})
#'    This is where x[2] = 1 - x[1]
#' 
#' @param entropy  The entropy value(s) (DEFAULT: NULL, i.e. select 3 at random)
#' @param plotIt  Set to true to plot the entropy curve (DEFAULT: FALSE)
#' @param numSamples  The number of samples to determine resolution (DEFAULT: 10000)
#' @return a data of microarray data
#' @export
#' @usage pVal <- reverseEntropy(entropy=NULL, plotIt=FALSE, numSamples=10000)
reverseEntropy <- function(entropy=NULL, plotIt=FALSE, numSamples=10000) {
	#Make the curve
	numSamples <- 10000
	pVals <- (1:numSamples) / (2*numSamples)
	entropys.lte.5 <- sapply(pVals, function(x) {-1*x*log(x) - (1-x)*log(1-x)})
	
	if (is.null(entropy)) {
		idxs <- sample(1:length(entropys.lte.5), 3)
		entropy <- entropys.lte.5[idxs]
		testPs <- pVals[idxs]
	}
	
	#Find the minimum differences
	idxs <- sapply(entropy, function(x) {which.min(abs(entropys.lte.5-x))})
	newPs <- pVals[idxs]
	
	if (plotIt == TRUE) {
		plot(pVals, entropys)
		title('Plot of p-Values from Entropy')
		points(newPs, entropys[idxs], pch=4, col='red', cex=3, lwd=3)
	}
	return(newPs)
}

#' Reverse engineer AL prediction results to generate AUC
#' 
#' @param entropy  The entropy value(s) (DEFAULT: NULL, i.e. select 3 at random)
#' @return a data of microarray data
#' @export
#' @usage pVal <- recreateAUC(fileName = './intermediate/ActiveLearning.results.2015-01-12.RData')
recreateAUC <- function(fileName = './intermediate/ActiveLearning.results.2015-01-12.RData') {
	varNames <- load(fileName)
	if (any(grepl('^predRes', varNames))) {
		predResults <- get(varNames[grepl('^predRes', varNames)])
	} else if (any(grepl('^alResults.cons', varNames))) {
		alResults.cons <- get(varNames[grepl('^alResults.cons', varNames)])
		predResults <- alResults.cons$predResults
	}
	pHEMI <- reverseEntropy(entropy=predResults$entropy)
	
	pHEMI[predResults$prediction == 'HEMI']  <- 1 - pHEMI[predResults$prediction == 'HEMI']
	
	###!!!###  NOW FIND THE CORRECT FILES
	pdfFile <- sub('\\.RData', '.AUC.pdf', fileName)
	pdf(pdfFile)
	titleStr <- sub('^.*\\.', '', sub('\\.AUC\\.pdf', '', pdfFile))
	#predStats <- unlist(getAUCandP(pClass=1-pHEMI, isClass=predResults$experiment=='NCAR', titleStr=titleStr))
	predStats <- unlist(getAUCandP(pClass=pHEMI, isClass=predResults$experiment=='HEMI', titleStr=titleStr))
	dev.off()
	
	#Resort the order
	curNames <- names(predStats)
	mccBool <- curNames=='mcc'
	curNames <- curNames[!mccBool]
	mccIdx.P <- which(curNames=='mcc.p')
	curNames <- c(curNames[1:(mccIdx.P-1)], 'mcc', curNames[mccIdx.P:length(curNames)])

	predStats <- predStats[curNames]
	return(predStats)
}

#' Reverse engineer AL prediction results to generate AUC
#' 
#' @param entropy  The entropy value(s) (DEFAULT: NULL, i.e. select 3 at random)
#' @return a data of microarray data
#' @export
#' @usage pValList <- findAllAUCs(pattern='^ActiveLearning\\.results.*\\.RData$')
findAllAUCs <- function(pattern='^ActiveLearning\\.forReal.*\\.RData$') {
	#fNames <- list.files(path='./intermediate', pattern=pattern, full.names = TRUE)
	fNames <- list.files(path='.', pattern=pattern, full.names = TRUE)
	
	pValList <- list()
	for (fName in fNames) {
		pValList[[fName]] <- recreateAUC(fileName = fName)
	}
	
	#Notes:
	#Classifier 1 = 2014-12-11
	#Internal Name: Method 1, GSA > 1
	#         auc    auc.pVal
	# 0.875000000 0.011415250
	#Huntington's Name: ActiveLearning.forReal.Huntingtons.NoHEMI.Method1.mProbes.pos1.2015-04-21
	predFile <- "./ActiveLearning.forReal.2014-12-11.RData"
	load(predFile)
	predRes <- alResults.cons$predResults
	
	hunt.predFile <- "./ActiveLearning.forReal.Huntingtons.NoHEMI.Method1.mProbes.pos1.2015-04-21.RData"
	load(hunt.predFile)
	predRes.hunt <- alResults.cons$predResults
	
	predRes <- rbind(predRes, predRes.hunt)
	fName <- 'Classifier1.PredAndHunt.RData'
	save(predRes, file=fName)
	accRow <- recreateAUC(fileName = fName)
	save(predRes, accRow, file=fName)
	pValList[['Classifier1.PredAndHunt']] <- accRow
	
	#
	#Classifier 2 = 2015-01-12
	#Internal Name: Method 2, GSA > 1
	#         auc    auc.pVal
	# 0.9027777778 0.0027059339
	#Huntington's Name: ActiveLearning.forReal.Huntingtons.NoHEMI.Method2.mProbes.pos1.2015-04-21
	predFile <- "./ActiveLearning.forReal.posGSA.method1.2015-01-12.RData"
	#Note, the name 'method1' was named by hand and is a misnomer
	tst <- load(predFile)
	predRes <- alResults.cons.posGSA.method1$predResults

	hunt.predFile <- "./ActiveLearning.forReal.Huntingtons.NoHEMI.Method2.mProbes.pos1.2015-04-21.RData"
	tst <- load(hunt.predFile)
	predRes.hunt <- alResults.cons$predResults

	predRes <- rbind(predRes, predRes.hunt)
	fName <- 'Classifier2.PredAndHunt.RData'
	save(predRes, file=fName)
	accRow <- recreateAUC(fileName = fName)
	save(predRes, accRow, file=fName)
	pValList[['Classifier2.PredAndHunt']] <- accRow
	
	
	#
	#Classifier 3 = 2015-01-08
	#Internal Name: Method 1, |GSA| > 1
	#         auc    auc.pVal
	# 0.888888889 0.011210211
	#Huntington's Name: ActiveLearning.forReal.Huntingtons.NoHEMI.Method1.mProbes.abs1.2015-04-24 (is 4-21 or 4-24 correct?)
	predFile <- "./ActiveLearning.forReal.absGSA.2015-01-08.RData"
	tst <- load(predFile)
	predRes <- alResults.cons.absGSA$predResults

	hunt.predFile <- "./ActiveLearning.forReal.Huntingtons.NoHEMI.Method1.mProbes.abs1.2015-04-24.RData"
	tst <- load(hunt.predFile)
	predRes.hunt <- alResults.cons$predResults

	predRes <- rbind(predRes, predRes.hunt)
	fName <- 'Classifier3.PredAndHunt.RData'
	save(predRes, file=fName)
	accRow <- recreateAUC(fileName = fName)
	save(predRes, accRow, file=fName)
	pValList[['Classifier3.PredAndHunt']] <- accRow

	#1) Agglomerate ALS predictions with huntingtons predictions.
	
	#To Do:
	#2) Make all-at-once predictions
	outFile <- paste('findAllAUCs', Sys.Date(), 'csv', sep='.')
	write.csv(pValList, file=outFile)
	return(pValList)
}

#' Perform a cross-validation of Karan's features
#'
#' @param pClass  The probability of being in 'class'
#' @param isClass  A boolean vect for if the example is actually 'class'
#' @export 
#' @usage preds <- getAUCandP(pClass, isClass, titleStr=NULL)
getAUCandP <- function(pClass, isClass, titleStr=NULL) {

	tps <- sum((pClass >= 0.5) & isClass)
	tns <- sum((pClass < 0.5) & !isClass)
	fps <- sum((pClass >= 0.5) & !isClass)
	fns <- sum((pClass < 0.5) & isClass)
	pVals <- getF1mcc(tps, fps, tns, fns)
	
	pVals$mcc.p <- cor.test(as.numeric(pClass >= 0.5), as.numeric(isClass))$p.value
	
	#Can actually use wilcox for everything
	#x1 = pClass[isClass==TRUE]; 
	#n1 = length(x1);            # prepare input data ...
	#x2 = pClass[isClass==FALSE]; 
	#n2 = length(x2);
	#wt <- wilcox.test(x1, x2, exact=0)
	#auc <- wt$statistic / (n1*n2)
	if (!all(isClass) & !all(!isClass)) {	
		rocStats <- suppressWarnings(roc.area(as.numeric(isClass), pClass))
		pVals$auc <- rocStats$A
		#pVals$auc.pVal <-  rocStats$p.value

		#Regression based p-value
		#http://stats.stackexchange.com/questions/75050/in-r-how-to-compute-the-p-value-for-area-under-roc
		diagnosticData <- data.frame(Group=isClass, continuousVar=pClass)
		GLM.1 <- glm(Group ~ continuousVar, family=binomial(logit), data=diagnosticData)
		GLM.2 <- glm(Group ~ 1, family=binomial(logit), data=diagnosticData) 
		at <- anova(GLM.2, GLM.1, test="Chisq")
		pVals$auc.pVal <- at[["Pr(>Chi)"]][2]

		if (!is.null(titleStr)) {
			curTitle <- paste('ALS Prediction (ANOVA p-Value ', signif(pVals$auc.pVal, digits = 3), ')', sep='')
			dataList <- list()
			dataList[[titleStr]] <- data.frame(predictions=pClass, labels=isClass)
			handel <- ROCon (dataList=dataList, titleStr=paste(curTitle, titleStr, sep='\n'), cex=1.0, lwd=4)
		}
	} else {
		pVals$auc <- as.numeric(NA)
		pVals$auc.pVal <- as.numeric(NA)
	}
	
	return(pVals)
}

#' Plot a list of data onto a ROC curve
#'
#' @param dataList A named list.  Each element of the list has $predictions and $labels
#' @param titleStr The title for the plot
#' @param cex The multiplier for the text size.  (DEFAULT: 1)
#' @param lwd The line width.  (DEFAULT: 4)
#' @export
#' @usage handel <- ROCon (dataList, titleStr, cex=1.0, lwd=4) 
ROCon <- function(dataList, titleStr, cex=1.0, lwd=4) {
	if (length(dataList) <= 3) {
		colVect <- c('red','blue','darkgreen')[1:length(dataList)]
	} else {
		colVect <- rainbow(length(dataList))
	} #This code section jsut picks out the color palette 

	ltys = rep(c(1,2,4,5,6), ceiling(length(dataList)/5))

	legendNames <- NULL
	for (i in 1:length(dataList)) {
		p1 <- dataList[[i]]$predictions
		l1 <- dataList[[i]]$labels
		keepBool <- !is.na(p1) & !is.na(l1)
		p1 <- p1[keepBool]
		l1 <- l1[keepBool]
		pred <- prediction(p1, l1)
		perf <- performance(pred, measure = "tpr", x.measure = "fpr")
		perf@'x.name' <- 'Sensitivity' 
		perf@'y.name' <- '1-Specificity'
		if (i == 1) {
			handle <- plot(perf, main=titleStr, col=colVect[i], lwd=lwd, cex.lab=cex, cex.axis=cex, cex.main=cex, cex.sub=cex, lty=ltys[i])
			#axis(1, cex.axis=cex)
			#axis(2, cex.axis=cex)
			#Note: Due to a bug in ROCR, cex.axis=cex doesn't seem to work, liekwise axis.false=FALSE
			abline(a=0, b=1) #slope=1, intercept=0
		} else {
			plot(perf, col=colVect[i], add=TRUE, lwd=lwd, lty=ltys[i])
		}

		auc <- performance(pred,"auc")
		legendName <- paste(names(dataList)[i], ' (AUC = ', sprintf('%.3f', round(auc@y.values[[1]],3)), ')', sep="")
		legendNames <- c(legendNames, legendName)
	}

	legend('bottomright', legend=legendNames, col=colVect, lty=ltys, lwd=lwd, cex=cex)

	return (handle)
} #End ROCon


#' Normalizes given test data to given reference data
#'   Author Karan Singh
#'
#' Creates a normalized and filtered version of testData matrix 
#'	whose vector distributions match those of the referenceData 
#'	vectors (training set) used to construct the classifier. 
#'    This method also removes columns in testData matrix that 
#'	are not in referenceData
#'
#' @param referenceData Reference set to be normalized against
#' @param testData Test set to be normalized
#'  
#' @return Returns normalized and filtered data set
#' @examples
#' normal <- normalizeTestData (testSet, trainingData)
normalizeTestData <- function(tData, referenceData) {	
	subset <- colnames(tData) %in% colnames(referenceData)
	tData <- tData[,subset]

	for (col in colnames(tData)) {
		normalCol <- normalize.quantiles.use.target(as.matrix(tData[,col]), referenceData[,col], copy=TRUE)
		tData[,col] <- normalCol
	}
	
	#tData <- foreach (col = colnames(tData), .combine=cbind) %dopar% {
	#	normalCol <- normalize.quantiles.use.target(as.matrix(tData[,col]), referenceData[,col], copy=TRUE)
	#	normalCol
	#}

	sameVal <- sapply(colnames(tData), function(x) length(unique(tData[,x])) < 2)
	tData <- tData[,!sameVal]
	return(tData)
}


#' Generate the list of differentially expressed genes
#'    Run for final analysis 12/17/14
#' @param 
#'
#' @export 
#' @usage makecMonkeyData() 
genDiffExpr <- function(dataSetName='ALS.47'){
	#Load the data   
	fsScores <- getFIRMAandGENE(dataSetName)
	
	#Differentially expressed	
	#pVals <- compareProbeCols(origData, testCols, controlCols, multProbeMethod=4, wilcox=TRUE)
	
	#Use mProbes selection
	fsScores.analyzed.exons <- analyseFSscores(fsScores = fsScores$geneExpr.vals, curPs=NULL, catagories=c('HEMI','NCAR'), pCut=0.05, numFeat=800, useWilcox=TRUE, method=3, minFeatures=100)
	fName = paste('./intermediate/', dataSetName, '.diffExon.mProbes.', Sys.Date(), '.csv', sep='')
	write.csv(fsScores.analyzed.exons, file=fName)
	
	geneExpr.vals <- fsScores$geneExpr.vals
	geneEx <- exonToGene(geneExpr.vals)
	curSymbols <- make.names(geneExpr.vals$symbol, unique=TRUE)
	
	geneExpr.vals$symbol <- as.character(geneExpr.vals$symbol)
	geneExpr.vals$symbol[grepl('^[0-9].*Rik$', geneExpr.vals$symbol)] <- sub('^', 'X', geneExpr.vals$symbol[grepl('^[0-9].*Rik$', geneExpr.vals$symbol)])
	gns <- rownames(geneEx)
	curSymbols <- make.names(geneExpr.vals$symbol, unique=TRUE)
	idxs <- match(gns, curSymbols)
	newHead <- geneExpr.vals[idxs ,1:grep('^cell$', colnames(geneExpr.vals))]
	geneEx <- cbind(newHead, geneEx)
	rownames(geneEx) <- gns
	#Add headers and so on back into the geneEx
	
	fsScores.analyzed.genes <- analyseFSscores(fsScores = geneEx, curPs=NULL, catagories=c('HEMI','NCAR'), pCut=0.05, numFeat=800, useWilcox=TRUE, method=3, minFeatures=100)
	fName = paste('./intermediate/', dataSetName, '.diffGene.mProbes.', Sys.Date(), '.csv', sep='')
	write.csv(fsScores.analyzed.genes, file=fName)
	
	return(list(fsScores.analyzed.exons=fsScores.analyzed.exons, fsScores.analyzed.genes=fsScores.analyzed.genes))
}


#' Load all of the cross-validation files in './intermediate', e.g. classifyDataSets.ALS.47.1.RData
#'   To compute some summary statistics on the pathways chosen.
#'
#' @param prefix  The file prefix (DEFAULT: 'ALS.47')
#' @param postFix  '' for GSA => 1, '-1' for GSA <= -1, 'useABScut' = useABScut (DEFAULT: '')
#' @export
#' @return  a list containing all of the mProbes lists that have been selected as well as summary statistics
#' @usage resultsList <- find.mProbes.range(prefix='ALS.47', postFix="")
find.mProbes.range <- function(prefix='ALS.47', postFix=""){
	mpFiles <- list.files(path = "./intermediate", pattern = "^classifyDataSets.ALS.47.*RData$")
	tags <- strsplit(mpFiles, '\\.')
	tags <- sapply(tags, function(x) {x[which(x=='RData')-1]})
	negs <- which(tags == "-1")
	abss <- which(tags == "useABScut")
	poss <- which(tags != "-1" & tags != "useABScut")
	if(postFix == "") {
		mpFiles <- mpFiles[poss]
	} else if(postFix == "-1") {
		mpFiles <- mpFiles[negs]
	} else if(postFix == "useABScut") {
		mpFiles <- mpFiles[abss]
	} #if(postFix == "") {
	
	expData <- getFIRMAandGENE('ALS.47')
	geneNames <- unique(expData$geneExpr.vals$symbol)
	genesetList <- getGenesets(geneNames=geneNames)
	
	summaryStats <- NULL
	summaryStats.genes <- NULL
	for (mpFile in mpFiles) {
		cat(mpFile, '\n')
		classifyDataSets <- get(load(paste('./intermediate/', mpFile, sep=""))[1])
		classifyDataSets <- classifyDataSets [grepl('train', names(classifyDataSets))]
		colCounts <- sapply(classifyDataSets, ncol)
		summaryStats <- rbind(colCounts, summaryStats)
		
		#geneSets.mProbes <- colnames(classifyDataSets$method2.mProbes.train.ALS.47.8)
		#geneSets.mProbes <- geneSets.mProbes[geneSets.mProbes != "class"]
		#genes.gs.mProbes <- colnames(classifyDataSets$method1.mProbes.train.ALS.47.8)
		#genes.gs.mProbes <- genes.gs.mProbes[genes.gs.mProbes != "class"]
		
		
		#Also look up the genes in the pathways selected by mProbes
		colCounts.genes <- colCounts
		for(expName in names(colCounts)) {
			if (grepl('^method2', expName)) {
				cat('\t',expName,'\n')
				geneSets <- colnames(classifyDataSets[[expName]])
				uGenes <- geneSetsToGenes(geneSets=geneSets, geneNames=geneNames, genesetList=genesetList)
				colCounts.genes[expName] <- length(uGenes)
			} #if (grepl('^method2', expName)) {
		} #for(expName in names(colCounts)) {
		summaryStats.genes <- rbind(colCounts.genes, summaryStats.genes)
	}
	return(list(summaryStats=summaryStats, summaryStats.genes=summaryStats.genes))
}

#' Return the genes in genesets
#'   Note, this will not take into acount genes that may have been filtered out
#'   when generating the genesets
#'
#' @param geneSets  
#' @param geneNames The list of possible gene names (DEFAULT: NULL, i.e. all genes ALS.47 )
#' @param genesetList The list of gene sets (DEFAULT: NULL, i.e. getGenesets(geneNames=geneNames) )
#'
#' @export
#' @usage uGenes <- geneSetsToGenes(geneSets, geneNames=NULL)
geneSetsToGenes <- function(geneSets, geneNames=NULL, genesetList=NULL) {
	if (is.null(geneNames)) {
		expData <- getFIRMAandGENE('ALS.47')
		geneNames <- unique(expData$geneExpr.vals$symbol)
	}
	
	if (is.null(genesetList)) {
		genesetList <- getGenesets(geneNames=geneNames)
	}
	curSets <- geneSets[geneSets %in% names(genesetList)]  #This should exclude 'class'
	geneList <- genesetList[curSets]
	
	uGenes <- unique(unlist(geneList))
	return(uGenes)
}


#' @param imFile The file produced by extract_inclusion_matrix.py (DEFAULT: '../BSCM/cmonkey-python/out.ALS.0.05/inclusionMatrix.031615.tsv' )
#'
#' @export
#' @usage uGenes <- testClusterClassifier(imFile = '../BSCM/cmonkey-python/out.ALS.0.05/inclusionMatrix.031615.tsv')
testClusterClassifier <- function(imFile = '../BSCM/cmonkey-python/out.ALS.0.05/inclusionMatrix.031615.tsv') {
	imatrix <- read.table(imFile, row.names=1, header=TRUE)
	classVect = as.character(grepl('HEMI', colnames(imatrix)))
	nmatrix <- matrix(0, nrow=nrow(imatrix), ncol=ncol(imatrix))
	nmatrix[imatrix == 'UP'] <- 1
	nmatrix[imatrix == 'DOWN'] <- -1
	rownames(nmatrix) <- rownames(imatrix)
	colnames(nmatrix) <- colnames(imatrix)
	
	nmatrix <- cbind(t(nmatrix), data.frame(Class=classVect))
	
	#Cross-validate
	#xvAcc <- xvalPreds(nmatrix, numsFeatures=ncol(nmatrix)-1, folds=10, verbose=FALSE, classType='SVM')
	#Need to rewrite xvAcc to work for a matrix.
	svm <- ksvm(nmatrix)  #xvalidated error = 0.249177 i.e. 75% acc.
	
	#mProbes
	
	#Cross-validate
	
	return(uGenes)
}

#' Look in the cross-validations to see how often give example(s) are predicted correctly.
#'
#' @param examples  The examples to investigate (DEFAULT: c('HB9EBs_HEMI123_9wk_24h', 'HB9EBs_HEMI124_9wk_24h') )
#'
#' @export
#' @usage predDF <- findAccuracyOnExamples(examples=c('HB9EBs_HEMI123_9wk_24h', 'HB9EBs_HEMI124_9wk_24h'))
findAccuracyOnExamples <- function(examples=c('HB9EBs_HEMI123_9wk_24h', 'HB9EBs_HEMI124_9wk_24h')) {
	xvFiles <- list.files(path='./intermediate/', pattern='^xvalResults.*\\.RData$')
	xvRunList <- list()	
	for (xvFile in xvFiles) {
		curXV <- load(paste('./intermediate/', xvFile, sep=''))
		for (meth in names(resultLists.svm)) {
			curPreds <- resultLists.svm[[meth]]
			idxlist <- unlist(lapply(examples, function(x) {grep(x, rownames(curPreds))}))
			if (length(idxlist) > 0) {
				xvRunList[[paste(xvFile, meth, sep='.')]] <- curPreds[idxlist,]
			}
		} #for (meth in name()) {
	} #for (xvFile in xvFiles) {
	
	#Generate an easy to read data frame
	curList <- list()
	for (curEx in examples) {
		curList[[curEx]] <- lapply(xvRunList, function(x) {x[grepl(curEx, rownames(x)), 'prediction']})
	}
	curDF <- t(do.call(rbind.data.frame, curList))
	curDF <- rbind(curDF, apply(predDF, 2, function(x) {mean(x=='HEMI')}))
	rownames(curDF)[nrow(curDF)] <- 'Mean.HEMI'
	
	write.csv(curDF, file=paste( 'xvPredsOn', Sys.Date(), 'csv', sep='.'))
	return(curDF)
}