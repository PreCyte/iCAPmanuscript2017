#'#'#'#'#'#'#'#'#'#'#'#'#'#'#'#'#'#'#'#'#'#'#'#'#'#'#'#'#'#'#'#'#'#'#'#
# Code used to develop Alzhiemer's Proof of concept classifier       	#
#   NOTE: Piggyback on development in ALSClassifier.R 								#
#                                                     								#
# @author Sam Danziger, PhD, G. Adam Whitney, PhD                     #
# @institution ISB / Seattle Biomed                 									#
#'#'#'#'#'#'#'#'#'#'#'#'#'#'#'#'#'#'#'#'#'#'#'#'#'#'#'#'#'#'#'#'#'#'#'#

#NOTE: This was cleaned up from ALZ.R in ~/gitScripts/iCAP/Library/Legacy/ALZ.R

source('ALSClassifier.R')
library(randomForest)
library(verification)
library(ROCR)
library(caret)
library(kernlab)
library(irr)
library(foreach)

doParallel::registerDoParallel(parallel:::detectCores())

#' Run FIRMA pipeline 
#'  NOTE::  GSA algorithm 'Restandardizes' rows are part of the perumation test.
#'          Therefore, the patterns detected might be very weak in raw expression value terms
#'            "restandardized version of the maxmean statistic, in which we center
#'            and scale the maxmean statistic by its mean and standard deviation under
#'            the row randomizations" -Bradley Efron and Rob Tibshirani. Tech report. August 2006 
#'
#' @param rebuild  Rebuild the annotation files (DEFAULT: TRUE)
#' @param recreateFIRMA  Set to TRUE to rebuild the FIRMA files (DEFAULT: FALSE)
#' @param nperms  The number of permutations for the GSA algorithm (DEFAULT: 300)
#' @param cutoff.GSA  The GSA cutoff (DEFAULT: 1)
#' @param dataSetName  The name of the dataset, alt 'ALZ.042214', 'ALZ.071713' (DEFAULT: 'ALZ.040714')
#' @param runExonDif  SEt to true to calculate non-FIRMA differential exon expression (DEFAULT: FALSE)
#' @param useDABG  Use the DABG filter (DEFAULT: TRUE)
#' @param removeNAs  Remove genes with NA annotation (DEFAULT: TRUE)
#' @export
#' @usage outFiles <- buildALZDenrich(rebuild=TRUE, recreateFIRMA=FALSE, nperms=300, cutoff.GSA=1, dataSetName='ALZ.040714', loadFIRMA=TRUE, runExonDif=FALSE, useDABG=TRUE, removeNAs=TRUE)
buildALZDenrich <- function(rebuild=TRUE, recreateFIRMA=FALSE, nperms=300, cutoff.GSA=1, dataSetName='ALZ.040714', loadFIRMA=TRUE, runExonDif=FALSE, useDABG=TRUE, removeNAs=TRUE) {
	fsScores <- getFIRMAandGENE(dataSetName=dataSetName, subSet=NULL, recreate=recreateFIRMA, chipType="HuEx-1_0-st-v2", tags=NULL, loadFIRMA=loadFIRMA)
	#exonExpr.vals <- fsScores$exonExpr.vals 
	# Maybe we'll need exonExpr later.
	
	geneExpr.vals <- fsScores$geneExpr.vals #get(load('./intermediate/HIVRAD.HuEx-1_0-st-v2.coreR2,A20070914,EP..geneExpr.RData')[1])
	if (loadFIRMA==TRUE) { firmaExpr.vals <- fsScores$exonExpr.vals } else { firmaExpr.vals <- NULL }
	
	#Print the raw files
	outFile.rawGeneExpr <- paste('./intermediate/', paste(dataSetName, "rawGeneExpr.csv",sep="."), sep="")
	if (!file.exists(outFile.rawGeneExpr)) {
		write.csv(geneExpr.vals, file=outFile.rawGeneExpr)
	}
	
	exonExpr.vals <- fsScores$exonExpr.vals #get(load('./intermediate/HIVRAD.HuEx-1_0-st-v2.coreR2,A20070914,EP..geneExpr.RData')[1])
	outFile.rawFirma <- paste('./intermediate/', paste(dataSetName, "rawFIRMA.csv",sep="."), sep="")
	if (!file.exists(outFile.rawFirma) && !is.null(exonExpr.vals)) {
		write.csv(exonExpr.vals, file=outFile.rawFirma)
	}
	
	rm(fsScores)
	gc(reset=TRUE)

	#Use DAGB filter if one is available
	outString <- dataSetName
	if (useDABG == TRUE) {
		outString <- paste(outString, 'useDABG', sep='.')
		psIDs <- loadDABGfilter(dataSetName=dataSetName, cutoff = 0.01)
	} else {
		psIDs <- geneExpr.vals$unitName
	}
	geneExpr.vals <- geneExpr.vals[geneExpr.vals$unitName %in% psIDs,]
	exonExpr.vals <- exonExpr.vals[exonExpr.vals$unitName %in% psIDs,]

	#Remove NAs 
	if (removeNAs == TRUE) {
		outString <- paste(outString, 'removeNA', sep='.')
		geneExpr.vals <- geneExpr.vals[!is.na(geneExpr.vals$symbol),]
		exonExpr.vals <- exonExpr.vals[!is.na(exonExpr.vals$symbol),]
	} 
	
	#Calculate experiment annotation enrichment
	geneData <- exonToGene(geneExpr.vals)	

	lut <- csfLookup()
	colNames.lut <- sub('\\.', '-', sapply(strsplit(colnames(geneData), '_'), function(x) {x[2]}))
	colnames(geneData) <- lut[colNames.lut]
	
	for (curExp in c('^mci', '^premci')) {  
		if (!any(grepl(curExp, tolower(colnames(geneData))))) {next;}
	
		x <- geneData[, grepl(curExp, tolower(colnames(geneData))) |  grepl('normal', tolower(colnames(geneData)))]
		#colnames(x) <- csfLookup()[sapply(strsplit(colnames(x), '_'), function(x) {x[2]})]
		y <- as.numeric(grepl('normal', tolower(colnames(x)))) + 1
		
		if (all(y==1)) {return()}  #If there are only controls, skip this loop

		outFile.geneExpr <- paste('./intermediate/', paste(outString, curExp, "geneExpr.csv",sep="."), sep="")
		if (!file.exists(outFile.geneExpr) | rebuild==TRUE) {
			write.csv(round(x), file=outFile.geneExpr)
		}

		outFile.geneEnrich <- paste('./intermediate/', paste(outString, curExp, "geneExpr.Enrich.csv",sep="."), sep="")
		if (!file.exists(outFile.geneEnrich) | rebuild==TRUE) {
			annotEnrich <- AnnotateGeneExpr(x=x, y=y, cutoff=cutoff.GSA, nperms=nperms)
			write.csv(annotEnrich, file=outFile.geneEnrich)
		}

		#Calculate the gene level expression changes
		normalNames <- colnames(x)[grepl('^normal', tolower(colnames(x)))]
		expNames <- colnames(x)[!grepl('^normal', tolower(colnames(x)))]
		testCols <- which(colnames(x) %in% expNames)
		controlCols <- which(colnames(x) %in% normalNames)
		outFile.gene <- paste('./intermediate/', paste(outString, curExp, "geneExpr.pVal.csv",sep="."), sep="")
		if (file.exists(outFile.gene)) {
			geneExpr.sum <- read.csv(outFile.gene)
		} else {
			#pVals <- compareProbeCols(geneExpr.vals, testCols, controlCols, multProbeMethod=1, wilcox=TRUE)
			pVals <- compareProbeCols(as.data.frame(x), testCols, controlCols, multProbeMethod=4, wilcox=FALSE)
			save(pVals, file=outFile.gene)
			geneExpr.sum <- cbind(geneExpr.vals[match(rownames(pVals), geneExpr.vals$symbol), c('symbol', 'symbol_description')], 
				data.frame(t.qVals=p.adjust(pVals$t.log.pVals, 'fdr'), log2ExpOnCont=log2(pVals$meanRatios), absLog2=abs(log2(pVals$meanRatios)), pVals))
			geneExpr.sum <- geneExpr.sum[order(geneExpr.sum$t.pVal),]
			write.csv(geneExpr.sum, file=outFile.gene)
		}

		for (exonType in c('exonExpr', 'firmaExpr')) {
			lut <- csfLookup()

			if ((exonType == 'exonExpr')) { 
				if(runExonDif == FALSE) { next; }
				origData <- cbind(data.frame(SystematicName=rownames(geneExpr.vals)), geneExpr.vals)
			} #if ((exonType == 'exonExpr')) { 
			if (exonType == 'firmaExpr') { 
				if(loadFIRMA == FALSE) { next; }
				origData <- cbind(data.frame(SystematicName=rownames(firmaExpr.vals)), firmaExpr.vals)
			}

			#Calculate the exon level expression changes
			colNames.lut <- sapply(strsplit(colnames(origData), '_'), function(x) {x[2]})
			normalNames <- names(lut)[grepl('^normal', tolower(lut))]
			expNames <- names(lut)[!grepl('^normal', tolower(lut))]
			testCols <- which(colNames.lut %in% expNames)
			controlCols <- which(colNames.lut %in% normalNames)
			outFile.exon <- paste('./intermediate/', paste(outString, curExp, exonType, "pVal.csv.gz",sep="."), sep="")
			outFile.exon.lte05 <- paste('./intermediate/', paste(outString, curExp, exonType, "pVal.lte05.csv",sep="."), sep="")
			if (file.exists(outFile.exon)) {
				geneExpr.sum <- read.csv(outFile.exon)
			} else {
				#This is took like two days to complete.  Have looked like the delay stitching results together.
				pVals <- compareProbeCols(origData, testCols, controlCols, multProbeMethod=4, wilcox=FALSE)
				save(pVals, file=sub('csv.gz$', 'RData', outFile.exon))
				#This will require a better LUT -- 4/28/14 - when this finishes, load file and fix symbol lookup
				#I should include unit.group.cell
				exonExpr.sum <- cbind(origData[match(rownames(pVals), origData$SystematicName), c('symbol', 'symbol_description', 'unitName', 'unit', 'group', 'cell')], 
					data.frame(t.qVals=p.adjust(pVals$t.log.pVals, 'fdr'), log2hemiOncar=log2(pVals$meanRatios)), absLog2=abs(log2(pVals$meanRatios)), pVals)
				exonExpr.sum <- exonExpr.sum[order(exonExpr.sum$absLog2, decreasing = TRUE),]
				write.csv(exonExpr.sum, file=gzfile(outFile.exon))

				exonExpr.sum.lte05 <- exonExpr.sum[exonExpr.sum$t.log.pVals <= 0.05, ]
				write.csv(exonExpr.sum.lte05, file=outFile.exon.lte05)
			} #if (file.exists(outFile.exon)) {
		} #for (exonType in c('exonExpr', 'firmaExpr')) {
	} # for (curExp in c('^MCI, ^preMCI')) {	
	
	#To do: return a list with all of the things calculated
	return(NULL)
}

#' Annotate the CSF data from a file
#'
#' @param fileName  Rebuild the annotation files (DEFAULT: TRUE)
#' @export
#' @usage expLUT <- csfLookup()
csfLookup <- function() {
	#pre MCI
	#ACBQ293
	##2AC/BQ071
	#BQ/2AC0048
	##2AC126


	#normal
	#AC249
	#BQ/AC0486
	#BQ6165
	#BQ11002

	#expLUT <- c(paste('preMCI', 1:4, sep=""), paste('normal', 1:4, sep=""))
	#names(expLUT) <- c('ACBQ293', '2ACBQ071', 'BQ2AC0048', '2AC126', 'AC249', 'BQAC0486', 'BQ6165', 'BQ11002')
	
	curTable <- read.csv('./ALZ.SampleLUT.csv')
	expLUT <- sub('\\.', '', make.names(curTable$Type, unique=TRUE))
	names(expLUT) <- gsub('[/#-]', '',curTable$ID) #Remove special characters
	
	#Manually add in the old-style names from 20130715
	#20130715_0001-2_MCIAD_iCellneuron_24h.CEL
	#20130715_0013-2_MCIAD_iCellneuron_24h.CEL
	#20130715_0081-2_MCIAD_iCellneuron_24h.CEL
	#20130715_0312-2_MCIAD_iCellneuron_24h.CEL
	#20130715_0016-2_Normal_iCellneuron_24h.CEL
	#20130715_0019-2_Normal_iCellneuron_24h.CEL
	#20130715_0031-2_Normal_iCellneuron_24h.CEL
	#20130715_0060-2_Normal_iCellneuron_24h.CEL
	
	nMax <- max(as.numeric(sub('^Normal', '', expLUT[grepl('^Normal', expLUT)])), na.rm=TRUE)
	mciMax <- max(as.numeric(sub('^MCI', '', expLUT[grepl('^MCI', expLUT)])), na.rm=TRUE)
	newLUT <- c(paste('MCI', (mciMax+1):(mciMax+4), sep=""), paste('Normal', (nMax+1):(nMax+4), sep=""))
	names(newLUT) <- c('0001-2','0013-2','0081-2','0312-2','0016-2','0019-2','0031-2','0060-2')
	
	expLUT <- c(expLUT, newLUT)
	
	return(expLUT)
}

#' Load the DAGB scores calculated by Affy power tools.
#'   Figure out which experiments do not have all Normal or MCIAD predictions that pass cutoff
#'
#' @param dataSetName  The name of the dataset, alt 'ALZ.042214', 'ALZ.071713' (DEFAULT: 'ALZ.040714')
#' @param cutoff  The DABG cutoff (DEFAULT: 0.01)
#' @export 
#' @usage psIDs <- loadDABGfilter(dataSetName='ALZ.040714', cutoff = 0.01)
loadDABGfilter <- function(dataSetName='ALZ.040714', cutoff = 0.01) {
	inPath <- paste('./rawData/', dataSetName, '/', sep="")
	inFiles <- list.files(inPath, pattern='RMA-EXON-CORE-DABG.*\\.TXT$', include.dirs = TRUE, full.names = TRUE)
	outFile <- paste('./intermediate/DABG.', dataSetName, '.', as.character(cutoff), '.RData', sep="")

	if (!file.exists(outFile)) {
		dagbs.p.all <- NULL
		#Load all the p-Values into a single matrix
		for (inFile in inFiles) {
			dagbs <- read.csv(gzfile(inFile), sep='\t')
			dagbs.p <- dagbs[,grepl('\\.p\\.value', colnames(dagbs))]
			rownames(dagbs.p) <- dagbs[,'Probe.Set.ID']

			lut <- csfLookup()
			colNames.lut <- sub( '-rma.*$', '', sub('\\.', '-', sapply(strsplit(colnames(dagbs.p), '_'), function(x) {x[2]})))
			colnames(dagbs.p) <- lut[colNames.lut]
			if (is.null(dagbs.p.all)) {
				dagbs.p.all <- dagbs.p
			} else {
				dagbs.p.all <- cbind(dagbs.p.all, dagbs.p[rownames(dagbs.p.all),])
			}
		} #for (inFile in inFiles) {	
		dabgs.p <- dagbs.p.all

		#Multiple testing correction is inappropriate due to the 'all'
		#  If this is too difficult for anything to pass, then I'm going to have to implement
		#  a more clever algorithm.  Perhaps >= 90% pass the cutoff, or the fdr cutoff?
		NormalBool <- grepl('Normal', colnames(dagbs.p))
		MciadBool <- !grepl('Normal', colnames(dagbs.p))#grepl('MCIAD', colnames(dagbs.p))
		nInc <- apply(dagbs.p[,NormalBool], 1, function(x) {all(x<=cutoff)})
		mInc <- apply(dagbs.p[,MciadBool], 1, function(x) {all(x<=cutoff)})
		psIDs <- rownames(dagbs.p)[nInc | mInc]
		save(psIDs, file=outFile)
	} else {
		psIDs <- get(load(outFile)[1])
	} #if (file.exists(outFile)) {
	
	return(psIDs)
}

#' Create a ratios matrix for cMonkey python
#'
#' @export 
#' @usage makecMonkeyData.ALZ() 
makecMonkeyData.ALZ <- function() {
	#fsScores <- getFIRMAandGENE(dataSetName='ALZ.10serum.043015', chipType="HuEx-1_0-st-v2")
	fsScores <- getFIRMAandGENE(dataSetName='ALZ.10serum.58', chipType="HuEx-1_0-st-v2", loadFIRMA = FALSE)
	exonExpr <- fsScores$geneExpr.vals #get(load('./intermediate/HIVRAD.HuEx-1_0-st-v2.coreR2,A20070914,EP..geneExpr.RData')[1])
	

	#exonExpr <- get(load('./intermediate/ALZ.all081914.HuEx-1_0-st-v2.coreR2,A20070914,EP..geneExpr.RData')[1])
	expr <- exonToGene(exonExpr)
	#Look up names
	expLUT <- csfLookup()
	expNames <- sub('_[^_]*$', '', sub('^[^_]*_', '', colnames(expr)))
	colnames(expr) <- paste(expNames, expLUT[expNames], sep='.')
	
	# 1. Scherf JM, Hu XS, Tepp WH, Ichtchenko K, Johnson EA, Pellett S: Analysis of 
	#Gene Expression in Induced Pluripotent Stem Cell-Derived Human Neurons Exposed 
	#to Botulinum Neurotoxin A Subtype 1 and a Type A Atoxic Derivative. 
	#PLoS One 2014, 9.
	#
	#  Used iCell Neurons like we did. 
	#
	#  GSE58149, 
	#   GSM1402476	Non-exposed cells, sample 1, 2 weeks
	#   *GSM1402477	Non-exposed cells, sample 1, 2 days
	#   GSM1402478	Non-exposed cells, sample 2, 2 weeks
	#   *GSM1402479	Non-exposed cells, sample 2, 2 days
	#   GSM1402480	Non-exposed cells, sample 3, 2 weeks
	#   *GSM1402481	Non-exposed cells, sample 3, 2 days
	#
	#  Why is this data log scale? Jenn thinks the affy analysis package does that automatically
	
	if (FALSE) { #USE GEO
		data2 <- getDataFromSoft(GEO='GSE58149')
		ref2.names <- toupper(c('GSM1402477', 'GSM1402479', 'GSM1402481'))
		ref2 <- data2$curData[, ref2.names]
		curNames <- sub('\\..*$', '', rownames(ref2))
		ref2 <- apply(ref2, 2, function(x) {tapply(x, curNames, mean, na.rm=TRUE)})
	} else { #Use the reference
		data2 <- getFIRMAandGENE(dataSetName='ALZ.reference', chipType="HuEx-1_0-st-v2", loadFIRMA = FALSE)
		exonExpr.ref <- data2$geneExpr.vals #get(load('./intermediate/HIVRAD.HuEx-1_0-st-v2.coreR2,A20070914,EP..geneExpr.RData')[1])
			
		#exonExpr <- get(load('./intermediate/ALZ.all081914.HuEx-1_0-st-v2.coreR2,A20070914,EP..geneExpr.RData')[1])
		ref2 <- exonToGene(exonExpr.ref)
	} #if (FALSE) { #USE GEO

	commonGenes <- rownames(expr)[rownames(expr) %in% rownames(ref2)]
	refData <- ref2[commonGenes,]
	exprData <- expr[commonGenes,]
	
	#Quantile normalize the relevant parts of the data sets based on expr
	refData.norm <- t(normalizeTestData(t(refData), t(exprData)))		
	exprData <- exprData[rownames(refData.norm),]
	
	refRats <- NULL
	refRats.names <- NULL
	for (i in 1:(ncol(refData.norm)-1)) {
		col1 <- refData.norm[,i]
		for (j in 2:ncol(refData.norm)) {
			col2 <- refData.norm[,j]
			curRat <- col1 / col2
			curName <- paste('refRat', i, j,sep=".")
			if (is.null(refRats)) {
				refRats <- curRat
				refRats.names <- curName
			} else {
				refRats <- cbind(refRats, curRat)
				refRats.names <- c(refRats.names, curName)
			}
		}	
	} #for (i in 1:(ncol(refData.norm)-1)) {
	colnames(refRats) <- refRats.names
	refRats <- log2(refRats)
	means <- apply(refRats, 1, mean, na.rm=TRUE)
	sds <- apply(refRats, 1, sd, na.rm=TRUE)
	
	exprRats <- NULL
	for (i in 1:ncol(refData.norm)) {
		newdata <- exprData / refData.norm[,i]
		if (is.null(exprRats)) {
			exprRats <- newdata
		} else {
			exprRats <- cbind(exprRats, newdata)
		}
	}
	colnames(exprRats) <- make.names(colnames(exprRats), unique=TRUE)
	
	#Calculate a background distribution from the references, & generate median values.
	
	nlZs <- apply(exprRats, 2, function(x) { (x - means)/sds })
		
	#z-cut: Use an FDR (qvalue) cutoff of 0.05
	pVals <- pnorm(abs(nlZs), mean = 0, sd = 1, lower.tail = FALSE, log.p = FALSE)
	qVals <- t(apply(pVals, 1, function(x) {p.adjust(x, method='fdr')}))
	qVals[is.na(qVals)] <- 2 #An NA can never pass a p-value cut
	
	allRats <- cbind(data.frame(X=rownames(exprRats)), log2(exprRats))
	for (pCut in c(1, 0.05, 0.01, 0.001, 0.0001)) {
		for (ratCut in c(0, 0.5, 1, 1.5, 2)) {		
			incMatrix <- qVals <= pCut
			toDel <- rowSums(incMatrix)==0

			folMatrix <- abs(allRats[,2:ncol(allRats)]) >= ratCut
			toDel2 <- rowSums(folMatrix)==0
			
			toKeep <- !toDel & !toDel2
			cat('pCut =', pCut, '; ratCut =', ratCut, '; Num Genes =', sum(toKeep), '\n')

			outString <- paste('AD.log2ratios.pCut', pCut, 'ratCut', ratCut, 'numExp', ncol(exprRats), 'tsv', sep=".")

			write.table(allRats[toKeep,], sep='\t', file=outString, quote=FALSE, row.names=FALSE)
		} #for (ratCut in c(0, 0.5, 1, 1.5, 2)) {
	} #for (pCut in c(1, 0.05, 0.01, 0.001, 0.0001)) {
	
	#Results on 04-22-15: (Experimental Reference)
	# 0.0001 ->  7,717 genes
	# 0.001  -> 12,284 genes
	# 0.01   -> 17,840 genes
	# 0.05   -> 22,001 genes
	# 1.00   -> 24,422 genes
	
	
	#Results on 10-16-14: (GEO reference)
	# 0.001 ->  5,738 genes
	# 0.01  ->  9,508 genes
	# 0.05  -> 13,678 genes
	# 1.00  -> 17,899 genes
	
	#Results on 05-22-15: (Experimental Reference)
	#pCut = 1 ; ratCut = 0 ; Num Genes = 24421
	#pCut = 1 ; ratCut = 0.5 ; Num Genes = 21563
	#pCut = 1 ; ratCut = 1 ; Num Genes = 8283
	#pCut = 1 ; ratCut = 1.5 ; Num Genes = 2602
	#pCut = 1 ; ratCut = 2 ; Num Genes = 857
	#pCut = 0.05 ; ratCut = 0 ; Num Genes = 22000
	#pCut = 0.05 ; ratCut = 0.5 ; Num Genes = 19142
	#pCut = 0.05 ; ratCut = 1 ; Num Genes = 6235
	#pCut = 0.05 ; ratCut = 1.5 ; Num Genes = 1975
	#pCut = 0.05 ; ratCut = 2 ; Num Genes = 782
	#pCut = 0.01 ; ratCut = 0 ; Num Genes = 17839
	#pCut = 0.01 ; ratCut = 0.5 ; Num Genes = 14981
	#pCut = 0.01 ; ratCut = 1 ; Num Genes = 4068
	#pCut = 0.01 ; ratCut = 1.5 ; Num Genes = 1386
	#pCut = 0.01 ; ratCut = 2 ; Num Genes = 581
	#pCut = 0.001 ; ratCut = 0 ; Num Genes = 12283
	#pCut = 0.001 ; ratCut = 0.5 ; Num Genes = 9426
	#pCut = 0.001 ; ratCut = 1 ; Num Genes = 1928
	#pCut = 0.001 ; ratCut = 1.5 ; Num Genes = 730
	#pCut = 0.001 ; ratCut = 2 ; Num Genes = 388
	#pCut = 1e-04 ; ratCut = 0 ; Num Genes = 7716
	#pCut = 1e-04 ; ratCut = 0.5 ; Num Genes = 4990
	#pCut = 1e-04 ; ratCut = 1 ; Num Genes = 714
	#pCut = 1e-04 ; ratCut = 1.5 ; Num Genes = 332
	#pCut = 1e-04 ; ratCut = 2 ; Num Genes = 212
}

#' Cross-validate classifier based on cMonkey biclusters
#'    NOTE:  Another alternative is to treat clusters like GSA associations#'
#'    'inclusionMatrix.numeric.tsv', 'out.AD.001.numeric.3motifs.tsv'
#'    
#' @param fName  The output from 'extract_inclusion_matrix.py' (DEFAULT: 'inclusionMatrix.numeric.tsv')
#'
#' @export
#' @usage rv <- loadClusters(fName = 'inclusionMatrix.numeric.tsv')
loadClusters <- function(fName = 'inclusionMatrix.232exp.numeric.tsv') {
	rv <- list()
	
	iMatrix <- read.table(fName, row.names=1, header=TRUE)
	
	#Combine replicates (i.e. same experiment / n different controls)
	expNames <- colnames(iMatrix)
	uExpNames <- sub('\\.[0-9]$', '', expNames)
	iMatrixCombo <- apply(iMatrix, 1, function(x) {
		tapply(x, uExpNames, mean, na.rm=TRUE)
	})
	
	classes <- classes.bin <- sub('[0-9]+$', '', sub('^[^\\.]*\\.', '', rownames(iMatrixCombo)))
	classes.bin[classes.bin=='preMCI' | classes.bin=='MCI'] <- 'Disease'
	
	rv <- list(clusterScores = iMatrixCombo, classes = classes, classes.bin = classes.bin)
	return(rv)
}

#' Call xvalidateKaranAndCluster multiple times
#'
#' @param numTimes  The number of times to run the n-fold xvalidation (DEFAULT: 10)
#' @param nFolds  The n-fold cross validation (DEFAULT: 10)
#' @export 
#' @usage saveFiles <- multiple10Fold(numTimes=25, nFolds=10)
multiple10Fold <- function(numTimes=25, nFolds=10) {
	saveFiles <- NULL
	for (n in 1:numTimes) {
		predRes <- try(xvalidateKaranAndCluster(nFolds=nFolds))
		predRes <- try(xvalidateKaranAndCluster.break(curFile=predRes$saveFile))
		saveFiles <- c(saveFiles, predRes$saveFile)
	}
	
	return(saveFiles)
}

#' Perform a cross-validation of Karan's features, the cluster based features, and both together
#'
#' @param nFolds  The number of folds to xvalidate.  Will max out at the number of examples (DEFAULT:10)
#' @param consensusPreds  Set to TRUE to use the consensus predictions of top features (DEFAULT: FALSE)
#' @export 
#' @usage predRes <- xvalidateKaranAndCluster(nFolds=10, consensusPreds=FALSE)
xvalidateKaranAndCluster <- function(nFolds=10, consensusPreds=FALSE, clusterFile="", useMprobes=FALSE) {
	fsScores <- getFIRMAandGENE(dataSetName='ALZ.10serum.58', subSet=NULL, recreate=FALSE, chipType="HuEx-1_0-st-v2", tags=NULL, loadFIRMA=FALSE)
	exonData <- fsScores$geneExpr.vals #get(load('./intermediate/HIVRAD.HuEx-1_0-st-v2.coreR2,A20070914,EP..geneExpr.RData')[1])
	exprData <- exonToGene(exonData)
	exprNames <- sub('_.*$', '', sub('^[^_]*_', '', colnames(exprData)))
	
	expLUT <- csfLookup()
	classes <- expLUT[match(exprNames, names(expLUT))]
	classes <- sub('[0-9]+$', '', classes)
	colnames(exprData) <- paste(colnames(exprData), classes, sep='.')
	classes.bin <- rep('Normal', ncol(exprData))
	names(classes.bin) <- colnames(exprData)
	classes.bin[!grepl('Normal$', names(classes.bin))] <- 'Disease'

	genesets <- getGenesets(geneNames=rownames(exprData), loadExtra=TRUE, getClusterStack=FALSE)	
	
	clusters <- loadClusters(fName = 'inclusionMatrix.232exp.numeric.tsv')
	cFeats <- clusters$clusterScores
	rownames(cFeats) <- sub('[0-9]+$', '', rownames(cFeats))
	rownames(cFeats) <- sub('X2', '2', rownames(cFeats))

	resList <- list()	
	if (FALSE) {	
		preds <- makeKaransFeatures(trainSet=exprData, trainClass=classes.bin, testSet=trainSet, genesets=genesets, consensus=consensusPreds, makePred=FALSE)
		#preds <- makeKaransFeatures(trainSet, trainClass, testSet, genesets=genesets, consensus=FALSE, makePred=FALSE)

		#Work with this tuning code during xvalidations!
		#bestmtry <- tuneRF(data$train[-13],data$train$income, ntreeTry=100, 
		#   stepFactor=1.5,improve=0.01, trace=TRUE, plot=TRUE, dobest=FALSE)

		#Fake Xvalidate karan's features to get stats
		curRF <- randomForest(x=t(exprData[preds$bestGenes,,drop=FALSE]), y=factor(classes.bin), ntree=10000)
		oob.pVals <- curRF$votes[,'Disease', drop=FALSE]
		oob.pVals.class <- sub('[0-9]+$', '', sub('^.*\\.', '', names(oob.pVals)))
		oob.predStats <- getAUCandP(oob.pVals, oob.pVals.class != 'Normal')
		resList[['oob.predStats']] <- unlist(oob.predStats)

		#STOP 06-16-15 at 5:00 PM

		#Xvalidate the cluster features.  Note, this is real xvalidation, opposed to Karan's features which are based known classes
		curRF.clust <- randomForest(x=clusters$clusterScores, y=factor(clusters$classes.bin), ntree=1000)
		oob.pVals.clust <- curRF.clust$votes[,'Disease', drop=FALSE]
		oob.pVals.clust.class <- sub('[0-9]+$', '', sub('^.*\\.', '', names(oob.pVals.clust)))
		oob.predStats.clust <- getAUCandP(oob.pVals.clust, oob.pVals.clust.class != 'Normal')
		resList[['oob.predStats.clust']] <- unlist(oob.predStats.clust)

		#Concatenate the features.
		kFeats <- t(exprData[preds$bestGenes,, drop=FALSE])
		rownames(kFeats) <- sub( '_[0-9]+\\.', '.', sub('^[^_]*_', '', rownames(kFeats)))
		allFeats <- cbind(kFeats, cFeats[rownames(kFeats),, drop=FALSE])
		allFeats.class <- sub('^.*\\.', '', rownames(allFeats))
		allFeats.class[allFeats.class!='Normal'] <- 'Disease'
		curRF.all <- randomForest(x=allFeats, y=factor(allFeats.class), ntree=1000)
		oob.pVals.all <- curRF.all$votes[,'Disease']
		oob.pVals.all.class <- sub('[0-9]+$', '', sub('^.*\\.', '', names(oob.pVals.all)))
		oob.predStats.all <- getAUCandP(oob.pVals.all, oob.pVals.all.class != 'Normal')
		resList[['oob.predStats.all']] <- unlist(oob.predStats.all)

		#Concatenate the features. Force balanced #s of predictions
		kFeats.bal <- t(exprData[names(preds$geneScores)[1:ncol(cFeats)],, drop=FALSE])
		rownames(kFeats.bal) <- sub( '_[0-9]+\\.', '.', sub('^[^_]*_', '', rownames(kFeats)))
		allFeats <- cbind(kFeats.bal, cFeats[rownames(kFeats),, drop=FALSE])
		allFeats.class <- sub('^.*\\.', '', rownames(allFeats))
		allFeats.class[allFeats.class!='Normal'] <- 'Disease'
		curRF.bal <- randomForest(x=allFeats, y=factor(allFeats.class), ntree=1000)
		oob.pVals.bal <- curRF.bal$votes[,'Disease', drop=FALSE]
		oob.pVals.bal.class <- sub('[0-9]+$', '', sub('^.*\\.', '', names(oob.pVals.bal)))
		oob.predStats.bal <- getAUCandP(oob.pVals.bal, oob.pVals.bal.class != 'Normal')
		resList[['oob.predStats.bal']] <- unlist(oob.predStats.bal)
	}
	
	#xvalidate
	
	#Fake svalidate to see if there's an improvement.
	
	#xvalidate for real with mProbes if there's an improvement
	if (nFolds > ncol(exprData)) { nFolds <- ncol(exprData) }
	folds <- cut(1:ncol(exprData), nFolds)
	uFolds <- unique(folds)
	folds <- folds[sample.int(length(folds))]
	
	#NOTE FOR AUC: KSVM builds a probabilistic model using 3-fold xvalidation to fit parameters on a Laplacian model
	#   This seems to be much more powerful that using the votes from the randomForest
	allPreds <- foreach (testFold = uFolds) %dopar% {
		#To Do: Add Karan alone and cluster alone predictions as well
	
		cat(which(uFolds == testFold), '\n')
		testBool <- folds == testFold
		trainBool <- ! testBool
		trainSet <- exprData[, trainBool, drop=FALSE]
		trainClass <- classes.bin[trainBool]
		testSet <- exprData[, testBool, drop=FALSE]
	
		preds <- makeKaransFeatures(trainSet, trainClass, testSet, genesets=genesets, consensus=consensusPreds, makePred=FALSE)
		#allPreds[[testFold]] <- preds
		kFeats <- t(trainSet[preds$bestGenes,, drop=FALSE])
		rownames(kFeats) <- sub( '_[0-9]+\\.', '.', sub('^[^_]*_', '', rownames(kFeats)))

		kFeats.test <- t(testSet[preds$bestGenes,, drop=FALSE])
		rownames(kFeats.test) <- sub( '_[0-9]+\\.', '.', sub('^[^_]*_', '', rownames(kFeats.test)))
	
		all.train <- cbind(kFeats, cFeats[rownames(kFeats),, drop=FALSE])
		all.test <- cbind(kFeats.test, cFeats[rownames(kFeats.test),, drop=FALSE])

		train.class <- sub('^.*\\.', '', rownames(all.train))
		train.class[train.class != 'Normal'] <- 'Disease'
	
		x <- all.train
		mtry <- tuneRF(x=x, y=factor(train.class), ntreeTry=500,  improve=0.01, plot=FALSE)
		mtry.num <- mtry[which.min(mtry[,'OOBError']),'mtry']
		curRF <- randomForest(x=x, y=factor(train.class), ntree=1000, mtry=mtry.num)
		
		x <- all.test
		curPreds <- predict(curRF, newdata=x, type='prob')
		
		#Supersample
		if (FALSE) {
			x <- all.train
			x.more <- x[train.class=='Normal',, drop=FALSE]
			x <- rbind(x, x.more)
			y <- c(train.class, rep('Normal', nrow(x.more)))
		
			mtry <- tuneRF(x=x, y=factor(y), ntreeTry=500,  improve=0.01, plot=FALSE)
			mtry.num <- mtry[which.min(mtry[,'OOBError']),'mtry']
			curRF.ss <- randomForest(x=x, y=factor(y), ntree=1000, mtry=mtry.num)
		
			x <- all.test
			curPreds.ss <- predict(curRF.ss, newdata=x, type='prob')
		} else {
			curPreds.ss <- NULL
		}
		
		#Regression model for p-Values
		x <- all.train
		y <- as.numeric(train.class!='Normal')
		
		mtry <- suppressWarnings(tuneRF(x=x, y=y, ntreeTry=500,  improve=0.01, plot=FALSE))
		mtry.num <- mtry[which.min(mtry[,'OOBError']),'mtry']
		curRF.reg <- suppressWarnings(randomForest(x=x, y=y, ntree=1000, mtry=mtry.num))
				
		x <- all.test
		curPreds.reg <- predict(curRF.reg, newdata=x)
		curPreds.reg[curPreds.reg < 0] <- 0
		curPreds.reg[curPreds.reg > 1] <- 1
		curPreds.reg <- data.frame(Disease=curPreds.reg, Normal=1-curPreds.reg)
		rownames(curPreds.reg) <- rownames(x)
		
		#Reweight
		if(FALSE) {	#DOES NOT WORK AS EXPECTED
			x <- all.train
			y <- factor(train.class)
			classwt <- table(y)
			classwt <- 1 - classwt/sum(classwt)
			classwt <- classwt^5

			mtry <- tuneRF(x=x, y=factor(y), ntreeTry=500,  improve=0.01, plot=FALSE, classwt=classwt)
			mtry.num <- mtry[which.min(mtry[,'OOBError']),'mtry']
			curRF.wt <- randomForest(x=x, y=factor(y), ntree=2000, mtry=mtry.num, classwt=classwt)

			x <- all.test
			curPreds.wt <- predict(curRF.wt, newdata=x, type='prob')
		} else {
			curPreds.wt <- NULL
		}

		#ksvm
		if (FALSE) {
			x <- all.train
			y <- factor(train.class)
			classwt <- table(y)
			classwt <- 1 - classwt/sum(classwt)
			curSVM <- ksvm(x=x, y=y, class.weights = classwt, kernel='polydot', prob.model = TRUE)

			x <- all.test
			curPreds.ksvm <- predict(curSVM, newdata=x, type='probabilities')
			rownames(curPreds.ksvm) <- rownames(x)
		} else {
			curPreds.ksvm <- NULL
		}

		
		#mProbes
		if (useMprobes == TRUE) {
			dataMat <- cbind(all.train, data.frame(class=train.class))
			curProbes <- mProbes(dataMat, classIdx = ncol(dataMat), 
					type = "rFerns", dynamicCutoff = Inf, ntree = 1000, mcWorkers=1)
			bestFeat <- names(curProbes[curProbes<1])

			x.mp <- all.train[, bestFeat]
			mtry <- tuneRF(x=x.mp, y=factor(train.class), ntreeTry=500,  improve=0.01, plot=FALSE)
			mtry.num <- mtry[which.min(mtry[,'OOBError']),'mtry']
			curRF.mp <- randomForest(x=x.mp, y=factor(train.class), ntree=1000, mtry=mtry.num)

			x.mp <- all.test[, bestFeat]
			curPreds.mp <- predict(curRF.mp, newdata=x.mp, type='prob')
		} else {
			curPreds.mp <- NULL
		}
		
		rv <- list(all.train=all.train, all.test=all.test, train.class=train.class, 
			curPreds.rf=curPreds, curPreds.ss=curPreds.ss, curPreds.mp=curPreds.mp, 
			curPreds.wt=curPreds.wt, curPreds.ksvm=curPreds.ksvm, curPreds.reg=curPreds.reg)
		return(rv)
	}
	
	idx <- 0
	#It might be good to include smoothing and the minimum number of features
	saveFile <- paste('preds.xvalidateKaransAndClust.nFolds', nFolds, Sys.Date(), idx, 'RData', sep='.')
	while(file.exists(saveFile)) {
		idx <- idx+1
		saveFile <- sub('[0-9]+\\.RData', paste(idx, 'RData', sep="."), saveFile)
	}
	save(allPreds, file=saveFile)
	
	cat('Optimizing Parameters for ML Methods with caret.\n')
	for (i in 1:length(allPreds)) {
		cat(i, '\n')
		all.train <- allPreds[[i]]$all.train
		all.test <- allPreds[[i]]$all.test
		train.class <- allPreds[[i]]$train.class
	
		#carat package for better? probabilistic RF model?
		#http://stats.stackexchange.com/questions/12425/creating-a-certainty-score-from-the-votes-in-random-forests
		#caret seems to make it possible to more agressively tune
		#caret methods: http://topepo.github.io/caret/modelList.html
		#Note: If included above, this will explode due to multiple threads called withing multiple threads.
		x <- all.train
		y <- factor(train.class)
		#Note: summaryFunction=twoClassSummary means 'AUC' for some reason
		trControl <- trainControl(method='cv',number=10, classProbs = TRUE, summaryFunction=twoClassSummary)
		rf.caret <- train(x=x, y=y, method='rf', metric='ROC', TuneLength=3, trControl=trControl)
		svm.caret <- train(x=x, y=y, method='svmPoly', metric='ROC', TuneLength=3, trControl=trControl)
		svm.rad.caret <- train(x=x, y=y, method='svmRadial', metric='ROC', TuneLength=3, trControl=trControl)

		#Predict class probabilities (i.e. 'certainty' scores)
		x <- all.test
		curPreds.rf.caret <- predict(rf.caret, x, "prob")
		curPreds.svm.caret <- predict(svm.caret, x, "prob")
		curPreds.svm.rad.caret <- predict(svm.rad.caret, x, "prob")
		rownames(curPreds.svm.caret) <- rownames(curPreds.svm.rad.caret) <- rownames(x)

		allPreds[[i]]$curPreds.rf.caret=curPreds.rf.caret
		allPreds[[i]]$curPreds.svm.caret=curPreds.svm.caret
		allPreds[[i]]$curPreds.svm.rad.caret=curPreds.svm.rad.caret
	}
	
	allPredList <- list()
	allPredList[['allPreds.rf']] <- do.call(rbind, lapply(allPreds, function(x) {x$curPreds.rf}))
	allPredList[['allPreds.ss']] <- do.call(rbind, lapply(allPreds, function(x) {x$curPreds.ss}))
	allPredList[['allPreds.wt']] <- do.call(rbind, lapply(allPreds, function(x) {x$curPreds.wt}))
	allPredList[['allPreds.ksvm']] <- do.call(rbind, lapply(allPreds, function(x) {x$curPreds.ksvm}))
	allPredList[['allPreds.reg']] <- do.call(rbind, lapply(allPreds, function(x) {x$curPreds.reg}))
	allPredList[['curPreds.rf.caret']] <- do.call(rbind, lapply(allPreds, function(x) {x$curPreds.rf.caret}))
	allPredList[['curPreds.svm.caret']] <- do.call(rbind, lapply(allPreds, function(x) {x$curPreds.svm.caret}))
	allPredList[['curPreds.svm.rad.caret']] <- do.call(rbind, lapply(allPreds, function(x) {x$curPreds.svm.rad.caret}))
	
	predResults <- list()
	for(curname in names(allPredList)) {
		allPreds.reg <- allPredList[[curname]]
		allClasses.reg <- sub('[0-9]+$', '', sub('^.*\\.', '', rownames(allPreds.reg)))
		allClasses.disease <- rep(FALSE, length(allClasses.reg))
		allClasses.disease[allClasses.reg != 'Normal'] <- TRUE
		names(allClasses.disease) <- rownames(allPreds.reg)

		allPreds.dPred <- allPreds.reg[,'Disease']
		names(allPreds.dPred) <- rownames(allPreds.reg)

		predResults[[curname]] <- unlist(getAUCandP(pClass=allPreds.dPred, isClass=allClasses.disease))
	}
	save(allPreds, predResults, saveFile, file=saveFile)

	return(list(allPreds=allPreds, predResults=predResults, saveFile=saveFile))
}

#' Calculate the ROC stats including AUC and associated p-Values
#'
#' @param pClass  The probability of being in 'class'
#' @param isClass  A boolean vect for if the example is actually 'class'
#' @param titleStr  Set to a string to plot the AUC (DEFAULT: NULL)
#' @param getMaxSpecificity  Set to true to calculate the maximum significant specificity (DEFAULT: FALSE)
#' @export 
#' @usage preds <- getAUCandP(pClass, isClass, titleStr=NULL, getMaxSpecificity=FALSE)
getAUCandP <- function(pClass, isClass, titleStr=NULL, getMaxSpecificity=FALSE) {

	tps <- sum((pClass >= 0.5) & isClass)
	tns <- sum((pClass < 0.5) & !isClass)
	fps <- sum((pClass >= 0.5) & !isClass)
	fns <- sum((pClass < 0.5) & isClass)
	pVals <- getF1mcc(tps, fps, tns, fns)
	
	pVals$mcc.p <- try(cor.test(as.numeric(pClass >= 0.5), as.numeric(isClass))$p.value, silent=TRUE)
	if (class(pVals$mcc.p) == "try-error") { pVals$mcc.p <- as.numeric(NA) }
	
	if(getMaxSpecificity == TRUE) {
		maxspec <- maxSpecificity(pClass, isClass)
		pVals$maxSpecificity <- maxspec$maxSpecificity
		pVals$maxSpec.p <- maxspec$bestMCC.p
	} else {
		pVals$maxSpecificity <- as.numeric(NA)
	}
	
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
		
		if(getMaxSpecificity == TRUE) {
			xcoord <- 1-maxspec$maxSpecificity
			lines(c(xcoord,xcoord), c(0,1), lty=5, lwd=2)
			#curStats <- maxspec$resDF[maxspec$bestIdx,]
			curLabel <- paste('Max Specificity = ', round(maxspec$maxSpecificity,3),
				', Cutoff = ', round(maxspec$bestCutoff,3), 
				'\nMCC = ', round(maxspec$bestMCC,3), 
				', P = ', round(maxspec$bestMCC.p,3), sep='')
			if (xcoord > 0.6) { adj=1; off=-0.025 } else { adj=0; off=0.025 }
			text(x=xcoord+off, y=0.15, labels=curLabel, adj=adj)	
		}
		
	} else {
		pVals$auc <- as.numeric(NA)
		pVals$auc.pVal <- as.numeric(NA)
	}
	
	return(pVals)
}


#' For a ROC curve, find the maximum specificity that's
#'   1) The maximum specificity that passes an MCC p-value cutoff of 0.05
#'   or, failing that
#'   2) Maximizes MCC * specificity
#'   NOTE: Specificity = TP / TP + FP where Positive means has disease
#'
#' @param pClass  The probability of being in 'class'
#' @param isClass  A boolean vect for if the example is actually 'class'
#' @export 
#' @usage spec <- maxSpecificity(pClass, isClass)
maxSpecificity <- function(pClass, isClass) {
	cutoffs <- unique(sort(pClass))
	oc <- c(cutoffs[2:length(cutoffs)], 1)
	cutoffs <- (cutoffs + oc)/2
	cutoffs <- cutoffs[1:(length(cutoffs)-1)]
	
	resDF <- NULL
	for (cutoff in cutoffs) {
		cBool <- pClass > cutoff
		tps <- sum(cBool & isClass)
		tns <- sum(!cBool & !isClass)
		fps <- sum(cBool & !isClass)
		fns <- sum(!cBool & isClass)
		
		pVals <- getF1mcc(tps, fps, tns, fns)			
		pVals$mcc.p <- cor.test(as.numeric(pClass > cutoff), as.numeric(isClass))$p.value
	
		statDF <- data.frame(cutoff=cutoff, specificity=pVals$specificity, mcc=pVals$mcc, mcc.p=pVals$mcc.p)
		resDF <- rbind(statDF, resDF)
	}
	
	isSig <- resDF$mcc.p <= 0.05
	if (any(isSig)) {
		bestIdx <- which(isSig)[1]
	} else {
		resDF$mccXspec <- resDF$specificity * resDF$mcc
		bestIdx <- which.max(resDF$mccXspec)
	}
	specificity <- resDF$specificity[bestIdx]
	idxs <- resDF$specificity == specificity
	bestMCC <- max(resDF[idxs, 'mcc'])
	bestMCC.p <- min(resDF[idxs, 'mcc.p'])
	bestIdx <- which((resDF[,'mcc'] == bestMCC) & (resDF[,'specificity'] == specificity))[1]
	bestCutoff <- resDF[bestIdx, 'cutoff']
	
	return(list(resDF=resDF, bestIdx=bestIdx, maxSpecificity=specificity, 
		bestMCC=bestMCC, bestMCC.p=bestMCC.p, bestCutoff=bestCutoff))
	
}

#' Perform a cross-validation of Karan's features
#'   Note: AnnotateGeneExpr may screw up when it's multithread.  
#' TODO: Stop minizing OOB error and instead maximize AUC or MCC
#'
#' @param trainSet  The training data
#' @param trainClass  The training class
#' @param testSet  The testing data
#' @param genesets  The gene sets (DEFAULT: NULL)
#' @param makePred  Set to TRUE to include predictions (DEFAULT: TRUE)
#' @param consensusPreds  Set to TRUE to use the consensu predictions of top features (DEFAULT: TRUE)
#' @param smoothErr  Smooth out the oob error beforing picking the number of features (DEFAULT: TRUE)
#' @export 
#' @usage karan.preds <- makeKaransFeatures(trainSet, trainClass, testSet, genesets=NULL, makePred=TRUE, consensusPreds=TRUE, smoothErr=TRUE)
makeKaransFeatures <- function(trainSet, trainClass, testSet, genesets=NULL, makePred=TRUE, consensusPreds=TRUE, smoothErr=TRUE) {
	#uclasses <- unique(trainClass)
	uc1.bool <- grepl("Disease", trainClass)
	uc2.bool <- !uc1.bool
	p.vals <- apply(trainSet, 1, function(x) {
		res <- try(t.test(x[uc1.bool], x[uc2.bool])$p.value, silent=TRUE)
		if (class(res) == "try-error") { 
			browser()
			#res <- 1
		}
		return(res)
	})
	
	foldRatios <- apply(trainSet, 1, function(x) {
		median(x[uc1.bool], na.rm=TRUE) / median(x[uc2.bool], na.rm=TRUE)
	})
	
	if (is.null(genesets)){
		genesets <- getGenesets(geneNames=rownames(trainSet), loadExtra=TRUE, getClusterStack=FALSE)
	}
	gsaScores <- AnnotateGeneExpr(x=trainSet, y = uc1.bool + 1, cutoff = 0, 
		nperms = 10, loadExtra = TRUE, genesets = genesets, getClusterStack = FALSE)
		
	gsa <- gsaScores[!is.na(gsaScores$gsaScore),]
	
	#There's got to be a faster way
	idx1 <- which(colnames(gsa) == 'gsaScore')
	idx2 <- which(colnames(gsa) == 'genes')
	gsaScores <- apply(gsa, 1, function(x) { 
		gsaGenes <- strsplit(as.character(x[idx2]), ';')[[1]]
		gsaScores.in <- rep(x[idx1], length(gsaGenes))
		names(gsaScores.in) <- gsaGenes
		return(gsaScores.in)
	})
	gs <- unlist(gsaScores)
	gns <- sub('\\.[0-9]+$', '', names(gs))
	gns <- sub('^.*\\. *', '', gns)
	gs2 <- tapply(gs, gns, max, na.rm=TRUE)
	gsaGenes <- as.numeric(gs2)
	names(gsaGenes) <- names(gs2)
	
	if (FALSE) {  #This is a very slow way to do the same thing
		gsaGenes <- list()
		for (i in 1:nrow(gsa)){
			curScore <- gsa$gsaScore[i]
			curGenes <- strsplit(as.character(gsa$genes[i]), '; ')[[1]]
			for (cg in curGenes) {
				if (is.null(gsaGenes[[cg]])) {
					gsaGenes[[cg]] <- curScore
				} else if (abs(gsaGenes[[cg]]) < abs(curScore)) {
					gsaGenes[[cg]] <- curScore
				}
			}
		} #for (i in 1:nrow(gsa)){
		gsaGenes <- unlist(gsaGenes)
	}
	
	xVal <- foldRatios
	yVal <- rep(0, length(xVal))
	names(yVal) <- names(xVal)
	yVal[names(gsaGenes)] <-  gsaGenes
	yVal <- abs(yVal[names(xVal)])
	zVal <- p.vals[names(xVal)]


	dataMat <- matrix(c(xVal, yVal, zVal), nrow = length(xVal), ncol = 3)
	
	rownames(dataMat) <- names(xVal)
	colnames(dataMat) <- c("Ratio", "ABS GSA Score", "pValues");

	if (FALSE) {
		if (includePlot) {
			plotPoints(dataMat);
		}
	}

	#Normalizes matrix
	normal <- sapply(colnames(dataMat), function(x) {
		val <- dataMat[,x]
		val <- (val - min(val))/(max(val) - min(val));
	});	

	geneScores <- normal[,1]^2 + normal[,2]^2 + ((1-normal[,3])^2)
	geneScores <- sort(geneScores, decreasing = T)
	
	
	#limitPos <- c(5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 
	#		55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 125, 150, 175,
	limitPos <- c(200, 225, 250, 275, 300, 325, 350, 375, 400, 425, 450, 475, 500,
			525, 550, 575, 600, 625, 650, 675, 700, 725, 750, 775, 800,
			825, 850, 875, 900, 925, 950, 975, 1000)
	
	
	#A better alternative may be to limit MSE for trinary predictions or maximize MCC
	#oob.MCC <- getF1mcc(tps=curRF$confusion[1,1], fps=curRF$confusion[2,1], 
	#			tns=curRF$confusion[2,2], fns=curRF$confusion[1,2])
	oobs <- NULL
	for (numFeat in limitPos) {
		topGenes <- names(geneScores)[1:numFeat]
		curRF <- randomForest(x=t(trainSet[topGenes,]), y=factor(trainClass))
		oobs <- c(oobs, mean(curRF$err.rate[,'OOB']))
	}
	
	#Smooth out bumps in the OOBS
	if (smoothErr == TRUE) {
		oobs <- smooth(oobs)
	}
	
	bestNumFeat <- limitPos[which.min(oobs)]
	bestGenes <- names(geneScores)[1:bestNumFeat]
	
	if(makePred == TRUE) {
		if (consensusPreds == FALSE) {
			#Work with this tuning code during xvalidations!
			#bestmtry <- tuneRF(data$train[-13],data$train$income, ntreeTry=100, 
		    	#   stepFactor=1.5,improve=0.01, trace=TRUE, plot=TRUE, dobest=FALSE)
	
		
			x <- t(trainSet[bestGenes,])
			curRF <- randomForest(x=x, y=factor(trainClass), ntree=5000)
	
			x <- t(testSet[bestGenes,, drop=FALSE])
			curPreds <- predict(curRF, newdata=x, type='prob') 
		} else {
			x <- cbind(t(trainSet[bestGenes,]),data.frame(Class=trainClass))
			x2 <- cbind(t(testSet[bestGenes,, drop=FALSE]),data.frame(Class='?'))
			#Past 200 or so, things slow down a lot, so count by 5s
			cutoff <- 300
			numsFeatures <- 1:length(bestGenes)
			if (length(bestGenes) > cutoff)  {
				numsFeatures <- 1:(cutoff-1)
				numsFeatures <- c(numsFeatures, seq(from = cutoff, to = length(bestGenes), by = 5))
			}
			conPreds.list <- makePreds(trainData=x, testData=x2, 
				numsFeatures = numsFeatures, onlySVM = FALSE,
				posClass = "Disease", negClass = "Normal", singleCore=TRUE)
			curPreds <- data.frame(Disease=conPreds.list$predTally$All$pHEMI,
						Normal=1-conPreds.list$predTally$All$pHEMI)
			rownames(curPreds) <- rownames(conPreds.list$predTally$All)
			
		}
	} else {
		curPreds <- NULL
	} #if(makePred == TRUE) {
	
	#To Do: return votes
	return(list(curPreds=curPreds, bestGenes=bestGenes, geneScores=geneScores))
}

#' Perform a cross-validation of Karan's features
#'   Note: AnnotateGeneExpr may screw up when it's multithread.  
#' TODO: Stop minizing OOB error and instead maximize AUC or MCC
#'
#' @param filePrefix  A regex used to select the 
#' @export 
#' @usage karan.preds <- tallyMultipleFold(filePrefix='preds.xvalidateKaransAndClust.nFolds.16.[\\.0-9\\-]*.karanClust.RData')
tallyMultipleFold <- function(filePrefix='preds.xvalidateKaransAndClust.nFolds.16.[\\.0-9\\-]*.karanClust.RData') {
	foldFiles <- list.files(path='.', pattern=filePrefix)
	
	#To Do: Make sure that each cross-validation is unique
	parts <- c('', 'cluster', 'karan')
	expNames <- c('curPreds.rf', 'curPreds.svm.rad.caret')
	predTallys <- list()
	for (curFile in foldFiles) {
		varNames <- load(curFile)
		allPreds <- get(varNames[grep('allPreds', varNames)[1]])
		predResults <- get(varNames[grep('predResults', varNames)[1]])
		
		for (part in parts) {
			for (curExp in expNames) {
				curDataName <- sub('^\\.', '', paste(part, curExp, sep="."))
				if (is.null(predTallys[[curDataName]])) {
					predTallys[[curDataName]] <- predResults[[curDataName]]
				} else {
					predTallys[[curDataName]] <- rbind(predTallys[[curDataName]], predResults[[curDataName]])
				}
			}
		} #for (part in parts) {
	} #for (curFile in foldFiles) {
	n <- nrow(predTallys[[1]])
	
	#For each classifier type
	#predTallyPvals <- list()
	curCols <- colnames(predTallys[[1]])
	prefixes <- c('rf$', 'svm.rad.caret$')
	allResList <- list()
	for (prefix in prefixes) {
		allResDF <- NULL
		for (curCol in curCols) {
			rfs <- names(predTallys)[grepl(prefix,names(predTallys))]
			#svms <- names(predTallys)[grepl('svm.rad.caret$',names(predTallys))]
			curRes <- sapply(predTallys[rfs], function(x) {x[,curCol]})
			ref <- colnames(curRes)[grepl('^curPreds', colnames(curRes))]
			refData <- curRes[,ref]
			test <- colnames(curRes)
			test <- test[test!=ref]
			medians <- apply(curRes, 2, median)
			means <- apply(curRes, 2, mean)
			pVals <- c(1, apply(curRes[,test], 2, function(x){wilcox.test(x, refData, paired=TRUE)$p.value}))
			resDF <- rbind(medians, means, pVals)
			rownames(resDF) <- paste(curCol, rownames(resDF), sep='.')
			allResDF <- rbind(allResDF, resDF)
		}
		allResList[[prefix]] <- allResDF
	} #for (prefix in prefixes) {
	
	outFile <- paste("tallyMultipleFold.wilcox", n, Sys.Date(), 'csv', sep='.')
	write.csv(allResList, file=outFile)
	
	return(allResList)
}

#' Look over the cross-validation runs and tally the genes that are selected.
#'   NOTE: All clusters are used in the classifier, but I could mProbes them to find the best
#'
#' @param filePrefix  A regex used to select the 
#' @export 
#' @usage features.sel <- tallyFeaturesSelected(filePrefix='preds.xvalidateKaransAndClust.nFolds.16.[\\.0-9\\-]*.karanClust.RData')
tallyFeaturesSelected <- function(filePrefix='preds.xvalidateKaransAndClust.nFolds.16.[\\.0-9\\-]*.karanClust.RData') {
	foldFiles <- list.files(path='.', pattern=filePrefix)
	
	#To Do: Make sure that each cross-validation is unique
	curFeatList <- NULL
	n <- length(foldFiles)
	allCounts <- 0
	for (curFile in foldFiles) {
		varNames <- load(curFile)
		allPreds <- get(varNames[grep('allPreds', varNames)[1]])
		predResults <- get(varNames[grep('predResults', varNames)[1]])
		allCounts <- allCounts + length(allPreds)
		
		curFeats <- unlist(lapply(allPreds, function(x) { colnames(x$all.test) }))
		curFeats <- curFeats[!grepl('^Cluster[0-9]+$', curFeats)]
		#curFeatList[[curFile]] <- sort(table(curFeats))
		curFeatList <- c(curFeatList, curFeats)
	} #for (curFile in foldFiles) {
	curFeatTable <- sort(table(curFeatList), decreasing = TRUE)
	curFeatDF <- data.frame(gene=names(curFeatTable), count=curFeatTable, 
		maxCount=allCounts, frac=curFeatTable/allCounts, 
		descr=addGeneDescr(names(curFeatTable)))
	
	#Add in some gene descriptions
	
	
	outFile <- paste("tallyMultipleFoldGenes", n, Sys.Date(), 'csv', sep='.')
	write.csv(curFeatDF, file=outFile)
	
	return(curFeatDF)
}

#' Look over the cross-validation runs and determine if MCI or preMCI predictions are better
#'
#' @param filePrefix  A regex used to select the 
#' @export 
#' @usage features.sel <- tallyMCIvPreMCIpreds(filePrefix='preds.xvalidateKaransAndClust.nFolds.16.[\\.0-9\\-]*.karanClust.RData')
tallyMCIvPreMCIpreds <- function(filePrefix='preds.xvalidateKaransAndClust.nFolds.16.[\\.0-9\\-]*.karanClust.RData') {
	foldFiles <- list.files(path='.', pattern=filePrefix)
	
	#To Do: Make sure that each cross-validation is unique
	expNames <- c('curPreds.rf', 'curPreds.svm.rad.caret')
	
	
	resList <- list()
	for (expName in expNames) {
		resList[[expName]] <- data.frame()
	}
	
	
	n <- length(foldFiles)
	for (curFile in foldFiles) {
		varNames <- load(curFile)
		allPreds <- get(varNames[grep('allPreds', varNames)[1]])
		predResults <- get(varNames[grep('predResults', varNames)[1]])
		
		for (i in 1:length(allPreds)) {
			x <- allPreds[[i]]
			for (expName in expNames) {
				resList[[expName]] <- rbind(resList[[expName]], as.data.frame(x[expName]))
			}
		}
	} #for (curFile in foldFiles) {
	
	accList <- list()
	for (expName in expNames) {
		dP <- colnames(resList[[expName]])[grepl('Disease$', colnames(resList[[expName]]))]
		resList[[expName]]$Pred.Disease <- resList[[expName]][,dP] >= 0.5
		resList[[expName]]$Type <- sub( '[0-9]+$', '', sub('^.*\\.', '', rownames(resList[[expName]])))
		resList[[expName]]$is.Disease <- resList[[expName]]$Type != 'Normal'
		mci.bool <- !resList[[expName]]$is.Disease | resList[[expName]]$Type == 'MCI'
		preMci.bool <- !resList[[expName]]$is.Disease | resList[[expName]]$Type == 'preMCI'
		accList[[paste(expName, 'MCI.acc', sep='.')]] <- getAUCandP(resList[[expName]][mci.bool,dP], resList[[expName]][mci.bool,'is.Disease'])
		accList[[paste(expName, 'MCI.acc', sep='.')]]$n <- length(resList[[expName]][mci.bool,dP])
		accList[[paste(expName, 'preMCI.acc', sep='.')]] <- getAUCandP(resList[[expName]][preMci.bool,dP], resList[[expName]][preMci.bool,'is.Disease'])
		accList[[paste(expName, 'preMCI.acc', sep='.')]]$n <- length(resList[[expName]][preMci.bool,dP])
	}
	
	accDF <- do.call(cbind, accList)
	
	colnames(accDF) <- sub('curPreds.', '', colnames(accDF))
	
	outFile <- paste("tallyMCIvPreMCIpreds", Sys.Date(), 'csv', sep='.')
	write.csv(accDF, file=outFile)
	
	return(accDF)
}


#' Get basic descriptions for a list of gene names
#'
#' @param geneNames  The list of gene names
#' @export 
#' @usage geneDescr <- addGeneDescr(geneNames)
addGeneDescr <- function(geneNames) {
	fsScores <- getFIRMAandGENE(dataSetName='ALZ.10serum.58', subSet=NULL, recreate=FALSE, chipType="HuEx-1_0-st-v2", tags=NULL, loadFIRMA=FALSE)
	gev <- fsScores$geneExpr.vals
	idxs <- match(geneNames, gev$symbol)
	geneDescr <- gev$symbol_description[idxs]
	names(geneDescr) <- geneNames
	return(geneDescr)
}	

#' Add the public data from GSE58149 and see if that improves the xvalidated specificity
#'   Note: First do the quick and dirty using the features from the 58
#'   Then do the full set.
#'
#' @param 
#' @export 
#' @usage geneDescr <- addPublicData()
addPublicData <- function() {

	#Load the 58 experiments.
	fsScores <- getFIRMAandGENE(dataSetName='ALZ.10serum.58', subSet=NULL, recreate=FALSE, chipType="HuEx-1_0-st-v2", tags=NULL, loadFIRMA=FALSE)
	exonData <- fsScores$geneExpr.vals #get(load('./intermediate/HIVRAD.HuEx-1_0-st-v2.coreR2,A20070914,EP..geneExpr.RData')[1])
	exprData <- exonToGene(exonData)
	exprNames <- sub('_.*$', '', sub('^[^_]*_', '', colnames(exprData)))
	
	expLUT <- csfLookup()
	classes <- expLUT[match(exprNames, names(expLUT))]
	classes <- sub('[0-9]+$', '', classes)
	colnames(exprData) <- paste(colnames(exprData), classes, sep='.')
	classes.bin <- rep('Normal', ncol(exprData))
	names(classes.bin) <- colnames(exprData)
	classes.bin[!grepl('Normal$', names(classes.bin))] <- 'Disease'

	genesets <- getGenesets(geneNames=rownames(exprData), loadExtra=TRUE, getClusterStack=FALSE)	

	resList <- list()			
	#FAKE RUN: NO XVALIDATION

	#Get Karan's features
	preds <- makeKaransFeatures(trainSet=exprData, trainClass=classes.bin, testSet=trainSet, genesets=genesets, consensus=FALSE, makePred=FALSE)
	
	#NOTE: It will be tricky to get bicluster scores for the new data
	
	#Get the cross-validated statistics & predictions.
	curRF <- randomForest(x=t(exprData[preds$bestGenes,,drop=FALSE]), y=factor(classes.bin), ntree=10000)
	oob.pVals <- curRF$votes[,'Disease', drop=TRUE]
	oob.pVals.class <- sub('[0-9]+$', '', sub('^.*\\.', '', names(oob.pVals)))
	oob.predStats <- getAUCandP(oob.pVals, oob.pVals.class != 'Normal')
	resList[['oob.predStats']] <- unlist(oob.predStats)

	#preds <- makeKaransFeatures(trainSet, trainClass, testSet, genesets=genesets, consensus=FALSE, makePred=FALSE)

	#Load the GSE58149 untreated samples
	data2 <- getDataFromSoft(GEO='GSE58149')
	ref2 <- data2$curData
	curNames <- sub('\\..*$', '', rownames(ref2))
	ref2 <- apply(ref2, 2, function(x) {tapply(x, curNames, mean, na.rm=TRUE)})
	
	probesets <- as.character(Table(GPLList(data2$curGEO)[[1]])$ID)
	GeneSymbol <- as.character(Table(GPLList(data2$curGEO)[[1]])$'Gene Symbol')
	names(GeneSymbol) <- probesets
	gns <- GeneSymbol[rownames(ref2)]
	GSE58149.all <- apply(ref2, 2, function(x) {tapply(x, gns, mean, na.rm=TRUE)})
	
	#Conormalize
	GSE58149.tst <- normalizeTestData(t(GSE58149.all), t(exprData))
	GSE58149.all <- t(GSE58149.tst)
	
	u.names <- toupper(c('GSM1402477', 'GSM1402479', 'GSM1402481', 'GSM1402476', 'GSM1402478', 'GSM1402480'))
	GSE58149.u <- GSE58149.all[,u.names]
	
	#Estimate using fake OOB
	commonGenes <- preds$bestGenes[preds$bestGenes %in% rownames(GSE58149.u)]
	uPlus58 <- cbind(exprData[commonGenes,,drop=FALSE], GSE58149.u[commonGenes,,drop=FALSE])
	newClass <- rep('Normal', ncol(GSE58149.u))
	names(newClass) <- colnames(GSE58149.u)
	classes.wu <- c(classes.bin, newClass)
	
	curRF.uw <- randomForest(x=t(uPlus58), y=factor(classes.wu), ntree=10000)
	oob.pVals <- curRF.uw$votes[,'Disease', drop=TRUE]
	oob.pVals.class <- sub('[0-9]+$', '', sub('^.*\\.', '', names(oob.pVals)))
	oob.predStats <- getAUCandP(oob.pVals, oob.pVals.class != 'Normal')
	resList[['oob.uPlus58']] <- unlist(oob.predStats)

	aPlus58 <- cbind(exprData[commonGenes,,drop=FALSE], GSE58149.all[commonGenes,,drop=FALSE])
	newClass <- rep('Normal', ncol(GSE58149.all))
	names(newClass) <- colnames(GSE58149.all)
	classes.wall <- c(classes.bin, newClass)
	
	curRF.wall <- randomForest(x=t(aPlus58), y=factor(classes.wall), ntree=10000)
	oob.pVals <- curRF.wall$votes[,'Disease', drop=TRUE]
	oob.pVals.class <- sub('[0-9]+$', '', sub('^.*\\.', '', names(oob.pVals)))
	oob.predStats <- getAUCandP(oob.pVals, oob.pVals.class != 'Normal')
	resList[['oob.aPlus58']] <- unlist(oob.predStats)

	#Cross-validate, 10-fold
	allData <- cbind(exprData[rownames(GSE58149.u),,drop=FALSE], GSE58149.u[,,drop=FALSE])
	
	nTimes <- 10
	
	rep10fold <- list()
	for (i in 1:nTimes) {
		#xvalidate for real with mProbes if there's an improvement
		nFolds <- 10
		folds <- cut(1:ncol(allData), nFolds)
		uFolds <- unique(folds)
		folds <- folds[sample.int(length(folds))]

		#NOTE FOR AUC: KSVM builds a probabilistic model using 3-fold xvalidation to fit parameters on a Laplacian model
		#   This seems to be much more powerful that using the votes from the randomForest
		allPreds <- foreach (testFold = uFolds) %dopar% {
			#To Do: Add Karan alone and cluster alone predictions as well
			cat(which(uFolds == testFold), '\n')
			testBool <- folds == testFold
			trainBool <- ! testBool
			trainSet <- allData[, trainBool, drop=FALSE]
			trainClass <- classes.wu[trainBool]
			testSet <- allData[, testBool, drop=FALSE]

			preds <- makeKaransFeatures(trainSet, trainClass, testSet, genesets=genesets, consensus=FALSE, makePred=TRUE)
			return(preds)
		}

		#Statistics
		pDisease <- NULL
		for (cp in allPreds) {
			pDisease <- c(pDisease, cp$curPreds[,'Disease'])
		}
		isClass <- grepl('MCI$', names(pDisease))
		preds <- getAUCandP(pDisease, isClass)

		rep10fold[[i]] <- unlist(preds)
	} #for (i in 1:nTimes) {
	
	save(rep10fold, file="rep10fold.073115.RData")
	
	return(rep10fold)
}


#' Load the clinical data
#'
#' @param 
#' @export 
#' @usage clinicalData <- AddClinicalData()
AddClinicalData <- function() {
	resList <- list()	
		
	allData <- read.csv('Clinical.AD080415.csv')
	
	#Turn it into a data matrix
	expNames <- allData$LabLabel
	expClass <- allData$Sample
	binClass <- as.character(expClass)
	binClass[binClass != 'Normal'] <- 'Disease'
	#Gender, Abeta, ptau
	APOEGen <- as.character(allData$APOEGen)
	APOEcount <- sapply(strsplit(APOEGen, ','), function(x) {sum(x == 4)})
	Gender <- allData$Gender
	dataMatrix <- data.frame(APOEcount=APOEcount, Gender=Gender)
	rownames(dataMatrix) <- expNames
	rownames(dataMatrix) <- gsub('[#-/]', '', rownames(dataMatrix))
	names(binClass) <- rownames(dataMatrix)
	
	#Filter to include only the correct 58
	fsScores <- getFIRMAandGENE(dataSetName='ALZ.10serum.58', subSet=NULL, recreate=FALSE, chipType="HuEx-1_0-st-v2", tags=NULL, loadFIRMA=FALSE)
	exonData <- fsScores$geneExpr.vals #get(load('./intermediate/HIVRAD.HuEx-1_0-st-v2.coreR2,A20070914,EP..geneExpr.RData')[1])
	exprData <- exonToGene(exonData)
	exprNames <- sub('_.*$', '', sub('^[^_]*_', '', colnames(exprData)))
	colnames(exprData) <- exprNames

	dataMatrix.apoe <- dataMatrix[exprNames,]
	binClass.apoe <- binClass[exprNames]


	#mTry is important because we only have two features
	curRF <- randomForest(x=dataMatrix.apoe, y=factor(binClass.apoe), 
		mtry =2, ntree=10000, importance=TRUE)
	oob.pVals <- curRF$votes[,'Disease', drop=TRUE]
	oob.pVals.class <- binClass[names(oob.pVals)]
	oob.predStats <- getAUCandP(oob.pVals, oob.pVals.class != 'Normal')
	resList[['APOE.Gender']] <- list(curRF=curRF, oob.predStats=unlist(oob.predStats))

	#NOTE: Classifier is horrible, try removing 'Gender'
	#> unlist(oob.predStats)
	#sensitivity specificity         fpr         fdr          f1         mcc
	#0.657894737 0.350000000 0.650000000 0.342105263 0.657894737 0.007894737
	#  agreement     p.value       mcc.p         auc    auc.pVal
	#0.551724138 0.179071651 0.953098206 0.640789474 0.134199310
	#> importance(curRF)
	#            Disease    Normal MeanDecreaseAccuracy MeanDecreaseGini
	#APOEcount  64.99904  64.93109             74.26110        4.9228227
	#Gender    -18.50353 -30.40305            -30.67835        0.9535036
	#
	dataMatrix.apoe.noGender <- dataMatrix.apoe[,'APOEcount',drop=FALSE]
	curRF <- randomForest(x=dataMatrix.apoe.noGender, y=factor(binClass.apoe), 
		mtry =1, ntree=10000)
	oob.pVals <- curRF$votes[,'Disease', drop=TRUE]
	oob.pVals.class <- binClass[names(oob.pVals)]
	pdf('APOE.Only.Classifier.pdf')
	oob.predStats <- getAUCandP(oob.pVals, oob.pVals.class != 'Normal', titleStr='APOE Only Classifier')
	dev.off()
	resList[['APOE']] <- list(curRF=curRF, oob.predStats=unlist(oob.predStats))

	#> unlist(oob.predStats)
	#sensitivity specificity         fpr         fdr          f1         mcc
	#0.657894737 0.750000000 0.250000000 0.166666667 0.735294118 0.387985285
	#  agreement     p.value       mcc.p         auc    auc.pVal
	#0.689655172 0.001116276 0.002618253 0.540789474 0.038355640

	#SVM?
	
	#Add in karan's and cluster's data
	expLUT <- csfLookup()
	classes <- expLUT[match(exprNames, names(expLUT))]
	classes <- sub('[0-9]+$', '', classes)
	colnames(exprData) <- paste(colnames(exprData), classes, sep='.')
	classes.bin <- rep('Normal', ncol(exprData))
	names(classes.bin) <- colnames(exprData)
	classes.bin[!grepl('Normal$', names(classes.bin))] <- 'Disease'

	genesets <- getGenesets(geneNames=rownames(exprData), loadExtra=TRUE, getClusterStack=FALSE)	

	clusters <- loadClusters(fName = 'inclusionMatrix.232exp.numeric.tsv')
	cFeats <- clusters$clusterScores
	rownames(cFeats) <- sub('[0-9]+$', '', rownames(cFeats))
	rownames(cFeats) <- sub('X2', '2', rownames(cFeats))

	#NOTE: Must be placed inside an xvalidation loop otherwise 
	#   will only give a theoretical maximum
	preds <- makeKaransFeatures(trainSet=exprData, trainClass=classes.bin, testSet=trainSet, genesets=genesets, consensus=consensusPreds, makePred=FALSE)
		
	#Concatenate the features.
	#  mProbes since cluster data is so large
	kFeats <- t(exprData[preds$bestGenes,, drop=FALSE])
	rownames(kFeats) <- sub( '_[0-9]+\\.', '.', sub('^[^_]*_', '', rownames(kFeats)))
	allFeats <- cbind(kFeats, cFeats[rownames(kFeats),, drop=FALSE])
	allFeats.mProbes.fdr.ferns <- mProbes(cbind(allFeats, 
		data.frame(class=allFeats.class)), ntree = 10000)
	#allFeats.mProbes.fdr <- allFeats.mProbes.fdr.rf <- mProbes(cbind(allFeats, data.frame(class=allFeats.class)), dynamicCutoff=2000)
	#keepFeats <- names(allFeats.mProbes.fdr)[allFeats.mProbes.fdr<1] #Note, rferns will not pass anything past .99
	#keepFeats <- names(allFeats.mProbes.fdr)[allFeats.mProbes.fdr<.99]
	ch <- hist(allFeats.mProbes.fdr.ferns) #Sturges' formula
	cutoff <- ch$breaks[length(ch$breaks)-1]
	keepFeats <- names(allFeats.mProbes.fdr.ferns)[allFeats.mProbes.fdr.ferns < cutoff]
	
	allFeats.class <- sub('^.*\\.', '', rownames(allFeats))
	allFeats.class[allFeats.class!='Normal'] <- 'Disease'
	curRF.all <- randomForest(x=allFeats[,keepFeats], y=factor(allFeats.class), ntree=1000)
	oob.pVals.all <- curRF.all$votes[,'Disease']
	oob.pVals.all.class <- sub('[0-9]+$', '', sub('^.*\\.', '', names(oob.pVals.all)))
	pdf('clustAndKaran.fake.pdf')
	oob.predStats.all <- getAUCandP(oob.pVals.all, oob.pVals.all.class != 'Normal',
		titleStr='Cluster and Karan, mProbes, NON-XVALIDATED OOB')
	dev.off()
	resList[['KaranAndCluster.fake']] <- list(curRF=curRF.all, oob.predStats=unlist(oob.predStats.all))

	#Concatenate with APOE
	apoe.sort <- dataMatrix.apoe.noGender[sub('\\..*$', '', rownames(allFeats)),]
	allFeats.apoe <- cbind(allFeats[,keepFeats], data.frame(apoe=apoe.sort))
	curRF.all <- randomForest(x=allFeats.apoe, y=factor(allFeats.class), ntree=1000, importance=TRUE)
	oob.pVals.all <- curRF.all$votes[,'Disease']
	oob.pVals.all.class <- sub('[0-9]+$', '', sub('^.*\\.', '', names(oob.pVals.all)))
	pdf('clustAndKaranAndAPOE.fake.pdf')
	oob.predStats.all <- getAUCandP(oob.pVals.all, oob.pVals.all.class != 'Normal',
		titleStr='Cluster and Karan and APOE, NON-XVALIDATED OOB')
	dev.off()
	resList[['KaranAndClusterAndAPOE.fake']] <- list(curRF=curRF.all, oob.predStats=unlist(oob.predStats.all))
	
	#Conclusion: mProbes seems to do a good job whith rFerns and Sturges
	nFolds <- 16
	folds <- cut(1:ncol(exprData), nFolds)
	uFolds <- unique(folds)
	folds <- folds[sample.int(length(folds))]

	#NOTE FOR AUC: KSVM builds a probabilistic model using 3-fold xvalidation to fit parameters on a Laplacian model
	#   This seems to be much more powerful that using the votes from the randomForest
	allPreds <- foreach (testFold = uFolds) %dopar% {
		curTestBool <- folds == testFold
		trainData.exp <- t(exprData[,!curTestBool])
		testData.exp <- t(exprData[,curTestBool])
		trainData.clust <- cFeats[rownames(trainData.exp),]
		testData.clust <- cFeats[rownames(testData.exp),]
	
		trainData.apoe <- dataMatrix.apoe.noGender[sub('\\..*$', '', rownames(trainData.exp)),,drop=FALSE]
		rownames(trainData.apoe) <- rownames(trainData.exp)
		testData.apoe <- dataMatrix.apoe.noGender[sub('\\..*$', '', rownames(testData.exp)),,drop=FALSE]
		rownames(testData.apoe) <- rownames(testData.exp)
		
		trainData.class <- classes.bin[rownames(trainData.exp)]
		testData.class <- classes.bin[rownames(testData.exp)]
		
		preds <- makeKaransFeatures(trainSet=t(trainData.exp), trainClass=trainData.class, testSet=testData.exp, genesets=genesets, consensus=FALSE, makePred=FALSE)
			
		#Concatenate the features.
		#  mProbes since cluster data is so large
		kFeats <- trainData.exp[, preds$bestGenes, drop=FALSE]
		rownames(kFeats) <- sub( '_[0-9]+\\.', '.', sub('^[^_]*_', '', rownames(kFeats)))
		allFeats <- cbind(kFeats, trainData.clust)
		allFeats.mProbes.fdr.ferns <- mProbes(cbind(allFeats, data.frame(class=trainData.class)), ntree = 10000, mcWorkers=1)
		ch <- hist(allFeats.mProbes.fdr.ferns) #Sturges' formula
		cutoff <- ch$breaks[length(ch$breaks)-1]
		keepFeats <- names(allFeats.mProbes.fdr.ferns)[allFeats.mProbes.fdr.ferns < cutoff]
		
		#Make RF and SVM Predictions
		x <- allFeats[,keepFeats]
		y <- factor(trainData.class)
		
		trControl <- trainControl(method='cv',number=10, classProbs = TRUE, summaryFunction=twoClassSummary)
		rf.reg <- randomForest(x=allFeats.apoe, y=factor(trainData.class), ntree=1000, importance=FALSE)
		rf.caret <- train(x=x, y=y, method='rf', metric='ROC', TuneLength=3, trControl=trControl)
		svm.caret <- train(x=x, y=y, method='svmPoly', metric='ROC', TuneLength=3, trControl=trControl)
		svm.rad.caret <- train(x=x, y=y, method='svmRadial', metric='ROC', TuneLength=3, trControl=trControl)

		#Predict class probabilities (i.e. 'certainty' scores)
		all.test <- cbind(testData.exp[, preds$bestGenes, drop=FALSE], testData.clust)
		x <- all.test[,keepFeats]
		curPreds.rf <- predict(rf.reg, x, "prob")
		curPreds.rf.caret <- predict(rf.caret, x, "prob")
		curPreds.svm.caret <- predict(svm.caret, x, "prob")
		curPreds.svm.rad.caret <- predict(svm.rad.caret, x, "prob")
		rownames(curPreds.svm.caret) <- rownames(curPreds.svm.rad.caret) <- rownames(x)

		#Make APOE alone predictions

		#Concatenate with APOE
		#Make RandomForest and SVM Predictions
		x <- cbind(allFeats[,keepFeats], trainData.apoe)
		y <- factor(trainData.class)
		
		trControl <- trainControl(method='cv',number=10, classProbs = TRUE, summaryFunction=twoClassSummary)
		rf.reg <- randomForest(x=allFeats.apoe, y=factor(trainData.class), ntree=1000, importance=FALSE)
		rf.caret <- train(x=x, y=y, method='rf', metric='ROC', TuneLength=3, trControl=trControl)
		svm.caret <- train(x=x, y=y, method='svmPoly', metric='ROC', TuneLength=3, trControl=trControl)
		svm.rad.caret <- train(x=x, y=y, method='svmRadial', metric='ROC', TuneLength=3, trControl=trControl)

		#Predict class probabilities (i.e. 'certainty' scores)
		all.test <- cbind(testData.exp[, preds$bestGenes, drop=FALSE], testData.clust)
		x <- cbind(all.test[,keepFeats], testData.apoe)
		curPreds.apoe.rf <- predict(rf.reg, x, "prob")
		curPreds.apoe.rf.caret <- predict(rf.caret, x, "prob")
		curPreds.apoe.svm.caret <- predict(svm.caret, x, "prob")
		curPreds.apoe.svm.rad.caret <- predict(svm.rad.caret, x, "prob")
		rownames(curPreds.apoe.svm.caret) <- rownames(curPreds.apoe.svm.rad.caret) <- rownames(x)

		#Save the results
		curPredList <- list(curPreds.rf=curPreds.rf, curPreds.rf.caret=curPreds.rf.caret, curPreds.svm.caret=curPreds.svm.caret, curPreds.svm.rad.caret=curPreds.svm.rad.caret,
			curPreds.apoe.rf=curPreds.apoe.rf, curPreds.apoe.rf.caret=curPreds.apoe.rf.caret, curPreds.apoe.svm.caret=curPreds.apoe.svm.caret, curPreds.apoe.svm.rad.caret=curPreds.apoe.svm.rad.caret)
		curIdx <- which(uFolds == testFold)
		saveFile <- paste('addAPOE.fold', curIdx, 'RData', sep='.')
		counter <- 0
		newSaveFile <- sub('RData', paste('rep', counter, 'RData', sep='.'), saveFile)
		while (file.exists(newSaveFile)) {
			counter <- counter + 1
			newSaveFile <- sub('RData', paste('rep', counter, 'RData', sep='.'), saveFile)
		}
		save(curPredList, file=newSaveFile)
		
		return(curPredList)
	} # allPreds <- foreach (testFold = uFolds) %dopar% {
	
	#Tally the results
	curNames <- names(allPreds[[1]])
	allPredList <- list()
	allStats <- NULL
	outFile <- paste('APOE.rocs', Sys.Date(), 'pdf', sep='.')
	pdf(outFile)
	for (curName in curNames) {
		allPredList[[curName]] <- foreach(x = allPreds, .combine=rbind) %do% {x[[curName]]}
		
		oob.pVals.all <- allPredList[[curName]][,'Disease']
		names(oob.pVals.all) <- rownames(allPredList[[curName]])
		oob.pVals.all.class <- sub('[0-9]+$', '', sub('^.*\\.', '', names(oob.pVals.all)))
		oob.predStats.all <- getAUCandP(oob.pVals.all, oob.pVals.all.class != 'Normal',
			titleStr=curName)
		allStats <- rbind(allStats, unlist(oob.predStats.all))
	}
	rownames(allStats) <- curNames
	dev.off()
	
	write.csv(allStats, file=sub('\\.pdf', '.csv', outFile))
}

#' Make predictions with and without the APOE status
#'
#' @param nFolds The number of folds. -1 means Leave-One-Out. (DEFAULT: 15)
#' @param upSample Set to true to try upsampling. (DEFAULT: FALSE)
#' @param useTri Set to true to also make trinary predictions. Overrides upSample. (DEFAULT: FALSE)
#' @export 
#' @usage allAPOEpreds <- MultipleNfoldWapoe(nFolds = 15, upSample = FALSE, useTri = FALSE)
MultipleNfoldWapoe <- function(nFolds = 15, upSample = FALSE, useTri = FALSE) {
	resList <- list()		
	allData <- read.csv('Clinical.AD080415.csv')

	if (useTri == TRUE) { upSamples <- FALSE }

	#Turn it into a data matrix
	expNames <- allData$LabLabel
	expClass <- allData$Sample
	binClass <- as.character(expClass)
	binClass[binClass != 'Normal'] <- 'Disease'
	#Gender, Abeta, ptau
	APOEGen <- as.character(allData$APOEGen)
	APOEcount <- sapply(strsplit(APOEGen, ','), function(x) {sum(x == 4)})
	Gender <- allData$Gender
	dataMatrix <- data.frame(APOEcount=APOEcount, Gender=Gender)
	rownames(dataMatrix) <- expNames
	rownames(dataMatrix) <- gsub('[#-/]', '', rownames(dataMatrix))
	names(binClass) <- rownames(dataMatrix)
	dataMatrix.apoe.noGender <- dataMatrix[,'APOEcount',drop=FALSE]
	
	#Filter to include only the correct 58
	fsScores <- getFIRMAandGENE(dataSetName='ALZ.10serum.58', subSet=NULL, recreate=FALSE, chipType="HuEx-1_0-st-v2", tags=NULL, loadFIRMA=FALSE)
	exonData <- fsScores$geneExpr.vals #get(load('./intermediate/HIVRAD.HuEx-1_0-st-v2.coreR2,A20070914,EP..geneExpr.RData')[1])
	exprData <- exonToGene(exonData)
	exprNames <- sub('_.*$', '', sub('^[^_]*_', '', colnames(exprData)))
	colnames(exprData) <- exprNames

	dataMatrix.apoe <- dataMatrix[exprNames,]
	binClass.apoe <- binClass[exprNames]

	expLUT <- csfLookup()
	classes <- expLUT[match(exprNames, names(expLUT))]
	classes <- sub('[0-9]+$', '', classes)
	colnames(exprData) <- paste(colnames(exprData), classes, sep='.')
	classes.bin <- rep('Normal', ncol(exprData))
	names(classes.bin) <- colnames(exprData)
	classes.bin[!grepl('Normal$', names(classes.bin))] <- 'Disease'

	genesets <- getGenesets(geneNames=rownames(exprData), loadExtra=TRUE, getClusterStack=FALSE)	

	clusters <- loadClusters(fName = 'inclusionMatrix.232exp.numeric.tsv')
	cFeats <- clusters$clusterScores
	rownames(cFeats) <- sub('[0-9]+$', '', rownames(cFeats))
	rownames(cFeats) <- sub('X2', '2', rownames(cFeats))
	
	#Conclusion: mProbes seems to do a good job whith rFerns and Sturges
	if (nFolds > ncol(exprData) | nFolds <= 0) {
		nFolds <- ncol(exprData)
		folds <- 1:ncol(exprData)
	} else {
		folds <- cut(1:ncol(exprData), nFolds)
	}
	uFolds <- unique(folds)
	folds <- folds[sample.int(length(folds))]

	#NOTE FOR AUC: KSVM builds a probabilistic model using 3-fold xvalidation to fit parameters on a Laplacian model
	#   This seems to be much more powerful that using the votes from the randomForest
	allPreds <- foreach (testFold = uFolds) %dopar% {
		curTestBool <- folds == testFold
		trainData.exp <- t(exprData[,!curTestBool,drop=FALSE])
		testData.exp <- t(exprData[,curTestBool,drop=FALSE])
		trainData.clust <- cFeats[rownames(trainData.exp),,drop=FALSE]
		testData.clust <- cFeats[rownames(testData.exp),,drop=FALSE]
	
		trainData.apoe <- dataMatrix.apoe.noGender[sub('\\..*$', '', rownames(trainData.exp)),,drop=FALSE]
		rownames(trainData.apoe) <- rownames(trainData.exp)
		testData.apoe <- dataMatrix.apoe.noGender[sub('\\..*$', '', rownames(testData.exp)),,drop=FALSE]
		rownames(testData.apoe) <- rownames(testData.exp)
		
		trainData.class <- classes.bin[rownames(trainData.exp)]
		testData.class <- classes.bin[rownames(testData.exp)]
		
		if (upSample == TRUE) {
			#Simplest method, just double the normals
			toDouble <- which(trainData.class == 'Normal')
			newList <- c(rownames(trainData.exp), names(toDouble))
			
			trainData.class <- trainData.class[newList]
			trainData.exp <- trainData.exp[newList,,drop=FALSE]
			trainData.clust <- trainData.clust[newList,,drop=FALSE]
			trainData.apoe <- trainData.apoe[newList,,drop=FALSE]
			
			names(trainData.class) <- make.names(names(trainData.class), unique=TRUE)
			rownames(trainData.exp) <- make.names(rownames(trainData.exp), unique=TRUE)
			rownames(trainData.clust) <- make.names(rownames(trainData.clust), unique=TRUE)
		}
		
		if (useTri == TRUE) {
			preds.tri <- list()
			comps <- list(c('MCI$', 'preMCI$'), c('preMCI$','Normal$'), c('MCI$','Normal$'))
			for (comp in comps) {
				incBool <- grepl(comp[1], rownames(trainData.exp)) | grepl(comp[2], rownames(trainData.exp))
				curTrain <- t(trainData.exp[incBool,,drop=FALSE])
				curClass <- sub('^.*\\.', '', colnames(curTrain))
				curClass[curClass == sub('.$','',comp[1])] <- 'Disease'
				preds.tri[[paste(comp, collapse='.')]] <- makeKaransFeatures(trainSet=curTrain, trainClass=curClass, testSet=testData.exp, genesets=genesets, consensus=FALSE, makePred=FALSE)$bestGenes
			}
			
			preds <- list()
			preds$bestGenes <- unique(unlist(lapply(preds.tri, function(x) {x})))
		} else {
			preds <- makeKaransFeatures(trainSet=t(trainData.exp), trainClass=trainData.class, testSet=testData.exp, genesets=genesets, consensus=FALSE, makePred=FALSE)
		}
		
		#Concatenate the features.
		#  mProbes since cluster data is so large
		kFeats <- trainData.exp[, preds$bestGenes, drop=FALSE]
		rownames(kFeats) <- sub( '_[0-9]+\\.', '.', sub('^[^_]*_', '', rownames(kFeats)))
		allFeats <- cbind(kFeats, trainData.clust)
		if (useTri == TRUE) {
			curFeats <- sub('^.*\\.', '', rownames(allFeats))
			allFeats.mProbes.fdr.ferns <- mProbes(cbind(allFeats, data.frame(class=curFeats)), ntree = 10000, mcWorkers=1)
			trainData.class <- curFeats
		} else {
			#really FWER not FDR.
			allFeats.mProbes.fdr.ferns <- mProbes(cbind(allFeats, data.frame(class=trainData.class)), ntree = 10000, mcWorkers=16)
		}
		ch <- hist(allFeats.mProbes.fdr.ferns) #Sturges' formula
		cutoff <- ch$breaks[length(ch$breaks)-1]
		keepFeats <- names(allFeats.mProbes.fdr.ferns)[allFeats.mProbes.fdr.ferns < cutoff]
		
		#Make RF and SVM Predictions
		x <- allFeats[,keepFeats,drop=FALSE]
		y <- factor(trainData.class)
		
		if (useTri == TRUE) {
			trControl <- trainControl(method='cv',number=10, classProbs = TRUE, allowParallel=FALSE)
			metric <- 'Kappa'
		} else {
			trControl <- trainControl(method='cv',number=10, classProbs = TRUE, summaryFunction=twoClassSummary, allowParallel=FALSE)
			metric <- 'ROC'
		}
		rf.reg <- randomForest(x=x, y=y, ntree=1000, importance=FALSE)
		rf.caret <- train(x=x, y=y, method='rf', metric=metric, TuneLength=3, trControl=trControl)
		svm.caret <- train(x=x, y=y, method='svmPoly', metric=metric, TuneLength=3, trControl=trControl)
		svm.rad.caret <- train(x=x, y=y, method='svmRadial', metric=metric, TuneLength=3, trControl=trControl)

		#Predict class probabilities (i.e. 'certainty' scores)
		all.test <- cbind(testData.exp[, preds$bestGenes, drop=FALSE], testData.clust)
		x <- all.test[,keepFeats,drop=FALSE]
		curPreds.rf <- predict(rf.reg, x, "prob")
		curPreds.rf.caret <- predict(rf.caret, x, "prob")
		curPreds.svm.caret <- predict(svm.caret, x, "prob")
		curPreds.svm.rad.caret <- predict(svm.rad.caret, x, "prob")
		rownames(curPreds.svm.caret) <- rownames(curPreds.svm.rad.caret) <- rownames(x)

		#Make APOE alone predictions
		x <- trainData.apoe
		y <- factor(trainData.class)

		#trControl <- trainControl(method='cv',number=10, classProbs = TRUE, summaryFunction=twoClassSummary)
		rf.reg <- randomForest(x=x, y=y, ntree=1000, importance=FALSE)
		rf.caret <- train(x=x, y=y, method='rf', metric=metric, TuneLength=3, trControl=trControl)
		svm.caret <- train(x=x, y=y, method='svmPoly', metric=metric, TuneLength=3, trControl=trControl)
		svm.rad.caret <- train(x=x, y=y, method='svmRadial', metric=metric, TuneLength=3, trControl=trControl)

		#Predict class probabilities (i.e. 'certainty' scores)
		x <- testData.apoe
		curPreds.apoeOnly.rf <- predict(rf.reg, x, "prob")
		curPreds.apoeOnly.rf.caret <- predict(rf.caret, x, "prob")
		curPreds.apoeOnly.svm.caret <- predict(svm.caret, x, "prob")
		curPreds.apoeOnly.svm.rad.caret <- predict(svm.rad.caret, x, "prob")
		rownames(curPreds.apoeOnly.svm.caret) <- rownames(curPreds.apoeOnly.svm.rad.caret) <- rownames(x)


		#Concatenate with APOE
		#Make RandomForest and SVM Predictions
		x <- cbind(allFeats[,keepFeats,drop=FALSE], trainData.apoe)
		y <- factor(trainData.class)
		
		rf.reg <- randomForest(x=x, y=y, ntree=1000, importance=FALSE)
		rf.caret <- train(x=x, y=y, method='rf', metric=metric, TuneLength=3, trControl=trControl)
		svm.caret <- train(x=x, y=y, method='svmPoly', metric=metric, TuneLength=3, trControl=trControl)
		svm.rad.caret <- train(x=x, y=y, method='svmRadial', metric=metric, TuneLength=3, trControl=trControl)

		#Predict class probabilities (i.e. 'certainty' scores)
		all.test <- cbind(testData.exp[, preds$bestGenes, drop=FALSE], testData.clust)
		x <- cbind(all.test[,keepFeats,drop=FALSE], testData.apoe)
		curPreds.apoe.rf <- predict(rf.reg, x, "prob")
		curPreds.apoe.rf.caret <- predict(rf.caret, x, "prob")
		curPreds.apoe.svm.caret <- predict(svm.caret, x, "prob")
		curPreds.apoe.svm.rad.caret <- predict(svm.rad.caret, x, "prob")
		rownames(curPreds.apoe.svm.caret) <- rownames(curPreds.apoe.svm.rad.caret) <- rownames(x)

		#Save the results
		curPredList <- list(curPreds.rf=curPreds.rf, curPreds.rf.caret=curPreds.rf.caret, curPreds.svm.caret=curPreds.svm.caret, curPreds.svm.rad.caret=curPreds.svm.rad.caret,
			curPreds.apoeOnly.rf=curPreds.apoeOnly.rf, curPreds.apoeOnly.rf.caret=curPreds.apoeOnly.rf.caret, curPreds.apoeOnly.svm.caret=curPreds.apoeOnly.svm.caret, curPreds.apoeOnly.svm.rad.caret=curPreds.apoeOnly.svm.rad.caret,
			curPreds.apoe.rf=curPreds.apoe.rf, curPreds.apoe.rf.caret=curPreds.apoe.rf.caret, curPreds.apoe.svm.caret=curPreds.apoe.svm.caret, curPreds.apoe.svm.rad.caret=curPreds.apoe.svm.rad.caret)
		curIdx <- which(uFolds == testFold)
		saveFile <- paste('./intermediate/addAPOE.fold', curIdx+1, nFolds, Sys.Date(),'RData', sep='.')
		if(upSample == TRUE) {
			saveFile <- sub('\\.fold', '.upSample.fold', saveFile)
		}
		if (useTri == TRUE) {
			saveFile <- sub('\\.fold', '.tri.fold', saveFile)
		}
		counter <- 0
		newSaveFile <- sub('RData', paste('rep', counter, 'RData', sep='.'), saveFile)
		while (file.exists(newSaveFile)) {
			counter <- counter + 1
			newSaveFile <- sub('RData', paste('rep', counter, 'RData', sep='.'), saveFile)
		}
		bestGenes <- preds$bestGenes
		save(bestGenes, curPredList, file=newSaveFile)
		
		return(curPredList)
	} # allPreds <- foreach (testFold = uFolds) %dopar% {
	
	#Comment out code to reload allPreds from files.  Used for debugging
	if (FALSE) {
		path <- './intermediate/'
		pattern <- 'addAPOE.tri.fold.[0-9]+.58.2015-08-19.rep.0.RData'
		fNames <- list.files(path=path, pattern=pattern, full.names = TRUE)
		allPreds <- list()
		for (fName in fNames) {
			load(fName)
			allPreds[[length(allPreds)+1]] <- curPredList
		} #for (fName in fNames) {
	} #if (FALSE)
	
	#Tally the results
	curNames <- names(allPreds[[1]])
	allPredList <- list()
	allStats <- NULL
	saveFile <- paste('./intermediate/APOE.rocs', nFolds, Sys.Date(), 'pdf', sep='.')
	if(upSample == TRUE) {
		saveFile <- sub('\\.rocs', '.rocs.upSample', saveFile)
	}
	if (useTri == TRUE) {
		saveFile <- sub('\\.rocs', '.rocs.tri', saveFile)
	}
		
	counter <- 0
	newSaveFile <- sub('pdf', paste('rep', counter, 'pdf', sep='.'), saveFile)
	while (file.exists(newSaveFile)) {
		counter <- counter + 1
		newSaveFile <- sub('pdf', paste('rep', counter, 'pdf', sep='.'), saveFile)
	}
	pdf(newSaveFile)
	for (curName in curNames) {
		allPredList[[curName]] <- foreach(x = allPreds, .combine=rbind) %do% {x[[curName]]}

		if (useTri == TRUE) {
			#For the 'tri' prediction, analysis is more complicated.
			###!!!###
			#Accuracy and Kappa would be a good place to start.  Can a get a p-Value for Cohen's kappa?
			probabilities <- allPredList[[curName]]
			expNames <- rownames(probabilities)
			classes <- sub('[0-9]+$', '', sub('^.*\\.', '', expNames))
			classMap <- c(Normal = 0, preMCI=0.5, MCI=1.0)
			oob.predStats.all <- curStats <- getTrinaryStats(probabilities, classes, classMap=classMap)
			oob.pVals.all <- 1-probabilities[,'Normal']
			oob.predStats.2 <- unlist(getAUCandP(oob.pVals.all, classes != 'Normal', titleStr=curName))
			oob.predStats.all <- c(oob.predStats.2, oob.predStats.all)
		} else {
			oob.pVals.all <- allPredList[[curName]][,'Disease']
			names(oob.pVals.all) <- rownames(allPredList[[curName]])
			oob.pVals.all.class <- sub('[0-9]+$', '', sub('^.*\\.', '', names(oob.pVals.all)))
			
			oob.predStats.all <- getAUCandP(oob.pVals.all, oob.pVals.all.class != 'Normal',
				titleStr=curName)
			oob.predStats.all <- unlist(oob.predStats.all)
		}
		allStats <- rbind(allStats, oob.predStats.all)
	}
	rownames(allStats) <- curNames
	dev.off()
	
	save(allStats, allPreds, allPredList, file=sub('\\.pdf', '.RData', newSaveFile))
	write.csv(allStats, file=sub('\\.pdf', '.csv', newSaveFile))
	
	return(allPredList)
}

#' Calculate statistics for trinary predictions.  
#'   For sure, accuracy, cohen's kappa (kappa2() {irr}) 
#'   Correlation should be possible if the classes have a clear rank
#'
#' @param probabilities  The number of times to run the n-fold xvalidation (DEFAULT: 'preds.xvalidateKaransAndClust.nFolds.16.2015-06.*\\.RData$')
#' @param classes  The number of times to run the n-fold xvalidation (DEFAULT: 'preds.xvalidateKaransAndClust.nFolds.16.2015-06.*\\.RData$')
#' @param classMap  An optional map for correlation, e.g. c(Normal = 0, preMCI=0.5, MCI=1.0) (DEFAULT: NULL)
#' @export 
#' @usage curStats <- getTrinaryStats(probabilities, classes, classMap=NULL)
getTrinaryStats <- function(probabilities, classes, classMap=NULL) {
	#Kappa statistics
	predictions <- apply(probabilities, 1, function(x){colnames(probabilities)[which.max(x)]})
	ratings <- data.frame(predictions=predictions, classes=classes)
	kappa.stats <- kappa2(ratings)
	kappaStat <- kappa.stats$value
	kappaStat.p <- kappa.stats$p.value
	
	#Agreement statistics
	agreement <- mean(predictions == classes)
	num <- round(agreement*nrow(ratings))
	size <- nrow(ratings)
	prob <- max(table(classes))/nrow(ratings)
	agreement.p <- pbinom(num, size=size, prob=prob, lower.tail=FALSE)
	
	res <- c(Agreement=agreement, Agreement.p=agreement.p, Kappa=kappaStat, Kappa.p=kappaStat.p)
	
	#Optional, correlation (assumes the classes are in increasing order)
	if(!is.null(classMap)) {
		class.num <- classMap[classes]
		prob.num <- classMap[predictions]
		cor.1 <- cor.test(class.num, prob.num)
		res <- c(res, Corr=cor.1$estimate, Corr.p=cor.1$p.value)
		
		prob.num.weight <- apply(probabilities, 1, function(x) {sum(classMap[names(x)]*x)})
		cor.2 <- cor.test(class.num, prob.num.weight)	
		res <- c(res, Corr.weight=cor.2$estimate, Corr.weight.p=cor.2$p.value)
	}
	
	return(res)
}


#' @param pattern  The number of times to run the n-fold xvalidation (DEFAULT: 'preds.xvalidateKaransAndClust.nFolds.16.2015-06.*\\.RData$')
#' @export 
#' @usage apoeData <- loadAPOE4()
loadAPOE4 <- function() {
	allData <- read.csv('Clinical.AD080415.csv')

	#Turn it into a data matrix
	expNames <- allData$LabLabel
	expClass <- allData$Sample
	binClass <- as.character(expClass)
	binClass[binClass != 'Normal'] <- 'Disease'
	#Gender, Abeta, ptau
	APOEGen <- as.character(allData$APOEGen)
	APOEcount <- sapply(strsplit(APOEGen, ','), function(x) {sum(x == 4)})
	Gender <- allData$Gender
	dataMatrix <- data.frame(APOEcount=APOEcount, Gender=Gender)
	rownames(dataMatrix) <- expNames
	rownames(dataMatrix) <- gsub('[#-/]', '', rownames(dataMatrix))
	names(binClass) <- rownames(dataMatrix)
	dataMatrix.apoe.noGender <- dataMatrix[,'APOEcount',drop=FALSE]
	
	return(dataMatrix.apoe.noGender)
}	

#' Run mProbes with the default parameters and Sturges selection.
#'
#' @param trainData  The data matrix with classes as the last column
#' @param mcWorkers  The number of works to use. (DEAFULT: getDoParWorkers())
#' @export 
#' @usage features <- runMprobes(trainData, mcWorkers=getDoParWorkers())
runMprobes <- function(trainData, mcWorkers=getDoParWorkers()) {
	allFeats.mProbes.fdr.ferns <- mProbes(trainData, ntree = 10000, mcWorkers=mcWorkers)
	ch <- hist(allFeats.mProbes.fdr.ferns) #Sturges' formula
	cutoff <- ch$breaks[length(ch$breaks)-1]
	keepFeats <- names(allFeats.mProbes.fdr.ferns)[allFeats.mProbes.fdr.ferns < cutoff]
	return(keepFeats)
}		

#' Reanalyze the 100 16-fold xvalidation runs
#'
#' @param pattern  The number of times to run the n-fold xvalidation (DEFAULT: 'preds.xvalidateKaransAndClust.nFolds.16.2015-06.*\\.RData$')
#' @export 
#' @usage saveFiles <- reanalyzeWithAPOE4(pattern='preds.xvalidateKaransAndClust.nFolds.16.2015-0[6-7].*\\.RData$')
reanalyzeWithAPOE4 <- function(pattern='preds.xvalidateKaransAndClust.nFolds.16.2015-0[6-7].*\\.RData$') {
	saveFiles <- list.files(path='.', pattern=pattern)
	saveFiles <- saveFiles[!grepl('karanClust', saveFiles)]
	for (curFile in saveFiles) {
		predRes <- try(addMprobesAPOE4(curFile=curFile))
		saveFiles <- c(saveFiles, predRes$saveFile)
	}
	
	return(saveFiles)
}

#' Reanalyze an element 100 16-fold cross-validation run with the new mProbes method 
#'   and with the APOE-4 data
#'
#' @param curFile  
#' @export 
#' @usage predRes <- addMprobesAPOE4(curFile='preds.xvalidateKaransAndClust.nFolds.16.2015-06-19.0.RData')
addMprobesAPOE4 <- function(curFile='preds.xvalidateKaransAndClust.nFolds.16.2015-06-19.0.RData') {
	varNames <- load(curFile)
	allPreds <- get(varNames[1])
	predResults <- get(varNames[2])
	
	#Load APOE
	apoeData <- loadAPOE4()
	
	allPredList <- list()
	for (i in 1:length(allPreds)) {
		cat(i, '\n')
		all.train.both <- allPreds[[i]]$all.train
		all.test.both <- allPreds[[i]]$all.test
		train.class <- allPreds[[i]]$train.class
		names(train.class) <- rownames(all.train.both)

		#Load APOE
		trainData.apoe <- apoeData[sub('\\..*$', '', rownames(all.train.both)),,drop=FALSE]
		rownames(trainData.apoe) <- rownames(all.train.both)
		testData.apoe <- apoeData[sub('\\..*$', '', rownames(all.test.both)),,drop=FALSE]
		rownames(testData.apoe) <- rownames(all.test.both)

		#Run mProbes
		curData <- cbind(all.train.both, data.frame(class=train.class))
		features <- runMprobes(curData)
		
		trControl <- trainControl(method='cv',number=10, classProbs = TRUE, summaryFunction=twoClassSummary)
		
		#Make SVM Predictions
		x <- all.train.both[,features,drop=FALSE]
		y <- factor(train.class)
		svm.caret <- train(x=x, y=y, method='svmPoly', metric='ROC', TuneLength=3, trControl=trControl)		
		
		x <- all.test.both[,features,drop=FALSE]
		curPreds.svm.caret <- predict(svm.caret, x, "prob")
		rownames(curPreds.svm.caret) <- rownames(x)
		allPredList[['iCAP']] <- rbind(allPredList[['iCAP']], curPreds.svm.caret)

		#Run APOE
		x <- trainData.apoe
		y <- factor(train.class)
		#svm.caret <- train(x=x, y=y, method='svmPoly', metric='ROC', TuneLength=3, trControl=trControl)		
		#c50.caret <- train(x=x, y=y, method='C5.0Tree', metric='ROC', TuneLength=3, trControl=trControl)		
		curRF <- randomForest(x=x, y=y)
		#curSVM <- ksvm(x=x, y=y, kernel='vanilladot', prob.model = TRUE)

		x <- testData.apoe
		#curPreds.c50 <- predict(c50.caret, x, "prob")
		curPreds.rf <- predict(curRF, x, "prob")
		#curPreds.svm <- predict(curSVM, x, "prob")
		rownames(curPreds.rf) <- rownames(x)
		
		allPredList[['APOE4']] <- rbind(allPredList[['APOE4']], curPreds.rf)
		
		#Run iCAP + APOE4
		x <- cbind(all.train.both[,features,drop=FALSE], trainData.apoe)
		y <- factor(train.class)
		svm.caret <- train(x=x, y=y, method='svmPoly', metric='ROC', TuneLength=3, trControl=trControl)		

		x <- cbind(all.test.both[,features,drop=FALSE], testData.apoe)
		curPreds.svm.caret <- predict(svm.caret, x, "prob")
		rownames(curPreds.svm.caret) <- rownames(x)

		allPredList[['iCAP.APOE4']] <- rbind(allPredList[['iCAP.APOE4']], curPreds.svm.caret)
	} #for (i in 1:length(allPreds)) {
	
	fName <- sub('RData', 'pdf', curFile)
	predResults <- list()
	pdf(fName)
	for(curname in names(allPredList)) {
		allPreds.reg <- allPredList[[curname]]
		allClasses.reg <- sub('[0-9]+$', '', sub('^.*\\.', '', rownames(allPreds.reg)))
		allClasses.disease <- rep(FALSE, length(allClasses.reg))
		allClasses.disease[allClasses.reg != 'Normal'] <- TRUE
		names(allClasses.disease) <- rownames(allPreds.reg)

		allPreds.dPred <- allPreds.reg[,'Disease']
		names(allPreds.dPred) <- rownames(allPreds.reg)

		predResults[[curname]] <- unlist(getAUCandP(pClass=allPreds.dPred, isClass=allClasses.disease, titleStr = curname))
	} #for(curname in names(allPredList)) {
	dev.off()
	
	saveFile <- sub('RData', 'apoe.RData', curFile)
	save(allPreds, predResults, allPredList, file=saveFile)
	
	return(list(allPreds=allPreds, predResults=predResults, saveFile=saveFile))
}

#' Reanalyze an element 100 LOO cross-validation run with the new mProbes method 
#'   and with the APOE-4 data
#'
#' @param curFile  
#' @export 
#' @usage predRes <- reanalyze58.loo(path='intermediate', pattern='addAPOE.fold.[0-9]+.58.2015-08-06.rep.0.RData')
reanalyze58.loo <- function(path='intermediate', pattern='addAPOE.fold.[0-9]+.58.2015-08-06.rep.0.RData') {
	allFiles <- list.files(path=path, pattern=pattern)
	iCAP <- NULL
	APOE <- NULL
	iCAP.APOE <- NULL
	for (curFile in allFiles) {
		tst <- load(paste(path,curFile,sep='/'))
		iCAP <- rbind(curPredList$curPreds.svm.caret, iCAP)
		APOE <- rbind(curPredList$curPreds.apoeOnly.svm.caret, APOE)
		iCAP.APOE <- rbind(curPredList$curPreds.apoe.svm.caret, iCAP.APOE)
	}
	
	loo58.stats <- list()
	subList <- list()
	outFile <- paste('iCAP.APOE.loo58', Sys.Date(), 'pdf', sep='.')
	pdf(outFile)
	for (curSetName in c('iCAP', 'APOE', 'iCAP.APOE')) {
		curSet <- get(curSetName)
		pClass <- curSet$Disease
		names(pClass) <- rownames(curSet)
		isClass <- !grepl('Normal', names(pClass))
		names(isClass) <- names(pClass)
		trash <- getAUCandP(pClass, isClass, titleStr = curSetName)
		loo58.stats[[curSetName]] <- unlist(getAUCandP(pClass, isClass, titleStr = curSetName, getMaxSpecificity = TRUE))
		
		mciBool <- !grepl('\\.preMCI$', names(pClass))
		subList[['MCI']][[curSetName]] <- rbind(unlist(getAUCandP(pClass[mciBool], isClass[mciBool])), subList[['MCI']][[curSetName]])
					
		premciBool <- !grepl('\\.MCI$', names(pClass))
		subList[['preMCI']][[curSetName]] <- rbind(unlist(getAUCandP(pClass[premciBool], isClass[premciBool])), subList[['preMCI']][[curSetName]])
	}
	dev.off()
	
	loo58.df <- do.call(rbind, loo58.stats)
	write.csv(loo58.df, file=sub('\\.pdf', '.csv', outFile))
	
	return(loo58.stats)
}

#' Reanalyze an element 100 16-fold cross-validation run with the new mProbes method 
#'   and with the APOE-4 data
#'
#' @param curFile  
#' @export 
#' @usage predRes <- reanalyze58.16x(pattern='^preds\\.xvalidateKaransAndClust\\.nFolds\\.16\\.2015.+\\.[0-9]+\\.apoe\\.RData$')
reanalyze58.16x <- function(pattern='^preds\\.xvalidateKaransAndClust\\.nFolds\\.16\\.2015.+\\.[0-9]+\\.apoe\\.RData$') {
	allFiles <- list.files(pattern=pattern)
	resList <- list()
	subList <- list(preMCI=list(), MCI=list())
	
	outFile <- paste('iCAP.APOE.100x16', Sys.Date(), 'pdf', sep='.')
	pdf(outFile)
	for (curFile in allFiles) {
		cat(which(allFiles==curFile), ', ', sep='')
		tst <- load(curFile)
		for (curSetName in c('iCAP', 'APOE4', 'iCAP.APOE4')) {
			curSet <- allPredList[[curSetName]]
			pClass <- curSet[,'Disease']
			names(pClass) <- rownames(curSet)
			isClass <- !grepl('Normal', names(pClass))
			names(isClass) <- names(pClass)
			titleStr <- paste(curSetName, which(allFiles==curFile), sep='.')
			trash <- getAUCandP(pClass, isClass, titleStr = titleStr)
			resList[[curSetName]] <- rbind(unlist(getAUCandP(pClass, isClass, titleStr = titleStr, getMaxSpecificity = TRUE)), resList[[curSetName]])
			
			mciBool <- !grepl('\\.preMCI$', names(pClass))
			subList[['MCI']][[curSetName]] <- rbind(unlist(getAUCandP(pClass[mciBool], isClass[mciBool])), subList[['MCI']][[curSetName]])
			
			premciBool <- !grepl('\\.MCI$', names(pClass))
			subList[['preMCI']][[curSetName]] <- rbind(unlist(getAUCandP(pClass[premciBool], isClass[premciBool])), subList[['preMCI']][[curSetName]])
		}	
	}
	cat('\n')
	dev.off()
	
	sortRows <- c("sensitivity","specificity","fpr","fdr", 
		"f1","agreement","p.value","mcc","mcc.p",
		"maxSpecificity","maxSpec.p","auc","auc.pVal")
	
	medians <- sapply(names(resList), function(x) {apply(resList[[x]], 2, median)})
	medians <- medians[sortRows,]

	medians.mci <- sapply(names(subList[['MCI']]), function(x) {apply(subList[['MCI']][[x]], 2, median)})
	write.csv(medians.mci, file=paste('Medians.100.MCI', Sys.Date(), 'csv', sep='.'))

	medians.premci <- sapply(names(subList[['preMCI']]), function(x) {apply(subList[['preMCI']][[x]], 2, median)})
	write.csv(medians.premci, file=paste('Medians.100.preMCI', Sys.Date(), 'csv', sep='.'))
	
	#Frac significant
	fracSigs <- list()
	for( pName in c("p.value","mcc.p","maxSpec.p","auc.pVal") ){
		fracSigs[[pName]] <- sapply(names(resList), function(x) {mean(resList[[x]][,pName] <= 0.05)})
	}
	fracSigs <- do.call(rbind, fracSigs)
	rownames(fracSigs) <- paste(rownames(fracSigs), 'fracSig', sep='.')
	
	#Now use a wilcox.test to compare iCAP.APOE4 to APOE4 and iCAP
	diffSummary <- NULL
	difList <- list()
	pValList <- list()
	ref <- resList[['iCAP.APOE4']][,sortRows]
	for (testName in c('APOE4','iCAP')) {
		tst <- resList[[testName]][,sortRows]
		difList[[testName]] <- curDif <- apply(ref, 2, median)-apply(tst, 2, median)
		pValList[[testName]] <- curPs <- sapply(colnames(ref), function(x){wilcox.test(ref[,x], tst[,x], paired=TRUE)$p.value})
		
		names(curPs) <- paste(names(curPs), 'pVal', sep='.')
		curDiff <- rbind(curDif, curPs)
		curName <- paste('iCAP.APOE4.vs', testName, rownames(curDiff), sep='.')
		rownames(curDiff) <- curName
		diffSummary <- rbind(curDiff, diffSummary)
	}
	
	#resList.df <- do.call(rbind, resList)
	#write.csv(resList.df, file=sub('\\.pdf', '.csv', outFile))
	rv <- list(resList=resList, medians=t(medians), diffSummary=diffSummary, subList=subList, medians.mci=t(medians.mci), medians.premci=t(medians.premci))
	save(rv, file=sub('\\.pdf', '.RData', outFile))
	
	write.csv(t(medians), file=sub('\\.pdf', 'medians.csv', outFile))
	write.csv(diffSummary, file=sub('\\.pdf', 'diffSummary.csv', outFile))
	
	return(rv)
}

#' Reanalyze an element 100 16-fold cross-validation run with the new mProbes method 
#'   and with the APOE-4 data
#'   Modified to try upsampling 
#'
#' @param curFile  
#' @export 
#' @usage predRes <- addMprobesAPOE4.upSample(curFile='preds.xvalidateKaransAndClust.nFolds.16.2015-06-19.0.RData')
addMprobesAPOE4.upSample <- function(curFile='preds.xvalidateKaransAndClust.nFolds.16.2015-06-19.0.RData') {
	varNames <- load(curFile)
	allPreds <- get(varNames[1])
	predResults <- get(varNames[2])
	
	#Load APOE
	apoeData <- loadAPOE4()
	
	allPredList <- list()
	for (i in 1:length(allPreds)) {
		cat(i, '\n')
		all.train.both <- allPreds[[i]]$all.train
		all.test.both <- allPreds[[i]]$all.test
		train.class <- allPreds[[i]]$train.class
		names(train.class) <- rownames(all.train.both)

		#Load APOE
		trainData.apoe <- apoeData[sub('\\..*$', '', rownames(all.train.both)),,drop=FALSE]
		rownames(trainData.apoe) <- rownames(all.train.both)
		testData.apoe <- apoeData[sub('\\..*$', '', rownames(all.test.both)),,drop=FALSE]
		rownames(testData.apoe) <- rownames(all.test.both)

		#Run mProbes
		curData <- cbind(all.train.both, data.frame(class=train.class))
		features <- runMprobes(curData)
		
		trControl <- trainControl(method='cv',number=10, classProbs = TRUE, summaryFunction=twoClassSummary)
		
		#Make SVM Predictions
		x <- all.train.both[,features,drop=FALSE]
		y <- factor(train.class)
		svm.caret <- train(x=x, y=y, method='svmPoly', metric='ROC', TuneLength=3, trControl=trControl)		
		
		x <- all.test.both[,features,drop=FALSE]
		curPreds.svm.caret <- predict(svm.caret, x, "prob")
		rownames(curPreds.svm.caret) <- rownames(x)
		allPredList[['iCAP']] <- rbind(allPredList[['iCAP']], curPreds.svm.caret)

		#Run APOE
		x <- trainData.apoe
		y <- factor(train.class)
		#svm.caret <- train(x=x, y=y, method='svmPoly', metric='ROC', TuneLength=3, trControl=trControl)		
		#c50.caret <- train(x=x, y=y, method='C5.0Tree', metric='ROC', TuneLength=3, trControl=trControl)		
		curRF <- randomForest(x=x, y=y)
		#curSVM <- ksvm(x=x, y=y, kernel='vanilladot', prob.model = TRUE)

		x <- testData.apoe
		#curPreds.c50 <- predict(c50.caret, x, "prob")
		curPreds.rf <- predict(curRF, x, "prob")
		#curPreds.svm <- predict(curSVM, x, "prob")
		rownames(curPreds.rf) <- rownames(x)
		
		allPredList[['APOE4']] <- rbind(allPredList[['APOE4']], curPreds.rf)
		
		#Run iCAP + APOE4
		x <- cbind(all.train.both[,features,drop=FALSE], trainData.apoe)
		y <- factor(train.class)
		svm.caret <- train(x=x, y=y, method='svmPoly', metric='ROC', TuneLength=3, trControl=trControl)		

		x <- cbind(all.test.both[,features,drop=FALSE], testData.apoe)
		curPreds.svm.caret <- predict(svm.caret, x, "prob")
		rownames(curPreds.svm.caret) <- rownames(x)

		allPredList[['iCAP.APOE4']] <- rbind(allPredList[['iCAP.APOE4']], curPreds.svm.caret)
	} #for (i in 1:length(allPreds)) {
	
	fName <- sub('RData', 'pdf', curFile)
	predResults <- list()
	pdf(fName)
	for(curname in names(allPredList)) {
		allPreds.reg <- allPredList[[curname]]
		allClasses.reg <- sub('[0-9]+$', '', sub('^.*\\.', '', rownames(allPreds.reg)))
		allClasses.disease <- rep(FALSE, length(allClasses.reg))
		allClasses.disease[allClasses.reg != 'Normal'] <- TRUE
		names(allClasses.disease) <- rownames(allPreds.reg)

		allPreds.dPred <- allPreds.reg[,'Disease']
		names(allPreds.dPred) <- rownames(allPreds.reg)

		predResults[[curname]] <- unlist(getAUCandP(pClass=allPreds.dPred, isClass=allClasses.disease, titleStr = curname))
	} #for(curname in names(allPredList)) {
	dev.off()
	
	saveFile <- sub('RData', 'apoe.RData', curFile)
	save(allPreds, predResults, allPredList, file=saveFile)
	
	return(list(allPreds=allPreds, predResults=predResults, saveFile=saveFile))
}
