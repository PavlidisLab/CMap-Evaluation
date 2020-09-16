# Cross-Dataset DE Agreement (Jaccard)

# Load Libraries
require(NCustomLibs)
require(data.table)
require(stringi)

# Declaring Global Variables
inDir = 'Data/Subset/'
outDir = 'Results/'

options(
	stringsAsFactors = FALSE,
	warn = 1
)

# Global Preparation ------------------------------------------------------
# Creation of Output Directory
dir.create(outDir, recursive = TRUE, showWarnings = FALSE)

# Custom Function
CreateTernaryVector = function(inputVector, type = c('FC', 'MZS')) {
	# Ternary Vector Representation of Post-Threshold Numeric Vector
	type = match.arg(type)
	
	finalVector = switch(
		type,
		FC = (abs(inputVector) >= 1) * sign(inputVector),
		MZS = (abs(inputVector) >= 2) * sign(inputVector)
	)
	return(finalVector)
}

PairwiseJaccard = function(vectorX, vectorY) {
	# Pairwise Jaccard Index
	sum.union = sum(as.logical(vectorX) | as.logical(vectorY), na.rm = TRUE)
	
	if (sum.union == 0) {
		return(0)
	}
	
	sum.intersect = sum((vectorX * vectorY) > 0, na.rm = TRUE)
	return(sum.intersect/sum.union)
}

# Main Subroutine -------------------------------------------------------
cat('Loading Data.\n')
tempPath = sprintf('%snewdose.subset.RDAT.XZ', inDir)
load(tempPath)

# Preparations
referenceCVal = cmapMetadata$new_cval
cmapMetadata = as.data.table(cmapMetadata)
lincs_l3Metadata = as.data.table(lincs_l3Metadata)
lincs_l5Metadata = as.data.table(lincs_l5Metadata)
geneMetadata = as.data.table(geneMetadata)

finalTable = data.table(
	TRT_CVAL = cmapMetadata[, sprintf('%sum %s in %s', numericDose * 1E6, trt_iname, cell_id)],
	TRT_INAME = cmapMetadata[, trt_iname],
	TRT_DOSE = cmapMetadata[, numericDose * 1E6],
	TRT_CELL = cmapMetadata[, cell_id],
	
	# CMap vs Level 3
	L3.XD.JC.ALL = as.numeric(NA),
	L3.XD.JC.LM = as.numeric(NA),
	
	# CMap vs Level 5
	L5.XD.JC.ALL = as.numeric(NA),
	L5.XD.JC.LM = as.numeric(NA),
	
	# Level 3 vs Level 5
	CTL.XD.JC.ALL = as.numeric(NA),
	CTL.XD.JC.LM = as.numeric(NA)
)

# Agreement Calculation
geneFilter = geneMetadata[as.logical(pr_is_lm), as.character(pr_gene_id)]
for (i in 1:length(referenceCVal)) {
	currentCVal = referenceCVal[i]
	
	# LINCS-MZS Sample Modifier
	sampleSize = lincs_l5Metadata[new_cval == currentCVal, trt_sample_size]
	
	# 1. Modified Jaccard - All Genes
	cmapVector = cmapFC[, cmapMetadata[new_cval == currentCVal, case_ID]]
	lincs3Vector = lincsFC[, lincs_l3Metadata[new_cval == currentCVal, case_ID]]
	lincs5Vector = lincsMZS[, lincs_l5Metadata[new_cval == currentCVal, case_ID]] * sampleSize
	
	cmapRankVector = CreateTernaryVector(inputVector = cmapVector, type = 'FC')
	lincs3RankVector = CreateTernaryVector(inputVector = lincs3Vector, type = 'FC')
	lincs5RankVector = CreateTernaryVector(inputVector = lincs5Vector, type = 'MZS')
	
	finalTable[i, `:=`(
		L3.XD.JC.ALL = PairwiseJaccard(vectorX = cmapRankVector, lincs3RankVector),
		L5.XD.JC.ALL = PairwiseJaccard(vectorX = cmapRankVector, lincs5RankVector),
		CTL.XD.JC.ALL = PairwiseJaccard(vectorX = lincs3RankVector, lincs5RankVector)
	)]
	
	# 2. Modified Jaccard - Landmark Genes
	cmapRankVector.gf = cmapRankVector[geneFilter]
	lincs3RankVector.gf = lincs3RankVector[geneFilter]
	lincs5RankVector.gf = lincs5RankVector[geneFilter]
	
	finalTable[i, `:=`(
		L3.XD.JC.LM = PairwiseJaccard(vectorX = cmapRankVector.gf, lincs3RankVector.gf),
		L5.XD.JC.LM = PairwiseJaccard(vectorX = cmapRankVector.gf, lincs5RankVector.gf),
		CTL.XD.JC.LM = PairwiseJaccard(vectorX = lincs3RankVector.gf, lincs5RankVector.gf)
	)]
}

cat('Saving Results.\n')
tempPath = sprintf('%smulti.pairwise.mainvar.DT.RDS.XZ', outDir)
XZSaveRDS(obj = finalTable, file = tempPath)