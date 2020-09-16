# Cross-Dataset DE Agreement (Spearman)

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

# Main Subroutine -------------------------------------------------------
cat('Loading Data.\n')
tempPath = sprintf('%snewdose.subset.RDAT.XZ', inDir)
load(tempPath)

# Preparations
referenceCVal = cmapMetadata$new_cval
cmapMetadata = as.data.table(cmapMetadata)
lincs_l3Metadata = as.data.table(lincs_l3Metadata)
lincs_l5Metadata = as.data.table(lincs_l5Metadata)

# Metadata Reformat
datasetID = lincs_l3Metadata[match(referenceCVal, new_cval), source_dataset]

newMetadata = data.table(
	REF.CVAL = referenceCVal,
	TRT_INAME = cmapMetadata[, trt_iname],
	TRT_DOSE = cmapMetadata[, as.numeric(gsub(pattern = 'M', replacement = '', x = trt_cdose)) * 1E6],
	TRT_CELL = cmapMetadata[, cell_id],
	TRT_DATASET = datasetID
)

# ----- All Genes Comparison -----
cat('All Gene Comparison.\n')
resultList = list()
resultList$NewCVal = referenceCVal
resultList$metadata = newMetadata

# 1. CMAP x LINCS (Level 3, FC)
resultList$CLFC = unlist(lapply(referenceCVal, function(i) {
	tempX = cmapFC[, cmapMetadata[new_cval == i, case_ID]]
	tempY = lincsFC[, lincs_l3Metadata[new_cval == i, case_ID]]
	return(PairwiseSpearman(tempX, tempY))
}), recursive = FALSE, use.names = FALSE)

# 2. CMAP x LINCS (Level 3, TS)
resultList$CLTS = unlist(lapply(referenceCVal, function(i) {
	tempX = cmapFC[, cmapMetadata[new_cval == i, case_ID]]
	tempY = lincsTS[, lincs_l3Metadata[new_cval == i, case_ID]]
	return(PairwiseSpearman(tempX, tempY))
}), recursive = FALSE, use.names = FALSE)

# 3. CMAP x LINCS (Level 5, MZS)
resultList$CLMZS = unlist(lapply(referenceCVal, function(i) {
	tempX = cmapFC[, cmapMetadata[new_cval == i, case_ID]]
	tempY = lincsMZS[, lincs_l5Metadata[new_cval == i, case_ID]]
	return(PairwiseSpearman(tempX, tempY))
}), recursive = FALSE, use.names = FALSE)

# 4. LINCS (Level 3, FC) x LINCS (Level 5, MZS)
resultList$LL = unlist(lapply(referenceCVal, function(i) {
	tempX = lincsFC[, lincs_l3Metadata[new_cval == i, case_ID]]
	tempY = lincsMZS[, lincs_l5Metadata[new_cval == i, case_ID]]
	return(PairwiseSpearman(tempX, tempY))
}), recursive = FALSE, use.names = FALSE)

cat('Saving Results.\n')
tempPath = sprintf('%sspearman.dist.all.newdose.RDS.XZ', outDir)
XZSaveRDS(obj = resultList, file = tempPath)
rm(resultList)

# ----- Landmark Genes Comparison -----
cat('Landmark Gene Comparison.\n')
resultList = list()
resultList$NewCVal = referenceCVal
resultList$metadata = newMetadata

# Subset Data By Gene
referenceGenes = as.character(geneMetadata$pr_gene_id[as.logical(geneMetadata$pr_is_lm)])
cmapFC = cmapFC[referenceGenes, ]
lincsFC = lincsFC[referenceGenes, ]
lincsTS = lincsTS[referenceGenes, ]
lincsMZS = lincsMZS[referenceGenes, ]
QuietGC()

# 1. CMAP x LINCS (Level 3, FC)
resultList$CLFC = unlist(lapply(referenceCVal, function(i) {
	tempX = cmapFC[, cmapMetadata[new_cval == i, case_ID]]
	tempY = lincsFC[, lincs_l3Metadata[new_cval == i, case_ID]]
	return(PairwiseSpearman(tempX, tempY))
}), recursive = FALSE, use.names = FALSE)

# 2. CMAP x LINCS (Level 3, TS)
resultList$CLTS = unlist(lapply(referenceCVal, function(i) {
	tempX = cmapFC[, cmapMetadata[new_cval == i, case_ID]]
	tempY = lincsTS[, lincs_l3Metadata[new_cval == i, case_ID]]
	return(PairwiseSpearman(tempX, tempY))
}), recursive = FALSE, use.names = FALSE)

# 3. CMAP x LINCS (Level 5, MZS)
resultList$CLMZS = unlist(lapply(referenceCVal, function(i) {
	tempX = cmapFC[, cmapMetadata[new_cval == i, case_ID]]
	tempY = lincsMZS[, lincs_l5Metadata[new_cval == i, case_ID]]
	return(PairwiseSpearman(tempX, tempY))
}), recursive = FALSE, use.names = FALSE)

# 4. LINCS (Level 3, FC) x LINCS (Level 5, MZS)
resultList$LL = unlist(lapply(referenceCVal, function(i) {
	tempX = lincsFC[, lincs_l3Metadata[new_cval == i, case_ID]]
	tempY = lincsMZS[, lincs_l5Metadata[new_cval == i, case_ID]]
	return(PairwiseSpearman(tempX, tempY))
}), recursive = FALSE, use.names = FALSE)

cat('Saving Results.\n')
tempPath = sprintf('%sspearman.dist.landmark.newdose.RDS.XZ', outDir)
XZSaveRDS(obj = resultList, file = tempPath)
rm(resultList)