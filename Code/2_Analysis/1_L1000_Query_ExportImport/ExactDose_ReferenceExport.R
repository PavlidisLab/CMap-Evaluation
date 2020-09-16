# L1000-Query Bulk Export Script (Exact Dose Comparison; Reference List)

# Load Libraries
require(NCustomLibs)
require(data.table)
require(stringi)

# Declaring Global Variables
inDir = 'Data/Subset/'
touchstoneDir = 'L1000_Reference/Touchstone/'
workDir = 'L1000_Reference/WorkList/'
touchstonePath = 'LINCS/GSE92742_Broad_LINCS_pert_info.txt'

options(
	stringsAsFactors = FALSE,
	warn = 1
)

# Global Preparation ------------------------------------------------------
# Creation of Output Directory
dir.create(touchstoneDir, recursive = TRUE, showWarnings = FALSE)
dir.create(workDir, recursive = TRUE, showWarnings = FALSE)

# Main Subroutine -------------------------------------------------------
cat('Loading Data.\n')
tempPath = sprintf('%snewdose.subset.RDAT.XZ', inDir)
load(tempPath)

cat('Loading LINCS Metadata.\n')
touchstoneMetadata = fread(touchstonePath, sep = '\t', header = TRUE)
touchstoneVector = touchstoneMetadata[as.logical(is_touchstone) & pert_type == 'trt_cp', stri_trans_toupper(pert_iname)]

# Preparations
rm(list = ls(pattern = '^lincs'))
rm(touchstoneMetadata)
QuietGC()

# Generate Mini-Metadata for CMAP, Threshold by L1000-Query Limit (10-150 Genes)
cmapMini = data.frame(
	queryName = sprintf(
		'%s_%s_%s',
		cmapMetadata$case_ID,
		cmapMetadata$cell_id,
		cmapMetadata$trt_iname
	),
	upCount = colSums(cmapFC >= 1, na.rm = TRUE),
	downCount = colSums(cmapFC <= -1, na.rm = TRUE),
	sum = colSums(abs(cmapFC) >= 1, na.rm = TRUE)
)
cmapMini$ratio = cmapMini$upCount/cmapMini$downCount
cmapMini$diffCount = abs(cmapMini$upCount - cmapMini$downCount)
cmapMini = cbind(
	cmapMetadata[, c('case_ID', 'cell_id', 'trt_iname')],
	cmapMini
)
cmapMini = as.data.table(cmapMini)
cmapMini = cmapMini[trt_iname %in% touchstoneVector]
setorder(cmapMini, -sum)
rm(touchstoneVector)

# Sequential Loading to LINCS and Creation of Output Folders
finalBatch = ceiling(nrow(cmapMini)/30)
for (currentBatch in 1:finalBatch) {
	# Create Output Folders
	dir.create(sprintf('%sBatch_%02i', touchstoneDir, currentBatch), recursive = TRUE, showWarnings = FALSE)
	
	# Create TSV Dump
	startIndex = 1 + (currentBatch - 1) * 30
	endIndex = min(currentBatch * 30, nrow(cmapMini))
	tempSubset = cmapMini[startIndex:endIndex, .(cell_id, trt_iname, queryName)]
	tempPath = sprintf('%sBatch_%02i.tsv', workDir, currentBatch)
	write.table(x = tempSubset, file = tempPath, quote = FALSE, sep = '\t', row.names = FALSE, col.names = TRUE)
}