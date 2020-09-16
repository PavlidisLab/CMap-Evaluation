# L1000-Query Bulk Export Script (Exact Dose Comparison; Hybrid-Limit)

# Load Libraries
require(NCustomLibs)
require(data.table)
require(stringi)

# Declaring Global Variables
inDir = 'Data/Subset/'
queryDir = 'L1000_Hybrid/Query/'
workDir = 'L1000_Hybrid/WorkList/'
logDir = 'L1000_Hybrid/WorkList/LOGS/'
touchstonePath = 'LINCS/GSE92742_Broad_LINCS_pert_info.txt'

options(
	stringsAsFactors = FALSE,
	warn = 1
)

# Global Preparation ------------------------------------------------------
# Creation of Output Directory
dir.create(queryDir, recursive = TRUE, showWarnings = FALSE)
dir.create(workDir, recursive = TRUE, showWarnings = FALSE)
dir.create(logDir, recursive = TRUE, showWarnings = FALSE)

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
selectionVector = (cmapMini$upCount >= 10) & (cmapMini$downCount >= 10)
cmapMini = cmapMini[selectionVector, ]
cmapMini = as.data.table(cmapMini)
cmapMini = cmapMini[trt_iname %in% touchstoneVector]
setorder(cmapMini, -sum)
rm(selectionVector, touchstoneVector)

# Sequential Loading to LINCS and Creation of Output Folders
finalBatch = ceiling(nrow(cmapMini)/30)
for (currentBatch in 1:finalBatch) {
	# Create Output Folders
	dir.create(sprintf('%sBatch_%02i', queryDir, currentBatch), recursive = TRUE, showWarnings = FALSE)
	
	# Create TSV Dump
	startIndex = 1 + (currentBatch - 1) * 30
	endIndex = min(currentBatch * 30, nrow(cmapMini))
	tempSubset = cmapMini[startIndex:endIndex, .(cell_id, trt_iname, queryName)]
	tempPath = sprintf('%sBatch_%02i.tsv', workDir, currentBatch)
	write.table(x = tempSubset, file = tempPath, quote = FALSE, sep = '\t', row.names = FALSE, col.names = TRUE)
	
	# Uploading Signature to LINCS
	for (currentIndex in startIndex:endIndex) {
		cat(sprintf('\rCurrently: %s\r', currentIndex))
		currentID = cmapMini[currentIndex, case_ID]
		currentName = cmapMini[currentIndex, sprintf('HYBRID.%s_%s_%s', cell_id, trt_iname, case_ID)]
		
		geneList = rownames(cmapFC)
		upMask = frankv(x = cmapFC[, currentID], order = -1, na.last = NA, ties.method = 'first') %in% (1:150)
		upGenes = geneList[(cmapFC[, currentID] >= 1) & upMask]
		downMask = frankv(x = cmapFC[, currentID], order = 1, na.last = NA, ties.method = 'first') %in% (1:150)
		downGenes = geneList[(cmapFC[, currentID] <= -1) & downMask]
		
		outPath = sprintf('%s/curlInstance.sh', workDir)
		QuietFileRemove(outPath)
		
		cat('curl -X POST --header "Content-Type: application/json" --header "Accept: application/json" --header "user_key: XXX" ', file = outPath, sep = '', append = FALSE)
		
		tempParam = '-d "{ \\"tool_id\\": \\"sig_gutc_tool\\", \\"data_type\\": \\"L1000\\", \\"name\\": \\"%s\\", \\"uptag-cmapfile\\": \\"TAG\\t\\t%s\\", \\"dntag-cmapfile\\": \\"TAG\\t\\t%s\\", \\"dataset\\": \\"Touchstone\\" }" '
		tempParam = sprintf(
			tempParam,
			currentName,
			paste0(upGenes, collapse = '\\t'),
			paste0(downGenes, collapse = '\\t')
		)
		cat(tempParam, file = outPath, sep = '', append = TRUE)
		cat('"https://api.clue.io/api/jobs"', file = outPath, sep = '', append = TRUE)
		
		logPath = sprintf('%s%03i_%s.LOG', logDir, currentIndex, currentName)
		
		system(sprintf('sh %s > %s 2> %s', outPath, logPath, logPath), timeout = 600)
		QuietFileRemove(outPath)
		Sys.sleep(5)
	}
	Sys.sleep(10)
}