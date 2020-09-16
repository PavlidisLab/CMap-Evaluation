# L1000-Query Result Extraction Script (Ignore Dose Comparison)

# Load Libraries
require(NCustomLibs)
require(data.table)
require(cmapR)
require(stringi)

# Declaring Global Variables
subsetDir = 'Data/Subset/'
gctxDir = 'L1000_Ignore/'
outDir = 'Results/'
touchstonePath = 'LINCS/GSE92742_Broad_LINCS_pert_info.txt'
cmapMetadataPath = 'Data/CMAP/metadata.RDS.XZ'

options(
	stringsAsFactors = FALSE,
	warn = 1
)

# Global Preparation ------------------------------------------------------
# Creation of Output Directory
dir.create(outDir, recursive = TRUE, showWarnings = FALSE)

# Main Subroutine -------------------------------------------------------
cat('Obtaining CMAP Drug List.\n')
cmapMetadata = readRDS(cmapMetadataPath)
cmapDrugList = unique(gsub(pattern = ' ', replacement = '-', x = stri_trans_toupper(cmapMetadata$trt_iname)))
rm(cmapMetadata)

cat('Generating Reference Metadata.\n')
tempPath = sprintf('%sfinal.subset.RDAT.XZ', subsetDir)
load(tempPath)

cat('Loading LINCS Metadata.\n')
touchstoneMetadata = fread(touchstonePath, sep = '\t', header = TRUE)
touchstoneVector = touchstoneMetadata[as.logical(is_touchstone) & pert_type == 'trt_cp', stri_trans_toupper(pert_iname)]

# Preparations
rm(list = ls(pattern = '^lincs'))
rm(geneMetadata, touchstoneMetadata)
QuietGC()

# Edge-Case Handling
cmapDrugList = c(cmapDrugList, 'SCOPOLAMINE')

# Generate Mini-Metadata for CMAP, Threshold by L1000-Query Limit (10-150 Genes)
cmapMini = data.frame(
	queryName = sprintf(
		'%s_%s_%s',
		cmapMetadata$case_ID,
		cmapMetadata$cell_id,
		cmapMetadata$trt_iname
	),
	trt_cdose_um = cmapMetadata$numericDose * 1E6,
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
selectionVector = between(cmapMini$upCount, 10, 150) & between(cmapMini$downCount, 10, 150)
cmapMini = cmapMini[selectionVector, ]
cmapMini = as.data.table(cmapMini)
cmapMini = cmapMini[trt_iname %in% touchstoneVector]
setorder(cmapMini, -diffCount)
batchInfo = unlist(lapply(1:20, function(i) {rep(i, 30)}), use.names = FALSE)
batchInfo = batchInfo[1:nrow(cmapMini)]
cmapMini$batchInfo = sprintf('Batch_%02i', batchInfo)
rm(selectionVector, touchstoneVector, cmapFC, batchInfo)
QuietGC()

cat('Iterative GCTX Parsing.\n')
resultList = lapply(1:nrow(cmapMini), function(i) {
	currentCase = cmapMini[i, ]
	
	# Parsing LINCS-GCT
	tempPath = sprintf(
		'%sTouchstone/%s/T.%s.gct',
		gctxDir,
		currentCase$batchInfo,
		currentCase$queryName
	)
	
	# Sanity Check
	if (!file.exists(tempPath)) {
		cat(sprintf('MISSING CONTROL: %s - %s\n', currentCase$batchInfo, currentCase$queryName))
	}
	
	suppressMessages({tempGCT = parse.gctx(fname = tempPath)})
	
	# Sanity Check
	actualCellLine = currentCase[, cell_id]
	actualCompound = currentCase[, trt_iname]
	
	# Edge-Case Handling
	if (actualCompound == 'SCOPOLAMINE-N-OXIDE') {
		actualCompound = 'SCOPOLAMINE'
	}
	
	gctCellLine = tempGCT@cdesc$cell_id
	gctCompound = stri_trans_toupper(tempGCT@cdesc$name)
	if (gctCellLine != actualCellLine) {
		cat(sprintf('CELL LINE MISMATCH - CONTROL: %s - %s\n', currentCase$batchInfo, currentCase$queryName))
		cat(sprintf('REF: %s - ACTUAL: %s\n', actualCellLine, gctCellLine))
	}
	if (gctCompound != actualCompound) {
		cat(sprintf('COMPOUND MISMATCH - CONTROL: %s - %s\n', currentCase$batchInfo, currentCase$queryName))
		cat(sprintf('REF: %s - ACTUAL: %s\n', actualCompound, gctCompound))
	}
	
	tempDT = data.table(
		compoundName = stri_trans_toupper(tempGCT@rdesc$name),
		tauScore = as.numeric(tempGCT@mat)
	)
	tempDT = tempDT[compoundName %in% cmapDrugList]
	tempDT = tempDT[(!is.na(tauScore))]
	setorder(tempDT, -tauScore)
	tempDT = unique(x = tempDT, fromLast = FALSE, by = 'compoundName')
	tempDT[, rank := frankv(x = tauScore, order = -1, na.last = NA, ties.method = 'dense')]
	setorder(tempDT, -rank)
	sliceDT = tempDT[compoundName == actualCompound]
	
	# Populate Results (LINCS)
	tempResult = currentCase[, .(trt_iname, trt_cdose_um, cell_id, upCount, downCount, sum, ratio, diffCount)]
	tempResult[, `:=` (
		lincsTauScore = sliceDT[1, tauScore],
		multiHit = nrow(sliceDT),
		final.lincsRank = sliceDT[1, rank],
		final.lincsNDrugs = nrow(tempDT),
		norm.lincsRank = (sliceDT[1, rank]/nrow(tempDT))
	)]
	
	# Parsing CMAP-GCT
	tempPath = sprintf(
		'%sQuery/%s/H.%s.gct',
		gctxDir,
		currentCase$batchInfo,
		currentCase$queryName
	)
	
	# Sanity Check
	if (!file.exists(tempPath)) {
		cat(sprintf('MISSING QUERY: %s - %s\n', currentCase$batchInfo, currentCase$queryName))
	}
	
	suppressMessages({tempGCT = parse.gctx(fname = tempPath)})
	
	# Sanity Check
	actualQuery = currentCase[, sprintf('%s_%s_%s', cell_id, trt_iname, case_ID)]
	
	gctCellLine = tempGCT@cdesc$cell_id
	gctQuery = gsub(pattern = '[.]*$', replacement = '', x = tempGCT@cdesc$name)
	gctQuery = paste0('^', gctQuery, collapse = '')
	if (gctCellLine != actualCellLine) {
		cat(sprintf('CELL LINE MISMATCH - QUERY: %s - %s\n', currentCase$batchInfo, currentCase$queryName))
		cat(sprintf('REF: %s - ACTUAL: %s\n', actualCellLine, gctCellLine))
	}
	if (!grepl(pattern = gctQuery, x = actualQuery)) {
		cat(sprintf('COMPOUND MISMATCH - QUERY: %s - %s\n', currentCase$batchInfo, currentCase$queryName))
		cat(sprintf('REF: %s - ACTUAL: %s\n', actualCompound, gctCompound))
	}
	
	tempDT = data.table(
		compoundName = stri_trans_toupper(tempGCT@rdesc$name),
		tauScore = as.numeric(tempGCT@mat)
	)
	tempDT = tempDT[compoundName %in% cmapDrugList]
	tempDT = tempDT[(!is.na(tauScore))]
	setorder(tempDT, -tauScore)
	tempDT = unique(x = tempDT, fromLast = FALSE, by = 'compoundName')
	tempDT[, rank := frankv(x = tauScore, order = -1, na.last = NA, ties.method = 'dense')]
	setorder(tempDT, -rank)
	sliceDT = tempDT[compoundName == actualCompound]
	
	# Populate Results (CMAP)
	tempResult[, `:=` (
		cmapTauScore = sliceDT[1, tauScore],
		final.cmapRank = sliceDT[1, rank],
		final.cmapNDrugs = nrow(tempDT),
		norm.cmapRank = (sliceDT[1, rank]/nrow(tempDT))
	)]
	
	return(tempResult)
})
resultDT = do.call('rbind', resultList)
rm(resultList)

cat('Saving Results.\n')
tempPath = sprintf('%sl1000.compounds.expanded.RDS.XZ', outDir)
XZSaveRDS(obj = resultDT, file = tempPath)
rm(resultDT)