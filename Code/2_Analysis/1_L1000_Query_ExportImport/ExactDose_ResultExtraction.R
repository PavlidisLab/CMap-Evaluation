# L1000-Query Result Extraction Script (Exact Dose Comparison)

# Load Libraries
require(NCustomLibs)
require(data.table)
require(cmapR)
require(stringi)

# Declaring Global Variables
subsetDir = 'Data/Subset/'
gctxDirList = c(
	'Reference' = 'L1000_Reference/Touchstone/',
	# 'FC-Only' = 'L1000_EX/Query_Harmonized/',
	'Hybrid' = 'L1000_Hybrid/Query/',
	'S20' = 'L1000_S20/Query/',
	'S300' = 'L1000_S300/Query/',
	'L100' = 'L1000_L100/Query/'
)
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
tempPath = sprintf('%snewdose.subset.RDAT.XZ', subsetDir)
load(tempPath)

cat('Loading LINCS Metadata.\n')
touchstoneMetadata = fread(touchstonePath, sep = '\t', header = TRUE)
touchstoneVector = touchstoneMetadata[as.logical(is_touchstone) & pert_type == 'trt_cp', stri_trans_toupper(pert_iname)]

# Preparations
rm(list = ls(pattern = '^lincs'))
rm(geneMetadata, touchstoneMetadata)
QuietGC()

# Edge-Case Handling
cmapDrugList = c(cmapDrugList, 'SCOPOLAMINE', 'PROSTAGLANDIN-E1')

# Generate Mini-Metadata for CMAP
cmapMetadata = as.data.table(cmapMetadata)
cmapMini = data.table(
	cmapMetadata[, .(case_ID, cell_id, trt_iname)],
	trt_cdose_um = cmapMetadata[, numericDose * 1E6],
	queryName = cmapMetadata[, sprintf('%s_%s_%s', case_ID, cell_id, trt_iname)],
	
	actual.upCount = colSums(cmapFC >= 1, na.rm = TRUE),
	actual.downCount = colSums(cmapFC <= -1, na.rm = TRUE),
	actual.sum = colSums(abs(cmapFC) >= 1, na.rm = TRUE)
)
cmapMini[, `:=`(
	actual.ratio = actual.upCount/actual.downCount,
	actual.diffCount = abs(actual.upCount - actual.downCount),
	
	fixed.upCount = pmin(150, actual.upCount),
	fixed.downCount = pmin(150, actual.downCount)
)]
cmapMini[, `:=`(
	fixed.sum = fixed.upCount + fixed.downCount,
	fixed.ratio = fixed.upCount/fixed.downCount,
	fixed.diffCount = abs(fixed.upCount - fixed.downCount),
	
	REFERENCE.TAU = as.numeric(NA),
	REFERENCE.RANK = as.numeric(NA),
	REFERENCE.DRUGS = as.numeric(NA),
	
	# FC_ONLY.TAU = as.numeric(NA),
	# FC_ONLY.RANK = as.numeric(NA),
	# FC_ONLY.DRUGS = as.numeric(NA),
	
	HYBRID.TAU = as.numeric(NA),
	HYBRID.RANK = as.numeric(NA),
	HYBRID.DRUGS = as.numeric(NA),
	
	S20.TAU = as.numeric(NA),
	S20.RANK = as.numeric(NA),
	S20.DRUGS = as.numeric(NA),
	
	S300.TAU = as.numeric(NA),
	S300.RANK = as.numeric(NA),
	S300.DRUGS = as.numeric(NA),
	
	L100.TAU = as.numeric(NA),
	L100.RANK = as.numeric(NA),
	L100.DRUGS = as.numeric(NA)
)]
cmapMini = cmapMini[trt_iname %in% touchstoneVector]
rm(touchstoneVector, cmapFC)
QuietGC()

cat('Iterative GCTX Parsing.\n')
for (currentDirIndex in 1:length(gctxDirList)) {
	currentDirType = names(gctxDirList)[currentDirIndex]
	currentDir = gctxDirList[currentDirIndex]
	
	workList = list.files(path = currentDir, recursive = TRUE, no.. = TRUE)
	for (currentItem in workList) {
		# Prepare Query Name and File Location
		currentQueryName = sub(pattern = 'Batch_[0-9]+/', replacement = '', x = currentItem)
		currentQueryName = sub(pattern = '.gct', replacement = '', x = currentQueryName)
		tempPath = sprintf('%s%s', currentDir, currentItem)
		
		# Sanity Check (File Existence)
		if (!file.exists(tempPath)) {
			cat(sprintf('MISSING: %s: %s\n', currentDirType , currentQueryName))
			next
		}
		
		suppressMessages({tempGCT = parse.gctx(fname = tempPath)})
		
		# Sanity Check Preparation (Cell Line/Compound Match)
		actualCellLine = cmapMini[queryName == currentQueryName, cell_id]
		actualCompound = cmapMini[queryName == currentQueryName, trt_iname]
		tempQueryString = cmapMini[queryName == currentQueryName, sprintf('%s_%s_%s', cell_id, trt_iname, case_ID)]
		actualQuery = switch(
			currentDirType,
			# 'FC-Only' = sprintf('EX.%s', tempQueryString),
			'Hybrid' = sprintf('HYBRID.%s', tempQueryString),
			'S20' = sprintf('S20.%s', tempQueryString),
			'S300' = sprintf('S300.%s', tempQueryString),
			'L100' = sprintf('L100.%s', tempQueryString),
			NA
		)
		
		# Edge-Case Handling
		actualCompound = switch(
			actualCompound,
			'SCOPOLAMINE-N-OXIDE' = 'SCOPOLAMINE',
			'ALPROSTADIL' = 'PROSTAGLANDIN-E1',
			actualCompound
		)
		
		# Actual Sanity Check
		gctCellLine = tempGCT@cdesc$cell_id
		if (gctCellLine != actualCellLine) {
			cat(sprintf('CELL LINE MISMATCH - %s: %s\n', currentDirType, currentQueryName))
			next
		}
		
		if (currentDirType == 'Reference') {
			gctCompound = stri_trans_toupper(tempGCT@cdesc$name)
			
			if (gctCompound != actualCompound) {
				cat(sprintf('COMPOUND MISMATCH - %s: %s\n', currentDirType, currentQueryName))
				next
			}
		} else {
			gctQuery = gsub(pattern = '[.]*$', replacement = '', x = tempGCT@cdesc$name)
			gctQuery = paste0('^', gctQuery, collapse = '')
			
			if (!grepl(pattern = gctQuery, x = actualQuery)) {
				cat(sprintf('COMPOUND MISMATCH - %s: %s\n', currentDirType, currentQueryName))
				next
			}
		}
		
		# Data Parsing
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
		
		# Populate Scores
		columnNameVector = switch(
			currentDirType,
			'Reference' = c('REFERENCE.TAU', 'REFERENCE.RANK', 'REFERENCE.DRUGS'),
			# 'FC-Only' = c('FC_ONLY.TAU', 'FC_ONLY.RANK', 'FC_ONLY.DRUGS'),
			'Hybrid' = c('HYBRID.TAU', 'HYBRID.RANK', 'HYBRID.DRUGS'),
			'S20' = c('S20.TAU', 'S20.RANK', 'S20.DRUGS'),
			'S300' = c('S300.TAU', 'S300.RANK', 'S300.DRUGS'),
			'L100' = c('L100.TAU', 'L100.RANK', 'L100.DRUGS')
		)
		finalValuesList = list(
			sliceDT[1, tauScore],
			sliceDT[1, rank],
			nrow(tempDT)
		)
		cmapMini[queryName == currentQueryName, (columnNameVector) := finalValuesList]
	}
}

cat('Saving Results.\n')
tempPath = sprintf('%sl1000.compounds.newdose.expanded.RDS.XZ', outDir)
XZSaveRDS(obj = cmapMini, file = tempPath)
