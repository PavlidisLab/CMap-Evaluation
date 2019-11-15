# Generate Table S1

# Load Libraries
require(NCustomLibs)
require(data.table)
require(stringi)

# Declaring Global Variables
resultDir = 'Data/Analysis_CMAP_LINCS/Results/'
outDir = 'Data/Analysis_CMAP_LINCS/Output/Supplement/'
touchstonePath = 'Data/LINCS/GSE92742-MISC/GSE92742_Broad_LINCS_pert_info.txt'

options(
	stringsAsFactors = FALSE,
	warn = 1
)

# Global Preparation ------------------------------------------------------
# Creation of Output Directory
dir.create(outDir, recursive = TRUE, showWarnings = FALSE)

# Main Subroutine -------------------------------------------------------
# Preparation
cat('Loading Pre-Calculated Results.\n')
tempPath = sprintf('%sspearman.dist.all.newdose.RDS.XZ', resultDir)
resultDT = readRDS(tempPath)$metadata

cat('Loading Touchstone Information.\n')
touchstoneMetadata = fread(touchstonePath, sep = '\t', header = TRUE)
touchstoneVector = touchstoneMetadata[as.logical(is_touchstone) & pert_type == 'trt_cp', stri_trans_toupper(pert_iname)]
rm(touchstoneMetadata)

resultDT[, TYPE := '']
resultDT[TRT_DATASET == 'GSE70138', TYPE := 'Phase 2']
resultDT[TRT_DATASET == 'GSE92742', TYPE := 'Phase 1 (Non-Touchstone)']
resultDT[TRT_DATASET == 'GSE92742' & TRT_INAME %in% touchstoneVector, TYPE := 'Phase 1 (Touchstone)']

finalDT = resultDT[, .(
	'Chemical Compound' = TRT_INAME,
	'Concentration (µM)' = TRT_DOSE,
	'Cell Line' = TRT_CELL,
	'CMap 2 Release' = TYPE
)]

tempPath = sprintf('%sTable_S1.CSV', outDir)
write.table(x = finalDT, file = tempPath, quote = TRUE, row.names = FALSE, col.names = TRUE, sep = ',')