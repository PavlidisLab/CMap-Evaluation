# Generate Table S3

# Load Libraries
require(NCustomLibs)
require(data.table)

# Declaring Global Variables
resultDir = 'Data/Analysis_CMAP_LINCS/Results/'
outDir = 'Data/Analysis_CMAP_LINCS/Output/Supplement/'

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
tempPath = sprintf('%sl1000.compounds.newdose.expanded.RDS.XZ', resultDir)
resultDT = readRDS(tempPath)

finalDT = resultDT[, .(
	'Chemical Compound' = trt_iname,
	'Concentration (µM)' = trt_cdose_um,
	'Cell Line' = cell_id
)]

tempPath = sprintf('%sTable_S3.CSV', outDir)
write.table(x = finalDT, file = tempPath, quote = TRUE, row.names = FALSE, col.names = TRUE, sep = ',')