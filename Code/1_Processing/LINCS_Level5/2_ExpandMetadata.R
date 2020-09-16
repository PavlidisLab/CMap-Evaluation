# LINCS Level 5: Metadata Expansion

# Load Libraries
require(data.table)
require(stringi)
require(NCustomLibs)

# Declaring Global Variables
inDir = c(
	'GSE92742' = 'LINCS/GSE92742-Level5/',
	'GSE70138' = 'LINCS/GSE70138-Level5/'
)

options(
	stringsAsFactors = FALSE,
	warn = 1
)

# Global Preparation ------------------------------------------------------
# Creation of Output Directory
outDir = 'Data/LINCS/Level5/'
dir.create(outDir, recursive = TRUE, showWarnings = FALSE)

# Iterative Processing -------------------------------------------------------

# Container For Final Metadata
finalMetadata = NULL

# ----- Generate Harmonized Metadata -----

cat('Generating Harmonized Metadata.\n')

for (currentIndex in 1:2) {
	currentCase = names(inDir)[currentIndex]
	currentDir = inDir[currentIndex]
	
	cat(sprintf('Loading Metadata - %s.\n', currentCase))
	tempPath = sprintf('%sSig_Info.txt', currentDir)
	metadata = fread(input = tempPath, sep = '\t', header = TRUE, data.table = FALSE)
	
	# Subset For Treatment Conditions
	metadata = metadata[grepl('^trt', metadata$pert_type), ]
	
	# Filter Out Unneeded Type
	metadata = metadata[!grepl('trt_sh.css', metadata$pert_type), ]
	
	# Cleanup Existing Metadata - Part 1
	metadata = metadata[!grepl('^-666', metadata$pert_iname), ]
	
	# Generate Final Sample-Level Metadata
	selectionVector = c(
		'sig_id',
		'distil_id'
	)
	metadata = metadata[, selectionVector]
	
	# Append Metadata to Finalized Version
	finalMetadata = rbind(finalMetadata, metadata)
	
	# Variable Cleanup
	rm(currentCase, currentDir)
	rm(selectionVector, metadata)
}

# ----- Finalizing -----
# Map IDs and Rename Columns
finalMetadata$case_ID = paste0('L5.CC.', 1:nrow(finalMetadata))
selectionVector = c(
	'case_ID',
	'sig_id',
	'distil_id'
)
finalMetadata = finalMetadata[, selectionVector]

# Saving Files
cat('Saving Files.\n')

# 1. Metadata
tempPath = sprintf('%smetadata.expansion.RDS.XZ', outDir)
XZSaveRDS(obj = finalMetadata, file = tempPath)