# CMAP-LINCS Data Subsetting (2)

# Load Libraries
require(NCustomLibs)
require(data.table)

# Declaring Global Variables
ioDir = 'Data/Subset/'

options(
	stringsAsFactors = FALSE,
	warn = 1
)

# Global Preparation ------------------------------------------------------
# Creation of Output Directory
dir.create(ioDir, recursive = TRUE, showWarnings = FALSE)

# Main Subroutine -------------------------------------------------------
cat('Loading Data.\n')
tempPath = sprintf('%spredose.subset.RDAT.XZ', ioDir)
load(tempPath)

# ----- Dose-Subset Step -----

cat('Subsetting By Common Dosage.\n')
# Generate Numerical Dosage Information
cmapMetadata$numericDose = as.numeric(gsub(pattern = 'M$', replacement = '', x = cmapMetadata$trt_cdose))
lincs_l3Metadata$numericDose = as.numeric(gsub(pattern = 'um$', replacement = '', x = lincs_l3Metadata$trt_cdose)) * 1E-6
lincs_l5Metadata$numericDose = as.numeric(gsub(pattern = 'um$|µM$', replacement = '', x = lincs_l5Metadata$trt_cdose)) * 1E-6

# Create New Combined-Value Column
cmapMetadata$new_cval = sprintf('%.2EM %s', cmapMetadata$numericDose, cmapMetadata$triplet)
lincs_l3Metadata$new_cval = sprintf('%.2EM %s', lincs_l3Metadata$numericDose, lincs_l3Metadata$triplet)
lincs_l5Metadata$new_cval = sprintf('%.2EM %s', lincs_l5Metadata$numericDose, lincs_l5Metadata$triplet)

# Select Common Combined-Value
commonCVal = intersect(
	intersect(cmapMetadata$new_cval, lincs_l3Metadata$new_cval),
	lincs_l5Metadata$new_cval
)

cmapMetadata = subset(cmapMetadata, new_cval %in% commonCVal)
lincs_l3Metadata = subset(lincs_l3Metadata, new_cval %in% commonCVal)
lincs_l5Metadata = subset(lincs_l5Metadata, new_cval %in% commonCVal)

cmapMetadata = as.data.table(cmapMetadata)
lincs_l3Metadata = as.data.table(lincs_l3Metadata)
lincs_l5Metadata = as.data.table(lincs_l5Metadata)

cmapMetadata = cmapMetadata[, .SD[1], by = new_cval]
lincs_l3Metadata = lincs_l3Metadata[, .SD[1], by = new_cval]
lincs_l5Metadata = lincs_l5Metadata[, .SD[1], by = new_cval]

cmapMetadata = as.data.frame(cmapMetadata)
lincs_l3Metadata = as.data.frame(lincs_l3Metadata)
lincs_l5Metadata = as.data.frame(lincs_l5Metadata)

cat('Subsetting Data.\n')

cmapFC = cmapFC[, cmapMetadata$case_ID]
lincsFC = lincsFC[, lincs_l3Metadata$case_ID]
lincsFDR = lincsFDR[, lincs_l3Metadata$case_ID]
lincsTS = lincsTS[, lincs_l3Metadata$case_ID]
lincsMZS = lincsMZS[, lincs_l5Metadata$case_ID]

cat('Saving Common-Dose Copy of Data/Metadata.\n')
tempPath = sprintf('%snewdose.subset.RDAT.XZ', ioDir)
tempList = c(
	'cmapMetadata',
	'lincs_l3Metadata',
	'lincs_l5Metadata',
	'geneMetadata',
	'cmapFC',
	'lincsFC',
	'lincsFDR',
	'lincsTS',
	'lincsMZS'
)
XZSave(list = tempList, file = tempPath, envir = .GlobalEnv)