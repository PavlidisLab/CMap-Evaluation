# CMAP-LINCS Data Subsetting (1)

# Load Libraries
require(NCustomLibs)
require(data.table)
require(stringi)

# Declaring Global Variables
cmapDir = 'Data/CMAP/'
lincs_l3Dir = 'Data/LINCS/Level3/'
lincs_l5Dir = 'Data/LINCS/Level5/'
genePath = 'LINCS/GSE92742-Level5/Gene_Info.txt'
outDir = 'Data/Subset/'

options(
	stringsAsFactors = FALSE,
	warn = 1
)

# Global Preparation ------------------------------------------------------
# Creation of Output Directory
dir.create(outDir, recursive = TRUE, showWarnings = FALSE)

# Main Subroutine -------------------------------------------------------
cat('Loading Metadata.\n')

# Load Metadata (CMAP)
tempPath = sprintf('%smetadata.RDS.XZ', cmapDir)
cmapMetadata = readRDS(tempPath)

# Load Metadata (LINCS-Level3)
tempPath = sprintf('%scondition.metadata.RDS.XZ', lincs_l3Dir)
lincs_l3Metadata = readRDS(tempPath)

# Load Metadata (LINCS-Level5)
tempPath = sprintf('%smetadata.RDS.XZ', lincs_l5Dir)
lincs_l5Metadata = readRDS(tempPath)

# Load Gene Metadata
geneMetadata = fread(input = genePath, sep = '\t', header = TRUE, data.table = FALSE)

# Correct For Drug Name Differences
cmapMetadata$old_trt_iname = cmapMetadata$trt_iname
cmapMetadata$trt_iname = gsub(pattern = ' ', replacement = '-', x = stri_trans_toupper(cmapMetadata$trt_iname))
lincs_l3Metadata$old_trt_iname = lincs_l3Metadata$trt_iname
lincs_l3Metadata$trt_iname = gsub(pattern = ' ', replacement = '-', x = stri_trans_toupper(lincs_l3Metadata$trt_iname))
lincs_l5Metadata$old_trt_iname = lincs_l5Metadata$trt_iname
lincs_l5Metadata$trt_iname = gsub(pattern = ' ', replacement = '-', x = stri_trans_toupper(lincs_l5Metadata$trt_iname))

cat('Subsetting By Similar Conditions.\n')
# NOTE: Using CellID-Drug-Time Combination
cmapMetadata$triplet = sprintf('%s %s %s', cmapMetadata$cell_id, cmapMetadata$trt_iname, cmapMetadata$trt_ctime)
lincs_l3Metadata$triplet = sprintf('%s %s %s', lincs_l3Metadata$cell_id, lincs_l3Metadata$trt_iname, lincs_l3Metadata$trt_ctime)
lincs_l5Metadata$triplet = sprintf('%s %s %s', lincs_l5Metadata$cell_id, lincs_l5Metadata$trt_iname, lincs_l5Metadata$trt_ctime)

commonTriplet = intersect(
	intersect(unique(cmapMetadata$triplet), unique(lincs_l3Metadata$triplet)),
	unique(lincs_l5Metadata$triplet)
)
cmapMetadata = subset(cmapMetadata, triplet %in% commonTriplet)
lincs_l3Metadata = subset(lincs_l3Metadata, triplet %in% commonTriplet)
lincs_l5Metadata = subset(lincs_l5Metadata, triplet %in% commonTriplet)

cat('Loading and Subsetting Data.\n')

# Load and Subset Data (CMAP)
tempPath = sprintf('%sdata.RDS.XZ', cmapDir)
cmapFC = readRDS(tempPath)
cmapFC = cmapFC[, cmapMetadata$case_ID]

# Load and Subset Data (LINCS-Level3)
tempPath = sprintf('%sfc.matrix.RDS.XZ', lincs_l3Dir)
lincsFC = readRDS(tempPath)
lincsFC = lincsFC[, lincs_l3Metadata$case_ID]

tempPath = sprintf('%sfdr.matrix.RDS.XZ', lincs_l3Dir)
lincsFDR = readRDS(tempPath)
lincsFDR = lincsFDR[, lincs_l3Metadata$case_ID]

tempPath = sprintf('%sts.matrix.RDS.XZ', lincs_l3Dir)
lincsTS = readRDS(tempPath)
lincsTS = lincsTS[, lincs_l3Metadata$case_ID]

# Load and Subset Data (LINCS-Level5)
tempPath = sprintf('%sdata.RDS.XZ', lincs_l5Dir)
lincsMZS = readRDS(tempPath)
lincsMZS = lincsMZS[, lincs_l5Metadata$case_ID]

# Subset by Common Genes & BING Genes
commonGenes = intersect(
	intersect(rownames(cmapFC), rownames(lincsTS)),
	intersect(rownames(lincsMZS), as.character(geneMetadata$pr_gene_id)[as.logical(geneMetadata$pr_is_bing)])
)
cmapFC = cmapFC[commonGenes, ]
lincsFC = lincsFC[commonGenes, ]
lincsTS = lincsTS[commonGenes, ]
lincsFDR = lincsFDR[commonGenes, ]
lincsMZS = lincsMZS[commonGenes, ]

# Synchronize NA
lincsFC[is.na(lincsTS)] = NA

# Synchronize and Subset Gene Metadata
geneMetadata = geneMetadata[match(commonGenes, geneMetadata$pr_gene_id), ]
rownames(geneMetadata) = NULL

cat('Saving Dose-Agnositc Copy of Data/Metadata.\n')
tempPath = sprintf('%spredose.subset.RDAT.XZ', outDir)
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
rm(tempList)
QuietGC()

# ----- Dose-Subset Step -----

cat('Subsetting By Strongest Dosage.\n')
# Generate Numerical Dosage Information
cmapMetadata$numericDose = as.numeric(gsub(pattern = 'M$', replacement = '', x = cmapMetadata$trt_cdose))
lincs_l3Metadata$numericDose = as.numeric(gsub(pattern = 'um$', replacement = '', x = lincs_l3Metadata$trt_cdose)) * 1E-6
lincs_l5Metadata$numericDose = as.numeric(gsub(pattern = 'um$|µM$', replacement = '', x = lincs_l5Metadata$trt_cdose)) * 1E-6

# Reorder by Dosage (Decreasing Order), Select First Entry
cmapMetadata = as.data.table(cmapMetadata)
lincs_l3Metadata = as.data.table(lincs_l3Metadata)
lincs_l5Metadata = as.data.table(lincs_l5Metadata)

setorder(cmapMetadata, -numericDose, case_ID)
setorder(lincs_l3Metadata, -numericDose, case_ID)
setorder(lincs_l5Metadata, -numericDose, case_ID)

cmapMetadata = cmapMetadata[, .SD[1], by = triplet]
lincs_l3Metadata = lincs_l3Metadata[, .SD[1], by = triplet]
lincs_l5Metadata = lincs_l5Metadata[, .SD[1], by = triplet]

cmapMetadata = as.data.frame(cmapMetadata)
lincs_l3Metadata = as.data.frame(lincs_l3Metadata)
lincs_l5Metadata = as.data.frame(lincs_l5Metadata)

cat('Subsetting Data.\n')

cmapFC = cmapFC[, cmapMetadata$case_ID]
lincsFC = lincsFC[, lincs_l3Metadata$case_ID]
lincsFDR = lincsFDR[, lincs_l3Metadata$case_ID]
lincsTS = lincsTS[, lincs_l3Metadata$case_ID]
lincsMZS = lincsMZS[, lincs_l5Metadata$case_ID]

cat('Saving Strongest-Dose Copy of Data/Metadata.\n')
tempPath = sprintf('%sfinal.subset.RDAT.XZ', outDir)
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