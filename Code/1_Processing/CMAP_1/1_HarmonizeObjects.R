# CMAP: Objects (Data + Metadata) Harmonization

# Load Libraries
require(data.table)
require(NCustomLibs)

# Declaring Global Variables
dataPath = 'CMAP/FC.tsv.xz'
metadataPath = 'CMAP/Metadata.tsv'

options(
	stringsAsFactors = FALSE,
	warn = 1
)

# Global Preparation ------------------------------------------------------
# Creation of Output Directory
outDir = 'Data/CMAP/'
dir.create(outDir, recursive = TRUE, showWarnings = FALSE)

# Main Processing -------------------------------------------------------
cat('Loading Metadata.\n')
metadata = fread(metadataPath, sep = '\t', header = TRUE, data.table = FALSE)

cat('Loading Data.\n')
cmapData = XZFread(dataPath, sep = '\t', header = TRUE)
tempRowNames = cmapData[, GENE]
cmapData[, GENE := NULL]
cmapData = as.matrix(cmapData)
rownames(cmapData) = tempRowNames
colnames(cmapData) = paste0('CMAP.', colnames(cmapData))
rm(tempRowNames)

# CC Selection and Consolidation
# Note: Refering to Concentration and Duration Columns by index, since the original naming is potentially code-breaking
cat('Consolidating Conditions.\n')
metadata$instance_id = paste0('CMAP.', metadata$instance_id)
metadata$pert_cdose = sprintf('%sM', metadata[, 5])
metadata$pert_ctime = sprintf('%sh', metadata[, 6])
metadata$condition = sprintf(
	'%s %s for %s in %s',
	metadata$pert_cdose,
	metadata$cmap_name,
	metadata$pert_ctime,
	metadata$cell
)

uniqueCondition = unique(metadata$condition)
finalMetadata = lapply(uniqueCondition, function(i) {
	tempSubset = subset(metadata, condition == i)[1, ]
	return(tempSubset[, c(
		'instance_id',
		'cell',
		'condition',
		'cmap_name',
		'pert_cdose',
		'pert_ctime',
		'vehicle'
	)])
})
finalMetadata = as.data.frame(do.call('rbind', finalMetadata))
colnames(finalMetadata) = c('case_ID', 'cell_id', 'trt_cval', 'trt_iname', 'trt_cdose', 'trt_ctime', 'ctl_iname')
finalData = cmapData[, finalMetadata$case_ID]

# Saving Files
cat('Saving Files.\n')

tempPath = sprintf('%smetadata.RDS.XZ', outDir)
XZSaveRDS(obj = finalMetadata, file = tempPath)

tempPath = sprintf('%sdata.RDS.XZ', outDir)
XZSaveRDS(obj = finalData, file = tempPath)