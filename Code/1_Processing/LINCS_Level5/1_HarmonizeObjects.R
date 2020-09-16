# LINCS Level 5: Objects (Data + Metadata) Harmonization

# Load Libraries
require(data.table)
require(stringi)
require(NCustomLibs)
require(cmapR)

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
	
	# Regenerate Malformed Columns (Phase 1 Metadata)
	if (currentIndex == 1) {
		# Fix Weird Dose Information
		tempList = strsplit(x = metadata$pert_dose, split = '\\|')
		metadata$pert_dose = sapply(tempList, function(i) {unique(as.numeric(i))})
		rm(tempList)
	}
	
	# Generate Missing Columns (Phase 2 Metadata)
	if (currentIndex == 2) {
		# Fix Missing Dose and Dose_Unit Columns
		metadata$pert_dose = as.numeric(gsub(pattern = ' um$', replacement = '', x = metadata$pert_idose))
		metadata$pert_dose_unit = rep('um', nrow(metadata))
		
		# Fix Missing Time and Time_Unit Columns
		metadata$pert_time = as.numeric(gsub(pattern = ' h$', replacement = '', x = metadata$pert_itime))
		metadata$pert_time_unit = rep('h', nrow(metadata))
	}
	
	# Cleanup Existing Metadata - Part 2
	metadata$pert_dose_unit[metadata$pert_dose == -666] = NA
	metadata$pert_dose[metadata$pert_dose == -666] = NA
	
	# Generate Combination Identifiers
	metadata$pert_cdose = sprintf(
		'%s%s',
		metadata$pert_dose,
		metadata$pert_dose_unit
	)
	metadata$pert_cdose[grepl('^NA', metadata$pert_cdose)] = ''
	
	metadata$pert_ctime = sprintf(
		'%s%s',
		metadata$pert_time,
		metadata$pert_time_unit
	)
	
	metadata$pert_cval = sprintf(
		'%s %s for %s in %s',
		metadata$pert_cdose,
		metadata$pert_iname,
		metadata$pert_ctime,
		metadata$cell_id
	)
	metadata$pert_cval = stri_trim_left(metadata$pert_cval)
	
	# Generate "Experiment Centre" Identifiers
	metadata$rna_centre = sapply(strsplit(x = metadata$sig_id, split = '_'), USE.NAMES = FALSE, FUN = function(i) {return(i[[1]])})
	
	# Generate Sample Size Counts
	metadata$trt_sample_size = sapply(strsplit(x = metadata$distil_id, split = '\\|'), USE.NAMES = FALSE, FUN = function(i) {return(length(i))})
	
	# Generate Final Sample-Level Metadata
	selectionVector = c(
		'sig_id',
		'rna_centre',
		'pert_iname',
		'pert_type',
		'pert_cdose',
		'pert_ctime',
		'pert_cval',
		'cell_id',
		'trt_sample_size'
	)
	metadata = metadata[, selectionVector]
	metadata$source_dataset = rep(currentCase, nrow(metadata))
	
	# Append Metadata to Finalized Version
	finalMetadata = rbind(finalMetadata, metadata)
	
	# Variable Cleanup
	rm(currentCase, currentDir)
	rm(selectionVector, metadata)
}

# ----- Generate Harmonized Data -----

cat('Generating Harmonized Data.\n')

# Obtain Reference Gene IDs
tempPath = sprintf('%sLevel5_Data.gctx', inDir['GSE92742'])
referenceRow = read.gctx.ids(gctx_path = tempPath, dimension = 'row')
dataList = lapply(inDir, function(tempDir) {
	tempPath = sprintf('%sLevel5_Data.gctx', tempDir)
	tempMatrix = parse.gctx(fname = tempPath)@mat[referenceRow, ]
	QuietGC()
	return(tempMatrix)
})
gctxMatrix = do.call('cbind', dataList)
rm(dataList)
QuietGC()

# Data Matrix Trimming
gctxMatrix = gctxMatrix[, finalMetadata$sig_id]
QuietGC()

# ----- Finalizing -----
# Reset IDs and Rename Columns
finalMetadata$case_ID = paste0('L5.CC.', 1:nrow(finalMetadata))
colnames(finalMetadata) = gsub(pattern = '^pert', replacement = 'trt', x = colnames(finalMetadata))
selectionVector = c(
	'case_ID',
	'source_dataset',
	'rna_centre',
	'cell_id',
	'trt_type',
	'trt_cval',
	'trt_iname',
	'trt_cdose',
	'trt_ctime',
	'trt_sample_size'
)
finalMetadata = finalMetadata[, selectionVector]
colnames(gctxMatrix) = finalMetadata$case_ID

# Saving Files
cat('Saving Files.\n')

# 1. Metadata
tempPath = sprintf('%smetadata.RDS.XZ', outDir)
XZSaveRDS(obj = finalMetadata, file = tempPath)

# 2. Data
tempPath = sprintf('%sdata.RDS.XZ', outDir)
XZSaveRDS(obj = gctxMatrix, file = tempPath)