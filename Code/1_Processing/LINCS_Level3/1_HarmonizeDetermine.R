# LINCS Level 3: Metadata Harmonization and Baseline Determination

# Load Libraries
require(data.table)
require(stringi)
require(NCustomLibs)
require(parallel)

# Declaring Global Variables
inDir = c(
	'GSE92742' = 'LINCS/GSE92742-Level3/',
	'GSE70138' = 'LINCS/GSE70138-06032017-Level3/'
)

options(
	stringsAsFactors = FALSE,
	warn = 1
)

# Global Preparation ------------------------------------------------------
# Creation of Output Directory
outDir = 'Data/LINCS/Level3/'
dir.create(outDir, recursive = TRUE, showWarnings = FALSE)

# Custom Functions
SelectionResolver = function(conditionDF, controlDF) {
	# Output Resolving Function
	finalSampleSize = min(conditionDF$sample_size, controlDF$sample_size, 20)
	tempVector = c(
		'source_dataset' = conditionDF$source_dataset,
		'rna_centre' = conditionDF$rna_centre,
		'cell_id' = conditionDF$cell_id,
		
		'trt_type' = conditionDF$pert_type,
		'trt_cval' = conditionDF$pert_cval,
		'trt_iname' = conditionDF$pert_iname,
		'trt_cdose' = conditionDF$pert_cdose,
		'trt_ctime' = conditionDF$pert_ctime,
		'trt_orig_sample_size' = conditionDF$sample_size,
		'trt_final_sample_size' = finalSampleSize,
		
		'ctl_type' = controlDF$pert_type,
		'ctl_cval' = controlDF$pert_cval,
		'ctl_iname' = controlDF$pert_iname,
		'ctl_cdose' = controlDF$pert_cdose,
		'ctl_ctime' = controlDF$pert_ctime,
		'ctl_orig_sample_size' = controlDF$sample_size,
		'ctl_final_sample_size' = finalSampleSize
	)
	return(tempVector)
}

# Iterative Processing -------------------------------------------------------
# Container For Final Metadata (Sample and Baseline-Matched)
finalSampleMetadata = NULL
finalMatchedMetadata = NULL

for (currentIndex in 1:2) {
	currentCase = names(inDir)[currentIndex]
	currentDir = inDir[currentIndex]
	
	cat(sprintf('\nLoading Metadata - %s.\n', currentCase))
	tempPath = sprintf('%sInst_Info.txt', currentDir)
	metadata = fread(input = tempPath, sep = '\t', header = TRUE, data.table = FALSE)
	
	# ----- Generate Harmonized Sample-Level Metadata -----
	
	cat('Generating Harmonized Sample-Level Metadata.\n')
	
	# Cleanup Existing Metadata
	metadata = metadata[!grepl('^-666', metadata$pert_iname), ]
	metadata$pert_dose_unit[metadata$pert_dose == -666] = NA
	metadata$pert_dose[metadata$pert_dose == -666] = NA
	
	# Harmonization of Column Names
	if (currentIndex == 2) {
		colnames(metadata)[colnames(metadata) == 'det_plate'] = 'rna_plate'
		colnames(metadata)[colnames(metadata) == 'det_well'] = 'rna_well'
	}
	
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
	metadata$rna_centre = sapply(strsplit(x = metadata$rna_plate, split = '_'), USE.NAMES = FALSE, FUN = function(i) {return(i[[1]])})
	
	# Generate Final Sample-Level Metadata
	selectionVector = c(
		'inst_id',
		'rna_plate',
		'rna_centre',
		'pert_iname',
		'pert_type',
		'pert_dose',
		'pert_dose_unit',
		'pert_time',
		'pert_cdose',
		'pert_ctime',
		'pert_cval',
		'cell_id'
	)
	sampleMetadata = metadata[, selectionVector]
	sampleMetadata$source_dataset = rep(currentCase, nrow(metadata))
	
	# Variable Cleanup
	rm(selectionVector, metadata)
	
	# ----- Generate Baseline-Matched Condition Metadata -----
	
	cat('Generating Baseline-Matched Condition Metadata.\n')
	
	# Match Contrast-Baseline
	selectionVector = c(
		'source_dataset',
		'pert_cval',
		'rna_centre',
		'pert_type',
		'cell_id',
		'pert_iname',
		'pert_ctime',
		'pert_cdose'
	)
	focusDF = unique(sampleMetadata[, selectionVector])
	
	# Populating Sample Counts
	sampleMetadata = as.data.table(sampleMetadata)
	setkey(sampleMetadata, rna_centre, pert_type, pert_cval)
	focusDF$sample_size = unlist(mclapply(1:nrow(focusDF), mc.preschedule = TRUE, mc.cores = 4, mc.cleanup = TRUE, FUN = function(i) {
		tempCondition = focusDF[i, ]
		return(sampleMetadata[.(tempCondition$rna_centre, tempCondition$pert_type, tempCondition$pert_cval), .N, nomatch = 0])
	}), recursive = FALSE, use.names = FALSE)
	sampleMetadata = as.data.frame(sampleMetadata)
	
	# Filter out Sample Size == 1
	focusDF = subset(focusDF, sample_size > 1)
	
	# Prepare Subsets
	casesDF = subset(focusDF, grepl('^trt', focusDF$pert_type))
	controlsDF = subset(focusDF, grepl('^ctl', focusDF$pert_type))
	
	tempList = mclapply(1:nrow(casesDF), mc.preschedule = TRUE, mc.cores = 32, mc.cleanup = TRUE, FUN = function(currentRowIndex) {
		currentCondition = casesDF[currentRowIndex, ]
		currentCentre = currentCondition$rna_centre
		currentType = currentCondition$pert_type
		currentCell = currentCondition$cell_id
		currentTime = currentCondition$pert_ctime
		currentDose = currentCondition$pert_cdose
		
		controlSubset = subset(controlsDF, rna_centre == currentCentre & cell_id == currentCell & pert_ctime == currentTime)
		
		# Compound Data Resolving
		if (currentType == 'trt_cp') {
			controlSubset = subset(controlSubset, pert_type %in% c('ctl_vehicle', 'ctl_untrt'))
			
			# No Matching Entry: SKIP
			if (nrow(controlSubset) == 0) {
				return(NULL)
			}
			
			if ('DMSO' %in% controlSubset$pert_iname) {
				# Baseline: DMSO
				controlSubset = subset(controlSubset, pert_iname == 'DMSO')
				finalOutput = controlSubset[order(controlSubset$sample_size, decreasing = TRUE), ][1, ]
				return(SelectionResolver(conditionDF = currentCondition, controlDF = finalOutput))
			}
			
			if (any(c('PBS', 'H2O', 'UnTrt') %in% controlSubset$pert_iname)) {
				# Baseline: PBS, H2O, Untreated Cells
				controlSubset = subset(controlSubset, pert_iname %in% c('PBS', 'H2O', 'UnTrt'))
				finalOutput = controlSubset[order(controlSubset$sample_size, decreasing = TRUE), ][1, ]
				return(SelectionResolver(conditionDF = currentCondition, controlDF = finalOutput))
				
			} else {
				# Baseline of Last Resort
				finalOutput = controlSubset[order(controlSubset$sample_size, decreasing = TRUE), ][1, ]
				return(SelectionResolver(conditionDF = currentCondition, controlDF = finalOutput))
			}
		}
		
		# # Ligand Data Resolving
		# if (currentType == 'trt_lig') {
		# 	controlSubset = subset(controlSubset, pert_type %in% c('ctl_vehicle', 'ctl_untrt'))
		# 	
		# 	# No Matching Entry: SKIP
		# 	if (nrow(controlSubset) == 0) {
		# 		return(NULL)
		# 	}
		# 	
		# 	if (any(c('PBS', 'H2O', 'UnTrt') %in% controlSubset$pert_iname)) {
		# 		# Baseline: PBS, H2O, Untreated Cells
		# 		controlSubset = subset(controlSubset, pert_iname %in% c('PBS', 'H2O', 'UnTrt'))
		# 		finalOutput = controlSubset[order(controlSubset$sample_size, decreasing = TRUE), ][1, ]
		# 		return(SelectionResolver(conditionDF = currentCondition, controlDF = finalOutput))
		# 	}
		# 	
		# 	if ('DMSO' %in% controlSubset$pert_iname) {
		# 		# Baseline: DMSO
		# 		controlSubset = subset(controlSubset, pert_iname == 'DMSO')
		# 		finalOutput = controlSubset[order(controlSubset$sample_size, decreasing = TRUE), ][1, ]
		# 		return(SelectionResolver(conditionDF = currentCondition, controlDF = finalOutput))
		# 		
		# 	} else {
		# 		# Baseline of Last Resort
		# 		finalOutput = controlSubset[order(controlSubset$sample_size, decreasing = TRUE), ][1, ]
		# 		return(SelectionResolver(conditionDF = currentCondition, controlDF = finalOutput))
		# 	}
		# }
		# 
		# # Genetic Modification Data Resolving
		# if (currentType %in% c('trt_sh', 'trt_xpr', 'trt_oe', 'trt_oe.mut')) {
		# 	controlSubset = subset(controlSubset, pert_type %in% c('ctl_vector', 'ctl_untrt'))
		# 	
		# 	# No Matching Entry: SKIP
		# 	if (nrow(controlSubset) == 0) {
		# 		return(NULL)
		# 	}
		# 	
		# 	if ('EMPTY_VECTOR' %in% controlSubset$pert_iname) {
		# 		# Baseline: Empty Vector
		# 		controlSubset = subset(controlSubset, pert_iname == 'EMPTY_VECTOR')
		# 		
		# 		if (currentDose %in% controlSubset$pert_cdose) {
		# 			finalOutput = subset(controlSubset, pert_cdose == currentDose)[1, ]
		# 		} else {
		# 			finalOutput = controlSubset[order(controlSubset$sample_size, decreasing = TRUE), ][1, ]
		# 		}
		# 		return(SelectionResolver(conditionDF = currentCondition, controlDF = finalOutput))
		# 	}
		# 	
		# 	if ('UnTrt' %in% controlSubset$pert_iname) {
		# 		# Baseline: Untreated Cells
		# 		controlSubset = subset(controlSubset, pert_iname == 'UnTrt')
		# 		finalOutput = controlSubset[order(controlSubset$sample_size, decreasing = TRUE), ][1, ]
		# 		return(SelectionResolver(conditionDF = currentCondition, controlDF = finalOutput))
		# 		
		# 	} else {
		# 		# Baseline of Last Resort
		# 		if (currentDose %in% controlSubset$pert_cdose) {
		# 			finalOutput = subset(controlSubset, pert_cdose == currentDose)[1, ]
		# 		} else {
		# 			finalOutput = controlSubset[order(controlSubset$sample_size, decreasing = TRUE), ][1, ]
		# 		}
		# 		return(SelectionResolver(conditionDF = currentCondition, controlDF = finalOutput))
		# 	}
		# }
	})
	
	# Assemble Condition Metadata
	cat('Assemble Condition Metadata.\n')
	matchedMetadata = as.data.frame(do.call('rbind', tempList))
	
	# Append Metadata to Finalized Version
	finalSampleMetadata = rbind(finalSampleMetadata, sampleMetadata)
	finalMatchedMetadata = rbind(finalMatchedMetadata, matchedMetadata)
	
	# Variable Cleanup
	rm(currentCase, currentDir)
	rm(selectionVector, focusDF, casesDF, controlsDF, tempList)
	rm(sampleMetadata, matchedMetadata)
}

# Appending Condition ID
finalMatchedMetadata = data.frame(
	case_ID = paste0('CC.', 1:nrow(finalMatchedMetadata)),
	finalMatchedMetadata
)

# Type-Switching
finalMatchedMetadata$trt_orig_sample_size = as.numeric(finalMatchedMetadata$trt_orig_sample_size)
finalMatchedMetadata$trt_final_sample_size = as.numeric(finalMatchedMetadata$trt_final_sample_size)

finalMatchedMetadata$ctl_orig_sample_size = as.numeric(finalMatchedMetadata$ctl_orig_sample_size)
finalMatchedMetadata$ctl_final_sample_size = as.numeric(finalMatchedMetadata$ctl_final_sample_size)

# Saving Files
cat('Saving Files.\n')

# 1. Sample Metadata
tempPath = sprintf('%ssample.metadata.RDS.XZ', outDir)
XZSaveRDS(obj = finalSampleMetadata, file = tempPath)

# 2. Matched Metadata
tempPath = sprintf('%scondition.metadata.RDS.XZ', outDir)
XZSaveRDS(obj = finalMatchedMetadata, file = tempPath)