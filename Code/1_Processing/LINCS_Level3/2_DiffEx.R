# LINCS Level 3: Differential Expression

# Load Libraries
require(cmapR)
require(data.table)
require(foreach)
require(doParallel)
require(NCustomLibs)

# Declaring Global Variables
inPath = c(
	'GSE92742' = 'LINCS/GSE92742-Level3/Level3_Data.gctx',
	'GSE70138' = 'LINCS/GSE70138-06032017-Level3/Level3_Data.gctx'
)
metadataDir = 'Data/LINCS/Level3/'

options(
	stringsAsFactors = FALSE,
	warn = 1
)

# Global Preparation ------------------------------------------------------
# Creation of Output Directory
outDir = 'Data/LINCS/Level3/'
dir.create(outDir, recursive = TRUE, showWarnings = FALSE)

# Creation of Interim RAM-backed Output Directory
ramDir = '/dev/shm/LINCS/'
dir.create(ramDir, recursive = TRUE, showWarnings = FALSE)

# Main Processing -------------------------------------------------------
cat('Loading Metadata.\n')
tempPath = sprintf('%ssample.metadata.RDS.XZ', metadataDir)
sampleDF = readRDS(tempPath)

tempPath = sprintf('%scondition.metadata.RDS.XZ', metadataDir)
matchedDF = readRDS(tempPath)

# Sample Metadata Trimming
uniqueCVals = unique(c(matchedDF$trt_cval, matchedDF$ctl_cval))
sampleDF = subset(sampleDF, pert_cval %in% uniqueCVals)
rm(uniqueCVals)

# Enable data.table optimization
sampleDF = as.data.table(sampleDF)
setkey(sampleDF, rna_centre, pert_type, pert_cval)

cat('Loading Data.\n')
referenceRow = read.gctx.ids(gctx_path = inPath['GSE92742'], dimension = 'row')
dataList = lapply(inPath, function(tempPath) {
	tempMatrix = parse.gctx(fname = tempPath)@mat[referenceRow, ]
	QuietGC()
	return(tempMatrix)
})
gctxMatrix = do.call('cbind', dataList)
rm(dataList)
QuietGC()

# Data Matrix Trimming
gctxMatrix = gctxMatrix[, sampleDF$inst_id]
QuietGC()

# ----- Differential Expression Analysis -----
# Calculate Differential Gene Expression By Condition
cat('Calculating Differential Expression.\n')

# Initial Housekeeping, Initiate Multi-Thread
.startTime = date()
currentSPID = SysGetSPID()
referenceColumn = matchedDF$case_ID
caseCount = length(matchedDF$case_ID)
geneCount = length(referenceRow)
registerDoParallel(cores = 8)

# Cleaning Cache
allPath = c(
	.TRDSPathGen(trdsDir = ramDir, spid = currentSPID, type = 'lincs_fc', idx = 1:caseCount),
	.TRDSPathGen(trdsDir = ramDir, spid = currentSPID, type = 'lincs_ts', idx = 1:caseCount),
	.TRDSPathGen(trdsDir = ramDir, spid = currentSPID, type = 'lincs_pv', idx = 1:caseCount),
	.TRDSPathGen(trdsDir = ramDir, spid = currentSPID, type = 'lincs_fdr', idx = 1:caseCount)
)
QuietFileRemove(paths = allPath, check = TRUE)

# Parallel Processing
cat('Begin DE Job: ')
tempOutput = foreach(i = 1:caseCount, .inorder = TRUE) %dopar% {
	# Status Update
	if (i %% 500 == 0) {
		cat(i, ' . ', sep = '')
	}
	
	# Obtain Matching-Sample Mapping
	currentTask = matchedDF[i, ]
	trtSamples = sampleDF[.(currentTask$rna_centre, currentTask$trt_type, currentTask$trt_cval), inst_id]
	ctlSamples = sampleDF[.(currentTask$rna_centre, currentTask$ctl_type, currentTask$ctl_cval), inst_id]
	
	# Sample-Size Control
	if (currentTask$trt_orig_sample_size != currentTask$trt_final_sample_size) {
		tempLogical = sample(x = 1:currentTask$trt_orig_sample_size, size = currentTask$trt_orig_sample_size, replace = FALSE)
		tempLogical = (tempLogical %in% 1:currentTask$trt_final_sample_size)
		trtSamples = trtSamples[tempLogical]
		rm(tempLogical)
	}
	
	if (currentTask$ctl_orig_sample_size != currentTask$ctl_final_sample_size) {
		tempLogical = sample(x = 1:currentTask$ctl_orig_sample_size, size = currentTask$ctl_orig_sample_size, replace = FALSE)
		tempLogical = (tempLogical %in% 1:currentTask$ctl_final_sample_size)
		ctlSamples = ctlSamples[tempLogical]
		rm(tempLogical)
	}
	
	# Obtain Data
	trtData = gctxMatrix[, trtSamples]
	ctlData = gctxMatrix[, ctlSamples]
	
	# Calculate: Fold Change
	tempFC = rowMeans(trtData) - rowMeans(ctlData)
	
	# Calculate: T-statistic, P-value, False Discovery Rate
	tempResult = lapply(1:geneCount, function(j) {
		tryCatch(
			{
				tempTest = t.test(x = trtData[j, ], y = ctlData[j, ], alternative = 'two.sided', paired = FALSE, var.equal = TRUE)
				tempVector = c(
					TS = as.numeric(tempTest$statistic),
					PV = tempTest$p.value
				)
				return(tempVector)
			},
			error = function(e) {
				return(c(TS = NA, PV = NA))
			}
		)
	})
	tempTS = unlist(lapply(tempResult, function(j) {return(j['TS'])}))
	tempPV = unlist(lapply(tempResult, function(j) {return(j['PV'])}))
	tempFDR = p.adjust(tempPV, method = 'fdr')
	
	# Save File
	tempPath = .TRDSPathGen(trdsDir = ramDir, spid = currentSPID, type = 'lincs_fc', idx = i)
	XZSaveRDS(obj = tempFC, file = tempPath)
	
	tempPath = .TRDSPathGen(trdsDir = ramDir, spid = currentSPID, type = 'lincs_ts', idx = i)
	XZSaveRDS(obj = tempTS, file = tempPath)
	
	tempPath = .TRDSPathGen(trdsDir = ramDir, spid = currentSPID, type = 'lincs_pv', idx = i)
	XZSaveRDS(obj = tempPV, file = tempPath)
	
	tempPath = .TRDSPathGen(trdsDir = ramDir, spid = currentSPID, type = 'lincs_fdr', idx = i)
	XZSaveRDS(obj = tempFDR, file = tempPath)
	
	# Variable Cleanup
	rm(currentTask, trtSamples, ctlSamples)
	rm(trtData, ctlData)
	rm(tempResult, tempFC, tempTS, tempPV, tempFDR)
	QuietGC()
	return(NULL)
}

# End Multi-Thread
registerDoSEQ()

# Variable Cleanup
rm(matchedDF, sampleDF, gctxMatrix, allPath, tempOutput)
QuietGC()

# Merging and Saving Matrices, RAM-Disk Cleanup
cat('\nMerging and Saving Matrices.\n')

# 1. Fold Change
tempMatrix = foreach(i = 1:caseCount, .inorder = TRUE, .combine = cbind, .maxcombine = 1000) %do% {
	tempPath = .TRDSPathGen(trdsDir = ramDir, spid = currentSPID, type = 'lincs_fc', idx = i)
	return(readRDS(tempPath))
}
colnames(tempMatrix) = referenceColumn
rownames(tempMatrix) = referenceRow

tempPath = sprintf('%sfc.matrix.RDS.XZ', outDir)
XZSaveRDS(obj = tempMatrix, file = tempPath)
rm(tempMatrix)
QuietGC()
QuietFileRemove(paths = .TRDSPathGen(trdsDir = ramDir, spid = currentSPID, type = 'lincs_fc', idx = 1:caseCount))

# 2. T-Statistic
tempMatrix = foreach(i = 1:caseCount, .inorder = TRUE, .combine = cbind, .maxcombine = 1000) %do% {
	tempPath = .TRDSPathGen(trdsDir = ramDir, spid = currentSPID, type = 'lincs_ts', idx = i)
	return(readRDS(tempPath))
}
colnames(tempMatrix) = referenceColumn
rownames(tempMatrix) = referenceRow

tempPath = sprintf('%sts.matrix.RDS.XZ', outDir)
XZSaveRDS(obj = tempMatrix, file = tempPath)
rm(tempMatrix)
QuietGC()
QuietFileRemove(paths = .TRDSPathGen(trdsDir = ramDir, spid = currentSPID, type = 'lincs_ts', idx = 1:caseCount))

# 3. P-value
tempMatrix = foreach(i = 1:caseCount, .inorder = TRUE, .combine = cbind, .maxcombine = 1000) %do% {
	tempPath = .TRDSPathGen(trdsDir = ramDir, spid = currentSPID, type = 'lincs_pv', idx = i)
	return(readRDS(tempPath))
}
colnames(tempMatrix) = referenceColumn
rownames(tempMatrix) = referenceRow

tempPath = sprintf('%spv.matrix.RDS.XZ', outDir)
XZSaveRDS(obj = tempMatrix, file = tempPath)
rm(tempMatrix)
QuietGC()
QuietFileRemove(paths = .TRDSPathGen(trdsDir = ramDir, spid = currentSPID, type = 'lincs_pv', idx = 1:caseCount))

# 4. False Discovery Rate
tempMatrix = foreach(i = 1:caseCount, .inorder = TRUE, .combine = cbind, .maxcombine = 1000) %do% {
	tempPath = .TRDSPathGen(trdsDir = ramDir, spid = currentSPID, type = 'lincs_fdr', idx = i)
	return(readRDS(tempPath))
}
colnames(tempMatrix) = referenceColumn
rownames(tempMatrix) = referenceRow

tempPath = sprintf('%sfdr.matrix.RDS.XZ', outDir)
XZSaveRDS(obj = tempMatrix, file = tempPath)
rm(tempMatrix)
QuietGC()
QuietFileRemove(paths = .TRDSPathGen(trdsDir = ramDir, spid = currentSPID, type = 'lincs_fdr', idx = 1:caseCount))

# Directory Cleanup
QuietFileRemove(ramDir)

# Print Timestamp
cat(sprintf('START TIME: %s\n', .startTime))
cat(sprintf('END TIME: %s\n\n', date()))