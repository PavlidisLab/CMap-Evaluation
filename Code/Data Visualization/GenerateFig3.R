# Generate Figure 3

# Load Libraries
require(NCustomLibs)
require(data.table)
require(ggplot2)
# require(cowplot)

# Theme Redefinition
theme_set(theme_bw(base_size = 13))

# Declaring Global Variables
dataDir = 'Data/Analysis_CMAP_LINCS/Subset/'
resultDir = 'Data/Analysis_CMAP_LINCS/Results/'
outDir = 'Data/Analysis_CMAP_LINCS/Output/'
outLogPath = sprintf('%sFig_3.LOG', outDir)

options(
	stringsAsFactors = FALSE,
	warn = 1
)

# Global Preparation ------------------------------------------------------
# Creation of Output Directory
dir.create(outDir, recursive = TRUE, showWarnings = FALSE)

# Custom Function (Rank Selection)
selectRank = function(x, y, type = c('max', 'median', 'min')) {
	type = match.arg(type)
	tempRank = frankv(x = x, order = 1, na.last = 'keep', ties.method = 'dense')
	
	if (type == 'max') {
		tempResult = y[which.max(tempRank)]
	}
	
	if (type == 'min') {
		tempResult = y[which.min(tempRank)]
	}
	
	if (type == 'median') {
		tempResult = y[which(tempRank == ceiling(length(tempRank)/2))][1]
	}
	
	return(tempResult)
}

# Main Subroutine -------------------------------------------------------
# Preparation
cat('Loading Data.\n')
exactEnv = new.env()
tempPath = sprintf('%snewdose.subset.RDAT.XZ', dataDir)
load(tempPath, envir = exactEnv)

cat('Loading Pre-Calculated Results.\n')
tempPath = sprintf('%sspearman.dist.landmark.newdose.RDS.XZ', resultDir)
exactResults = readRDS(tempPath)

# Renaming Matrix Columns for Faster Manipulation
colnames(exactEnv$cmapFC) = exactEnv$cmapMetadata$new_cval
colnames(exactEnv$lincsFC) = exactEnv$lincs_l3Metadata$new_cval
colnames(exactEnv$lincsMZS) = exactEnv$lincs_l5Metadata$new_cval

# Renaming Metadata Rows for Faster Manipulation
rownames(exactEnv$cmapMetadata) = exactEnv$cmapMetadata$new_cval
rownames(exactEnv$lincs_l3Metadata) = exactEnv$lincs_l3Metadata$new_cval
rownames(exactEnv$lincs_l5Metadata) = exactEnv$lincs_l5Metadata$new_cval

# ----- Generate Plots -----
resultTypeVector = c('max', 'median', 'min')
resultTypeIDVector = c('Highest', 'Median', 'Lowest')
dataTypeVector = c('lincsFC', 'lincsMZS')
dataTypeIDVector = c('CMap 2-FC', 'CMap 2-MZS')
dataResultIDVector = c('CLFC', 'CLMZS')

# Initialize Plot PDF -----
tempPath = sprintf('%sFig_3.PDF', outDir)
pdf(file = tempPath, height = 10, width = 8)

# 1. Exact Dose -----
resultDFList = lapply(1:3, function(resultIndex) {
	tempList = lapply(1:2, function(dataIndex) {
		tempValues = exactResults[[dataResultIDVector[dataIndex]]]
		tempTarget = selectRank(x = tempValues, y = exactResults$NewCVal, type = resultTypeVector[resultIndex])
		nGenes = nrow(exactEnv$cmapFC)
		
		tempDF = data.frame(
			CMAP = exactEnv$cmapFC[, tempTarget],
			LINCS = exactEnv[[dataTypeVector[dataIndex]]][, tempTarget],
			LINCS_TYPE = rep(dataTypeIDVector[dataIndex], nGenes),
			GENE_TYPE = ifelse(as.logical(exactEnv$geneMetadata$pr_is_lm), 'Landmark', 'Inferred'),
			RESULT_TYPE = rep(resultTypeIDVector[resultIndex], nGenes)
		)
		
		# Report on Selected Combination and Spearman Correlation
		tempAllSpearman = PairwiseSpearman(tempDF$CMAP, tempDF$LINCS)
		subsetDF = subset(tempDF, GENE_TYPE == 'Landmark')
		tempSubSpearman = PairwiseSpearman(subsetDF$CMAP, subsetDF$LINCS)
		tempLINCS = switch(dataIndex, exactEnv$lincs_l3Metadata, exactEnv$lincs_l5Metadata)
		tempLINCS = switch(
			dataIndex,
			sprintf('LINCS: %s; %s\n', tempLINCS[tempTarget, 'trt_cval'], tempLINCS[tempTarget, 'ctl_cval']),
			sprintf('LINCS: %s\n', tempLINCS[tempTarget, 'trt_cval'])
		)
		cat(
			sprintf('TYPE: (%s - %s)\n', resultTypeIDVector[resultIndex], dataTypeIDVector[dataIndex]),
			sprintf('CMAP: %s\n', exactEnv$cmapMetadata[tempTarget, 'trt_cval']),
			tempLINCS,
			sprintf('All Gene Spearman: %.2f\n', tempAllSpearman),
			sprintf('Measured Gene Spearman: %.2f\n\n', tempSubSpearman),
			file = outLogPath,
			sep = '',
			append = TRUE
		)
		
		return(tempDF)
	})
	tempDF = do.call('rbind', tempList)
	return(tempDF)
})
resultDF = do.call('rbind', resultDFList)
resultDF$RESULT_TYPE = factor(resultDF$RESULT_TYPE, resultTypeIDVector)
rm(resultDFList)

tempGG = ggplot(resultDF, aes(x = LINCS, y = CMAP, colour = GENE_TYPE, alpha = GENE_TYPE)) +
	geom_point(show.legend = c('colour' = TRUE, 'alpha' = FALSE), size = 1) +
	theme(axis.text.x = element_text(angle = 45, hjust = 1),
				panel.spacing.x = unit(25, 'pt'), panel.spacing.y = unit(25, 'pt')) +
	facet_grid(rows = vars(RESULT_TYPE), cols = vars(LINCS_TYPE)) +
	scale_colour_manual(values = c('Landmark' = 'Red', 'Inferred' = 'Dark Grey'), name = 'Gene Type') +
	scale_alpha_manual(values = c('Landmark' = 1, 'Inferred' = 0.2)) +
	xlab('CMap 2 Differential Expression') +
	ylab('CMap 1 Differential Expression') +
	labs(colour = 'Gene Type')

plot(tempGG)
rm(resultDF)

# Close PDF -----
dev.off()