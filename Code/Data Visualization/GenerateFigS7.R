# Generate Figure S7

# Load Libraries
require(NCustomLibs)
require(ggplot2)
require(cowplot)

# Declaring Global Variables
resultPath = 'Data/Analysis_CMAP_LINCS/Results/IntraCmp/EXP/Packed/linked.109drugs.allGenes.DT.RDS.XZ'
outDir = 'Data/Analysis_CMAP_LINCS/Output/'
outLogPath = sprintf('%sFig_S7.LOG', outDir)

options(
	stringsAsFactors = FALSE,
	warn = 1
)

# Global Preparation ------------------------------------------------------
# Creation of Output Directory
dir.create(outDir, recursive = TRUE, showWarnings = FALSE)

# Custom Function
AddItemToList = function(inputList, inputObject) {
	inputList[[length(inputList) + 1]] = inputObject
	return(inputList)
}

CorTest = function(vectorX, vectorY) {
	return(
		cor.test(x = vectorX, y = vectorY, alternative = 'two.sided', method = 'spearman')$p.value
	)
}

MinAndMax = function(inputVector) {
	return(c(
		min(inputVector, na.rm = TRUE), max(inputVector, na.rm = TRUE)
	))
}

# Main Subroutine -------------------------------------------------------
# Preparation
resultDT = readRDS(resultPath)
resultDT = resultDT[!(is.na(REFERENCE_RANK) | is.na(RETRIEVAL_RANK))]

# ----- Generate Plots -----

# Scatterplots -----
tempPlotList = list()

# 1. CMAP 1 Recall x CMAP Inter-DE
tempGG = ggplot(resultDT, aes(x = CROSS_CMAP, y = RETRIEVAL_RANK)) +
	geom_point(alpha = 0.2) +
	xlim(c(-0.2, 0.6)) +
	ylim(c(0, 600)) +
	xlab('Cross-Dataset DE Agreement') +
	ylab('Retrieval Rank\n(CMap 1 Signature)')
tempPlotList = AddItemToList(tempPlotList, tempGG)

# 2. CMAP 2 Recall x CMAP Inter-DE
tempGG = ggplot(resultDT, aes(x = CROSS_CMAP, y = REFERENCE_RANK)) +
	geom_point(alpha = 0.2) +
	xlim(c(-0.2, 0.6)) +
	ylim(c(0, 600)) +
	xlab('Cross-Dataset DE Agreement') +
	ylab('Retrieval Rank\n(Self-Query)')
tempPlotList = AddItemToList(tempPlotList, tempGG)

# Summary Report
cat(
	sprintf('CMAP Recall x Cross-CMAP: %.3f\n', resultDT[, PairwiseSpearman(CROSS_CMAP, RETRIEVAL_RANK)]),
	sprintf('[P-VAL]: %.3e\n\n', resultDT[, CorTest(CROSS_CMAP, RETRIEVAL_RANK)]),
	
	sprintf('LINCS Recall x Cross-CMAP: %.3f\n', resultDT[, PairwiseSpearman(CROSS_CMAP, REFERENCE_RANK)]),
	sprintf('[P-VAL]: %.3e\n\n', resultDT[, CorTest(CROSS_CMAP, REFERENCE_RANK)]),
	
	sprintf('# CMAP Counts: %s\n', resultDT[!is.na(CROSS_CMAP), .N]),
	sprintf('# LINCS Counts: %s\n', resultDT[!is.na(CROSS_CMAP), .N]),
	
	'CMAP RECALL Value Spread:\n',
	paste0(sprintf('%.3f', resultDT[, MinAndMax(RETRIEVAL_RANK)]), collapse = '; '),
	
	'\nLINCS RECALL Value Spread:\n',
	paste0(sprintf('%.3f', resultDT[, MinAndMax(REFERENCE_RANK)]), collapse = '; '),
	
	'\nCross-CMAP Value Spread:\n',
	paste0(sprintf('%.3f', resultDT[, MinAndMax(CROSS_CMAP)]), collapse = '; '),
	
	file = outLogPath,
	sep = '',
	append = TRUE
)

# Cowplot Output -----
tempCPlot = plot_grid(
	plotlist = tempPlotList,
	nrow = 1,
	ncol = 2,
	labels = LETTERS[1:2]
)
tempPath = sprintf('%sFig_S7.pdf', outDir)
save_plot(tempPath, plot = tempCPlot, ncol = 2, nrow = 1)