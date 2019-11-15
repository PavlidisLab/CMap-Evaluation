# Generate Figure S3

# Load Libraries
require(NCustomLibs)
require(data.table)
require(ggplot2)
require(cowplot)

# Declaring Global Variables
resultDir = 'Data/Analysis_CMAP_LINCS/Results/'
outDir = 'Data/Analysis_CMAP_LINCS/Output/'
outLogPath = sprintf('%sFig_S3.LOG', outDir)

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

# Main Subroutine -------------------------------------------------------
# Preparation
cat('Loading Pre-Calculated Results.\n')
tempPath = sprintf('%sl1000.compounds.expanded.RDS.XZ', resultDir)
resultDT = readRDS(tempPath)

tempPath = sprintf('%sl1000.compounds.newdose.expanded.RDS.XZ', resultDir)
newdose.resultDT = readRDS(tempPath)

# ----- Generate Plots -----

# Scatterplots -----
# 1. From 588 Conditions
tempPlotList = list()
tempDT = resultDT[, .(
	trt_iname,
	cell_id,
	# upCount,
	# downCount,
	# sum,
	# l10_ratio = log10(ratio),
	CMAP = final.cmapRank,
	LINCS = final.lincsRank
)]

# Summary Report
wilcoxCMAP = wilcox.test(formula = CMAP ~ cell_id, data = tempDT, alternative = 'two.sided')
wilcoxLINCS = wilcox.test(formula = LINCS ~ cell_id, data = tempDT, alternative = 'two.sided')
cat(
	sprintf('Wilcox: CMAP ~ Cell: %.2E\n', wilcoxCMAP$p.value),
	sprintf('Wilcox: LINCS ~ Cell: %.2E\n', wilcoxLINCS$p.value),
	
	sprintf('Median Rank CMAP - MCF7: %s\n', tempDT[cell_id == 'MCF7', median(CMAP)]),
	sprintf('Median Rank CMAP - PC3: %s\n', tempDT[cell_id == 'PC3', median(CMAP)]),
	
	sprintf('Median Rank LINCS - MCF7: %s\n', tempDT[cell_id == 'MCF7', median(LINCS)]),
	sprintf('Median Rank LINCS - PC3: %s\n', tempDT[cell_id == 'PC3', median(LINCS)]),
	
	# sprintf('Spearman: CMAP ~ Sum: %s\n', PairwiseSpearman(tempDT$CMAP, tempDT$sum)),
	# sprintf('Spearman: LINCS ~ Sum: %s\n', PairwiseSpearman(tempDT$LINCS, tempDT$sum)),
	
	# sprintf('Spearman: CMAP ~ Up Count: %s\n', PairwiseSpearman(tempDT$CMAP, tempDT$upCount)),
	# sprintf('Spearman: LINCS ~ Up Count: %s\n', PairwiseSpearman(tempDT$LINCS, tempDT$upCount)),
	
	# sprintf('Spearman: CMAP ~ Down Count: %s\n', PairwiseSpearman(tempDT$CMAP, tempDT$downCount)),
	# sprintf('Spearman: LINCS ~ Down Count: %s\n', PairwiseSpearman(tempDT$LINCS, tempDT$downCount)),
	
	# sprintf('Spearman: CMAP ~ Ratio: %s\n', PairwiseSpearman(tempDT$CMAP, tempDT$l10_ratio)),
	# sprintf('Spearman: LINCS ~ Ratio: %s\n', PairwiseSpearman(tempDT$LINCS, tempDT$l10_ratio)),
	
	sprintf('Spearman: CMAP ~ Sum [P-value]: %.2E\n', cor.test(tempDT$CMAP, tempDT$sum, alternative = 'two.sided', method = 'spearman')$p.value),
	
	file = outLogPath,
	sep = '',
	append = TRUE
)

tempGG = ggplot(tempDT, aes(x = sum, y = CMAP)) +
	geom_point(alpha = 0.3) +
	geom_smooth(method = 'lm', se = FALSE, colour = 'Red') +
	xlab('# DE Genes') +
	ylab('Retrieval Rank\n(CMap 1 Signature)')
tempPlotList = AddItemToList(tempPlotList, tempGG)

# 2. From 109 Conditions
tempDT = newdose.resultDT[, .(
	trt_iname,
	cell_id,
	# upCount = fixed.upCount,
	# downCount = fixed.downCount,
	# sum = fixed.sum,
	# l10_ratio = log10(fixed.ratio),
	CMAP = HYBRID.RANK,
	LINCS = REFERENCE.RANK
)]

tempGG = ggplot(tempDT, aes(x = sum, y = CMAP)) +
	geom_point(alpha = 0.3) +
	geom_smooth(method = 'lm', se = FALSE, colour = 'Red') +
	xlab('# DE Genes') +
	ylab('Retrieval Rank\n(CMap 1 Signature)')
tempPlotList = AddItemToList(tempPlotList, tempGG)

# Cowplot Output
tempCPlot = plot_grid(
	plotlist = tempPlotList,
	ncol = 2,
	nrow = 1,
	labels = LETTERS[1:2]
)
tempPath = sprintf('%sFig_S3.pdf', outDir)
save_plot(tempPath, plot = tempCPlot, ncol = 2, nrow = 1)