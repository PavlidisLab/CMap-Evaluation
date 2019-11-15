# Generate Figure S8

# Load Libraries
require(NCustomLibs)
require(data.table)
require(ggplot2)
require(cowplot)

# Declaring Global Variables
resultDir = 'Data/Analysis_CMAP_LINCS/Results/IntraCmp/DE/'
outDir = 'Data/Analysis_CMAP_LINCS/Output/'
outLogPath = sprintf('%sFig_S8.LOG', outDir)

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
tempPath = sprintf('%sintra.de.cmap.cordt.RDS.XZ', resultDir)
cmapResult = readRDS(tempPath)

tempPath = sprintf('%sintra.de.lincs.l3.cordt.RDS.XZ', resultDir)
lincs_3Result = readRDS(tempPath)

tempPath = sprintf('%sintra.de.lincs.l5.cordt.RDS.XZ', resultDir)
lincs_5Result = readRDS(tempPath)

# ----- Generate Plots -----
tempPlotList = list()

# Use Max Value as Summary
cmapResultFinal = cmapResult[, .(
	SPEAR = max(PWSPEAR.LM),
	JAC = max(PWJAC.LM)
), by = REF_CVAL]
lincs_3ResultFinal = lincs_3Result[, .(
	SPEAR = max(PWSPEAR.LM),
	JAC = max(PWJAC.LM)
), by = REF_CVAL]
lincs_5ResultFinal = lincs_5Result[, .(
	SPEAR = max(PWSPEAR.LM),
	JAC = max(PWJAC.LM)
), by = REF_CVAL]

# 1. Spearman Density Plot -----
tempDF = data.frame(
	MAX.SP = c(
		cmapResultFinal[, SPEAR],
		lincs_3ResultFinal[, SPEAR],
		lincs_5ResultFinal[, SPEAR]
	),
	GROUP = c(
		rep('CMap 1', nrow(cmapResultFinal)),
		rep('CMap 2-FC', nrow(lincs_3ResultFinal)),
		rep('CMap 2-MZS', nrow(lincs_5ResultFinal))
	)
)

tempGG = ggplot(tempDF, aes(x = MAX.SP, group = GROUP, color = GROUP)) +
	stat_density(geom = 'line', position = 'identity', show.legend = FALSE) +
	xlim(c(-0.7, 1.0)) +
	xlab('Within-Dataset DE Agreement\n(Spearman)') +
	ylab('Kernel Density') +
	scale_colour_manual(
		values = c('CMap 1' = 'Red', 'CMap 2-FC' = 'Blue', 'CMap 2-MZS' = 'Dark Green'),
		guide = guide_legend(title = 'Dataset')
	)
tempPlotList = AddItemToList(tempPlotList, tempGG)

# Summary Report
cat('----- Spearman Correlation -----\n', file = outLogPath, sep = '', append = TRUE)
for (resultIndex in c('cmapResultFinal', 'lincs_3ResultFinal', 'lincs_5ResultFinal')) {
	tempDT = .GlobalEnv[[resultIndex]]
	cat(
		sprintf('TYPE-Histogram: %s; COUNTS: %s\n', resultIndex, nrow(tempDT)),
		paste0(sprintf('%.3f', fivenum(tempDT[, SPEAR])), collapse = '; '),
		sprintf('\nMEAN: %.3f\n', mean(tempDT[, SPEAR], na.rm = TRUE)),
		file = outLogPath,
		sep = '',
		append = TRUE
	)
}

# 2. Jaccard Density Plot -----
tempDF = data.frame(
	MAX.JAC = c(
		cmapResultFinal[, JAC],
		lincs_3ResultFinal[, JAC],
		lincs_5ResultFinal[, JAC]
	),
	GROUP = c(
		rep('CMap 1', nrow(cmapResultFinal)),
		rep('CMap 2-FC', nrow(lincs_3ResultFinal)),
		rep('CMap 2-MZS', nrow(lincs_5ResultFinal))
	)
)

tempGG = ggplot(tempDF, aes(x = MAX.JAC, group = GROUP, color = GROUP)) +
	stat_density(geom = 'line', position = 'identity') +
	xlim(c(0, 0.8)) +
	xlab('Within-Dataset DE Agreement\n(Jaccard)') +
	ylab('Kernel Density') +
	scale_colour_manual(
		values = c('CMap 1' = 'Red', 'CMap 2-FC' = 'Blue', 'CMap 2-MZS' = 'Dark Green'),
		guide = guide_legend(title = 'Dataset')
	)
commonLegend = get_legend(tempGG)
tempPlotList = AddItemToList(tempPlotList, tempGG + theme(legend.position = 'none'))

# Summary Report
cat('\n----- Jaccard Index -----\n', file = outLogPath, sep = '', append = TRUE)
for (resultIndex in c('cmapResultFinal', 'lincs_3ResultFinal', 'lincs_5ResultFinal')) {
	tempDT = .GlobalEnv[[resultIndex]]
	cat(
		sprintf('TYPE-Histogram: %s; COUNTS: %s\n', resultIndex, nrow(tempDT)),
		paste0(sprintf('%.3f', fivenum(tempDT[, JAC])), collapse = '; '),
		sprintf('\nMEAN: %.3f\n', mean(tempDT[, JAC], na.rm = TRUE)),
		file = outLogPath,
		sep = '',
		append = TRUE
	)
}

# Cowplot Output -----
tempCPlot = plot_grid(
	plotlist = tempPlotList,
	labels = LETTERS[1:2],
	nrow = 1,
	ncol = 2
)
tempCPlot = plot_grid(
	tempCPlot,
	commonLegend,
	rel_widths = c(1, .25)
)
tempPath = sprintf('%sFig_S8.pdf', outDir)
save_plot(tempPath, plot = tempCPlot, ncol = 2.2, nrow = 1)