# Generate Figure S2

# Load Libraries
require(NCustomLibs)
require(data.table)
require(ggplot2)
require(cowplot)

# Declaring Global Variables
resultDir = 'Data/Analysis_CMAP_LINCS/Results/'
outDir = 'Data/Analysis_CMAP_LINCS/Output/'
outLogPath = sprintf('%sFig_S2.LOG', outDir)

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
tempPath = sprintf('%sl1000.compounds.newdose.expanded.RDS.XZ', resultDir)
resultsDT = readRDS(tempPath)

# ----- Generate Plots -----
tempPlotList = list()

tempDT = data.table(
	RANKS = c(
		resultsDT[, REFERENCE.RANK],
		resultsDT[, HYBRID.RANK],
		resultsDT[, S300.RANK],
		resultsDT[, S20.RANK],
		resultsDT[, L100.RANK]
	),
	GROUP = c(
		rep('CMap 2 (Control)', nrow(resultsDT)),
		rep('CMap 1 [Hybrid]', nrow(resultsDT)),
		rep('CMap 1 [Top-300, BING]', nrow(resultsDT)),
		rep('CMap 1 [Top-20, BING]', nrow(resultsDT)),
		rep('CMap 1 [Top-100, LM]', nrow(resultsDT))
	)
)

ggColourScale = scale_color_manual(
	values = c(
		'CMap 2 (Control)' = 'Blue',
		'CMap 1 [Hybrid]' = 'Red',
		'CMap 1 [Top-300, BING]' = 'Orange',
		'CMap 1 [Top-20, BING]' = 'Purple 4',
		'CMap 1 [Top-100, LM]' = 'Green'
	), guide = guide_legend(title = 'Dataset [Threshold]')
)

# 1. Frequency Plot -----
tempGG = ggplot(tempDT, aes(x = RANKS, group = GROUP, colour = GROUP)) +
	geom_freqpoly(closed = 'right', show.legend = FALSE, binwidth = 20) +
	xlab('Compound Retrieval Rank') +
	ylab('# Signatures') +
	ggColourScale +
	geom_vline(xintercept = 68, colour = 'Grey', linetype = 2)
tempPlotList = AddItemToList(tempPlotList, tempGG)

# 2. ECDF Plot -----
tempGG = ggplot(tempDT, aes(x = RANKS, group = GROUP, colour = GROUP, )) +
	geom_step(stat = 'ecdf') +
	xlab('Compound Retrieval Rank') +
	ylab('ECDF') +
	ggColourScale +
	geom_vline(xintercept = 68, colour = 'Grey', linetype = 2)
commonLegend = get_legend(tempGG)
tempPlotList = AddItemToList(tempPlotList, tempGG + theme(legend.position = 'none'))

# Summary Report
for (resultIndex in c('REFERENCE.RANK', 'HYBRID.RANK', 'S300.RANK', 'S20.RANK', 'L100.RANK')) {
	tempDT = resultsDT[[resultIndex]]
	cat(
		sprintf('TYPE-Histogram: %s; COUNTS: %s\n', resultIndex, length(tempDT)),
		sprintf('MEDIAN: %.1f\n', median(tempDT, na.rm = TRUE)),
		file = outLogPath,
		sep = '',
		append = TRUE
	)
}

# Cowplot Output -----
tempCPlot = plot_grid(
	plotlist = tempPlotList,
	nrow = 1,
	ncol = 2,
	labels = LETTERS[1:2]
)
tempCPlot = plot_grid(
	tempCPlot,
	commonLegend,
	rel_widths = c(1, .3)
)
tempPath = sprintf('%sFig_S2.pdf', outDir)
save_plot(tempPath, plot = tempCPlot, ncol = 2.2, nrow = 1)