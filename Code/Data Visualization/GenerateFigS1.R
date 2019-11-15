# Generate Figure S1

# Load Libraries
require(NCustomLibs)
require(data.table)
require(ggplot2)
require(cowplot)

# Declaring Global Variables
resultDir = 'Data/Analysis_CMAP_LINCS/Results/'
outDir = 'Data/Analysis_CMAP_LINCS/Output/'

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
resultDT = readRDS(tempPath)

# ----- Generate Plots -----

# Histogram + ECDF -----
tempPlotList = list()
tempDT = data.table(
	GROUP = c(
		rep('CMap 1', nrow(resultDT)),
		rep('CMap 2', nrow(resultDT))
	),
	VAL = c(
		resultDT[, HYBRID.RANK],
		resultDT[, REFERENCE.RANK]
	)
)

# 1. Frequency Plot
tempGG = ggplot(tempDT, aes(x = VAL, group = GROUP, colour = GROUP)) +
	geom_freqpoly(closed = 'right', show.legend = FALSE, binwidth = 20) +
	xlab('Compound Retrieval Rank') +
	ylab('# Signatures') +
	scale_colour_manual(
		values = c('CMap 1' = 'Red', 'CMap 2' = 'Blue'),
		guide = guide_legend(title = 'Dataset')
	) +
	geom_vline(xintercept = 68, colour = 'Grey', linetype = 2)
tempPlotList = AddItemToList(tempPlotList, tempGG)

# 2. ECDF Plot
tempGG = ggplot(tempDT, aes(x = VAL, group = GROUP, colour = GROUP)) +
	geom_step(stat = 'ecdf') +
	xlab('Compound Retrieval Rank') +
	ylab('ECDF') +
	scale_colour_manual(
		values = c('CMap 1' = 'Red', 'CMap 2' = 'Blue'),
		guide = guide_legend(title = 'Dataset'),
		labels = c('CMap 1' = 'CMap 1', 'CMap 2' = 'CMap 2 (Control)')
	) +
	geom_vline(xintercept = 68, colour = 'Grey', linetype = 2)
commonLegend = get_legend(tempGG)
tempPlotList = AddItemToList(tempPlotList, tempGG + theme(legend.position = 'none'))

# Cowplot Output
tempCPlot = plot_grid(
	plotlist = tempPlotList,
	nrow = 1,
	ncol = 2,
	labels = LETTERS[1:2]
)
tempCPlot = plot_grid(
	tempCPlot,
	commonLegend,
	rel_widths = c(1, .25)
)
tempPath = sprintf('%sFig_S1.pdf', outDir)
save_plot(tempPath, plot = tempCPlot, ncol = 2, nrow = 1)