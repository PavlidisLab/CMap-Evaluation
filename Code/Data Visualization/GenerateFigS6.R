# Generate Figure S6

# Load Libraries
require(NCustomLibs)
require(ggplot2)
require(cowplot)

# Declaring Global Variables
resultDir = 'Data/Analysis_CMAP_LINCS/Results/'
outDir = 'Data/Analysis_CMAP_LINCS/Output/'
outLogPath = sprintf('%sFig_S6.LOG', outDir)

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
tempPath = sprintf('%smulti.pairwise.mainvar.DT.RDS.XZ', resultDir)
resultsDT = readRDS(tempPath)

# ----- Generate Plots -----

# Spearman Histograms -----
tempPlotList = list()

plotTypeVector = c(
	'FC-ALL' = 'L3.XD.JC.ALL',
	'MZS-ALL' = 'L5.XD.JC.ALL',
	'FC-LM' = 'L3.XD.JC.LM',
	'MZS-LM' = 'L5.XD.JC.LM'
)

for (plotIndex in 1:4) {
	tempDF = data.frame(x = resultsDT[[plotTypeVector[plotIndex]]])
	tempBreaks = seq(from = 0, to = 0.32, by = 0.01)
	tempGG = ggplot(tempDF, aes(x = x)) +
		geom_histogram(closed = 'right', breaks = tempBreaks, colour = 'Black') +
		xlab('Jaccard Index') +
		ylab('# Conditions') +
		xlim(c(0, 0.32))
	
	tempPlotList = AddItemToList(tempPlotList, tempGG)
	
	# Summary Report
	cat(
		sprintf('TYPE: %s\n', names(plotTypeVector)[plotIndex]),
		paste0(sprintf('%.3f', fivenum(tempDF$x)), collapse = ';  '),
		sprintf('\nMEAN: %.3f\n', mean(tempDF$x, na.rm = TRUE)),
		'\n\n',
		file = outLogPath,
		sep = '',
		append = TRUE
	)
}

# Cowplot Output -----
tempCPlot = plot_grid(
	plotlist = tempPlotList,
	nrow = 2,
	ncol = 2,
	labels = LETTERS[1:4]
)
tempPath = sprintf('%sFig_S6.pdf', outDir)
save_plot(tempPath, plot = tempCPlot, ncol = 2, nrow = 2)