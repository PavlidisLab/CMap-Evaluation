# Generate Figure S4

# Load Libraries
require(NCustomLibs)
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
tempPath = sprintf('%sl1000.compounds.expanded.RDS.XZ', resultDir)
resultDT.588 = readRDS(tempPath)

tempPath = sprintf('%sl1000.compounds.newdose.expanded.RDS.XZ', resultDir)
resultDT.97 = readRDS(tempPath)

# ----- Generate Plots -----
tempPlotList = list()

# 1. 588 Conditions Plot
tempGG = ggplot(resultDT.588, aes(x = final.lincsRank, y = final.cmapRank)) +
	geom_point(alpha = 0.2) +
	geom_vline(xintercept = 68, colour = 'Red', linetype = 2) +
	geom_hline(yintercept = 68, colour = 'Red', linetype = 2) +
	xlab('Self-Query Retrieval Rank\n(N = 588)') +
	ylab('Retrieval Rank\n(CMap 1 Signature)')
tempPlotList = AddItemToList(tempPlotList, tempGG)

# 2. 97 Conditions Plot
tempGG = ggplot(resultDT.97, aes(x = REFERENCE.RANK, y = HYBRID.RANK)) +
	geom_point(alpha = 0.2) +
	geom_vline(xintercept = 68, colour = 'Red', linetype = 2) +
	geom_hline(yintercept = 68, colour = 'Red', linetype = 2) +
	xlab('Self-Query Retrieval Rank\n(N = 97; "Hybrid Threshold")') +
	ylab('Retrieval Rank\n(CMap 1 Signature)')
tempPlotList = AddItemToList(tempPlotList, tempGG)

# Cowplot Output -----
tempCPlot = plot_grid(
	plotlist = tempPlotList,
	nrow = 1,
	ncol = 2,
	labels = LETTERS[1:2]
)
tempPath = sprintf('%sFig_S4.pdf', outDir)
save_plot(tempPath, plot = tempCPlot, ncol = 2, nrow = 1)