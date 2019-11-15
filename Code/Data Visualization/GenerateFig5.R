# Generate Figure 5

# Load Libraries
require(NCustomLibs)
require(ggplot2)
require(cowplot)

# Declaring Global Variables
resultPath = 'Data/Analysis_CMAP_LINCS/Results/IntraCmp/EXP/Packed/complete.intra.dataset.lincs.l5.dt.RDS.XZ'
outDir = 'Data/Analysis_CMAP_LINCS/Output/'
outLogPath = sprintf('%sFig_5.LOG', outDir)

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

DoseReproducibilityFigure = function(chemicalList) {
	tempList = list()
	for (currentChemical in chemicalList) {
		tempDT = resultDT[TRT_INAME == currentChemical, .(DOSE = log10(TRT_DOSE), COR = DE.PWS.ALL, CELL = TRT_CELL)]
		tempGG = ggplot(tempDT, aes(x = DOSE, y = COR, group = CELL, colour = CELL)) +
			geom_point(alpha = 0.2) +
			geom_smooth(method = 'loess', se = FALSE, span = 1) +
			xlab(expression('log'[10]*' Concentration')) +
			ylab('All Pairwise Correlations') +
			scale_colour_manual(
				values = c('MCF7' = 'Red', 'PC3' = 'Blue'),
				guide = guide_legend(title = 'Cell Line')
			)
		tempList = AddItemToList(tempList, tempGG)
	}
	return(tempList)
}

# Main Subroutine -------------------------------------------------------
# Preparation
resultDT = readRDS(resultPath)

# ----- Generate Plots -----

# Spearman Density Plot -----
# ** All CMap 2-MZS Entries
chemicalList = c('VORINOSTAT', 'TRICHOSTATIN-A', 'GELDANAMYCIN', 'WORTMANNIN')

# Iterative Processing
tempPlotList = DoseReproducibilityFigure(chemicalList)

# Cowplot Output -----
tempCPlot = plot_grid(
	plotlist = tempPlotList,
	nrow = 2,
	ncol = 2,
	labels = LETTERS[1:4]
)
tempPath = sprintf('%sFig_5.pdf', outDir)
save_plot(tempPath, plot = tempCPlot, ncol = 2, nrow = 2)