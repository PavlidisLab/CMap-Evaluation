# Generate Figure S12

# Load Libraries
require(NCustomLibs)
require(ggplot2)
require(cowplot)

# Declaring Global Variables
resultPath = 'Data/Analysis_CMAP_LINCS/Results/IntraCmp/DE/intra.de.cmap.cordt.RDS.XZ'
outDir = 'Data/Analysis_CMAP_LINCS/Output/'

options(
	stringsAsFactors = FALSE,
	warn = 1
)

# Global Preparation ------------------------------------------------------
# Creation of Output Directory
dir.create(outDir, recursive = TRUE, showWarnings = FALSE)

# Main Subroutine -------------------------------------------------------
# Preparation
resultDT = readRDS(resultPath)

# ----- Generate Plots -----

# Spearman Density Plot -----
# Valproic Acid
tempDT = resultDT[TRT_INAME == 'VALPROIC-ACID' & TRT_CELL == 'MCF7', .(DOSE = log10(TRT_DOSE), COR = PWSPEAR.ALL)]
tempSumDT = tempDT[, .(MCOR = mean(COR)), by = DOSE]
tempSumDT[, `:=`(
	xStart = DOSE - 0.1,
	xEnd = DOSE + 0.1
)]
tempGG = ggplot(tempDT, aes(x = DOSE, y = COR)) +
	geom_point(alpha = 0.2) +
	geom_segment(
		mapping = aes(x = xStart, xend = xEnd, y = MCOR, yend = MCOR),
		data = tempSumDT, inherit.aes = FALSE, colour = 'Red'
	) +
	xlab(expression('log'[10]*' Concentration')) +
	ylab('All Pairwise Correlations')

# Cowplot Output -----
tempPath = sprintf('%sFig_S12.pdf', outDir)
save_plot(tempPath, plot = tempGG)