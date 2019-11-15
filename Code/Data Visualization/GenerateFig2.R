# Generate Figure 2

# Load Libraries
require(NCustomLibs)
require(ggplot2)
require(cowplot)

# Declaring Global Variables
resultDir = 'Data/Analysis_CMAP_LINCS/Results/'
outDir = 'Data/Analysis_CMAP_LINCS/Output/'
outLogPath = sprintf('%sFig_2.LOG', outDir)

options(
	stringsAsFactors = FALSE,
	warn = 1
)

# Global Preparation ------------------------------------------------------
# Creation of Output Directory
dir.create(outDir, recursive = TRUE, showWarnings = FALSE)

# Main Subroutine -------------------------------------------------------
# Preparation
cat('Loading Pre-Calculated Results.\n')
tempPath = sprintf('%sspearman.dist.all.newdose.RDS.XZ', resultDir)
allResults = readRDS(tempPath)

tempPath = sprintf('%sspearman.dist.landmark.newdose.RDS.XZ', resultDir)
landmarkResults = readRDS(tempPath)

# ----- Generate Plots -----

# Spearman Histograms -----
tempPlotList = list()

for (resultType in c('allResults', 'landmarkResults')) {
	for (dataType in c('CLFC', 'CLMZS')) {
		tempDF = data.frame(x = .GlobalEnv[[resultType]][[dataType]])
		tempBreaks = hist(tempDF$x, breaks = 'FD')$breaks
		tempGG = ggplot(tempDF, aes(x = x)) +
			geom_histogram(closed = 'right', breaks = tempBreaks, colour = 'Black') +
			xlab('Spearman (r)') +
			ylab('# Conditions') +
			xlim(c(-0.25, 0.75))
		
		tempPlotList[[sprintf('%s-%s', resultType, dataType)]] = tempGG
		
		# Summary Report
		cat(
			sprintf('TYPE: (%s - %s)\n', resultType, dataType),
			paste0(sprintf('%.3f', fivenum(tempDF$x)), collapse = ';  '),
			sprintf('\nMEAN: %.3f\n', mean(tempDF$x, na.rm = TRUE)),
			'\n\n',
			file = outLogPath,
			sep = '',
			append = TRUE
		)
	}
}

# Cowplot Output -----
tempCPlot = plot_grid(
	plotlist = tempPlotList,
	nrow = 2,
	ncol = 2,
	labels = LETTERS[1:4]
)
tempPath = sprintf('%sFig_2.pdf', outDir)
save_plot(tempPath, plot = tempCPlot, ncol = 2, nrow = 2)