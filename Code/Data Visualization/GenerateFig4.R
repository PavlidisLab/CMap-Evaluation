# Generate Figure 4

# Load Libraries
require(NCustomLibs)
require(data.table)
require(ggplot2)
require(cowplot)

# Declaring Global Variables
resultDir = 'Data/Analysis_CMAP_LINCS/Results/IntraCmp/DE/'
outDir = 'Data/Analysis_CMAP_LINCS/Output/'
outLogPath = sprintf('%sFig_4.LOG', outDir)

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
tempPath = sprintf('%sintra.de.cmap.cordt.RDS.XZ', resultDir)
cmapResult = readRDS(tempPath)

tempPath = sprintf('%sintra.de.lincs.l3.cordt.RDS.XZ', resultDir)
lincs_3Result = readRDS(tempPath)

tempPath = sprintf('%sintra.de.lincs.l5.cordt.RDS.XZ', resultDir)
lincs_5Result = readRDS(tempPath)

# ----- Generate Plots -----

# Use Max Value as Summary
cmapResult = cmapResult[, max(PWSPEAR.ALL), by = REF_CVAL]
lincs_3Result = lincs_3Result[, max(PWSPEAR.ALL), by = REF_CVAL]
lincs_5Result = lincs_5Result[, max(PWSPEAR.ALL), by = REF_CVAL]

# Spearman Density Plot -----
tempDF = data.frame(
	MAX.SP = c(
		cmapResult$V1,
		lincs_3Result$V1,
		lincs_5Result$V1
	),
	GROUP = c(
		rep('CMap 1', nrow(cmapResult)),
		rep('CMap 2-FC', nrow(lincs_3Result)),
		rep('CMap 2-MZS', nrow(lincs_5Result))
	)
)

tempGG = ggplot(tempDF, aes(x = MAX.SP, group = GROUP, color = GROUP)) +
	stat_density(geom = 'line', position = 'identity') +
	geom_vline(xintercept = 0.9, colour = 'Grey', linetype = 2) +
	xlim(c(-0.7, 1.0)) +
	xlab('Within-Dataset DE Agreement') +
	ylab('Kernel Density') +
	scale_colour_manual(
		values = c('CMap 1' = 'Red', 'CMap 2-FC' = 'Blue', 'CMap 2-MZS' = 'Dark Green'),
		guide = guide_legend(title = 'Dataset')
	)

# Summary Report
for (resultIndex in c('cmapResult', 'lincs_3Result', 'lincs_5Result')) {
	tempDF = .GlobalEnv[[resultIndex]]
	cat(
		sprintf('TYPE-Histogram: %s; COUNTS: %s\n', resultIndex, nrow(tempDF)),
		paste0(sprintf('%.3f', fivenum(tempDF$V1)), collapse = '; '),
		sprintf('\nMEAN: %.3f\n', mean(tempDF$V1, na.rm = TRUE)),
		file = outLogPath,
		sep = '',
		append = TRUE
	)
}

# Cowplot Output -----
tempPath = sprintf('%sFig_4.pdf', outDir)
save_plot(tempPath, plot = tempGG, base_width = 7)