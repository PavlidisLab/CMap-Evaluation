# Generate Figure 6

# Load Libraries
require(NCustomLibs)
require(data.table)
require(ggplot2)
require(cowplot)

# Declaring Global Variables
resultDir = 'Data/Analysis_CMAP_LINCS/Results/IntraCmp/EXP/Packed/'
outDir = 'Data/Analysis_CMAP_LINCS/Output/'
outLogPath = sprintf('%sFig_6.LOG', outDir)

options(
	stringsAsFactors = FALSE,
	warn = 1
)

# Global Preparation ------------------------------------------------------
# Creation of Output Directory
dir.create(outDir, recursive = TRUE, showWarnings = FALSE)

# Main Subroutine -------------------------------------------------------
# Preparation
tempPath = sprintf('%scomplete.intra.dataset.lincs.l5.dt.RDS.XZ', resultDir)
lincs5DT = as.data.table(readRDS(tempPath))
setorder(lincs5DT, -DE.PWS.ALL)
final.lincs5DT = lincs5DT[, .SD[1], by = REF_CVAL]

tempPath = sprintf('%scomplete.intra.dataset.lincs.l3.treatment.dt.RDS.XZ', resultDir)
lincs3DT_Treat = as.data.table(readRDS(tempPath))
setorder(lincs3DT_Treat, -DE.PWS.ALL)
final.lincs3DT_Treat = lincs3DT_Treat[, .SD[1], by = REF_CVAL]

tempPath = sprintf('%scomplete.intra.dataset.lincs.l3.control.dt.RDS.XZ', resultDir)
lincs3DT_Control = as.data.table(readRDS(tempPath))
setorder(lincs3DT_Control, -DE.PWS.ALL)
final.lincs3DT_Control = lincs3DT_Control[, .SD[1], by = REF_CVAL]

# ----- Generate Plots -----

# Kernel Density Plot -----
final.tempDT = data.table(
	MAX.SP = c(
		final.lincs5DT$WTH.MINMAX.ALL,
		final.lincs3DT_Treat$WTH.MINMAX.ALL,
		final.lincs3DT_Control$WTH.MINMAX.ALL
	),
	GROUP = c(
		rep('MZS-Treatment', nrow(final.lincs5DT)),
		rep('FC-Treatment', nrow(final.lincs3DT_Treat)),
		rep('FC-Control', nrow(final.lincs3DT_Control))
	)
)
tempGG = ggplot(final.tempDT, aes(x = MAX.SP, group = GROUP, color = GROUP)) +
	stat_density(geom = 'line', position = 'identity') +
	geom_vline(xintercept = 0.98, colour = 'Grey', linetype = 2) +
	xlim(c(0.3, 1.0)) +
	xlab('"In-Replicate" Correlation') +
	ylab('Kernel Density') +
	scale_colour_manual(
		values = c('MZS-Treatment' = 'Red', 'FC-Treatment' = 'Blue', 'FC-Control' = 'Dark Green'),
		guide = guide_legend(title = 'CMap 2 Samples')
	)

# Summary Report
for (resultIndex in c('MZS-Treatment', 'FC-Treatment', 'FC-Control')) {
	tempDT = final.tempDT[GROUP == resultIndex]
	tempDensity = density(tempDT$MAX.SP)
	cat(
		sprintf('TYPE-Histogram: %s; COUNTS: %s\n', resultIndex, nrow(tempDT)),
		paste0(sprintf('%.3f', fivenum(tempDT$MAX.SP)), collapse = '; '),
		sprintf('\nMean: %.3f\n', mean(tempDT$MAX.SP)),
		sprintf('\nMode: %.3f\n', tempDensity$x[which.max(tempDensity$y)]),
		file = outLogPath,
		sep = '',
		append = TRUE
	)
}

# Cowplot Output -----
tempPath = sprintf('%sFig_6.pdf', outDir)
save_plot(tempPath, plot = tempGG, base_width = 7)