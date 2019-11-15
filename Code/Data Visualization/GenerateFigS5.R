# Generate Figure S5

# Load Libraries
require(NCustomLibs)
require(data.table)
require(stringi)
require(ggplot2)
require(cowplot)

# Declaring Global Variables
subsetPath = 'Data/Analysis_CMAP_LINCS/Subset/newdose.subset.RDAT.XZ'
resultDir = 'Data/Analysis_CMAP_LINCS/Results/'
infoPath = 'Data/LINCS/GSE92742-MISC/GSE92742_Broad_LINCS_pert_info.txt'
outDir = 'Data/Analysis_CMAP_LINCS/Output/'
outLogPath = sprintf('%sFig_S5.LOG', outDir)

options(
	stringsAsFactors = FALSE,
	warn = 1
)

# Global Preparation ------------------------------------------------------
# Creation of Output Directory
dir.create(outDir, recursive = TRUE, showWarnings = FALSE)

# Main Subroutine -------------------------------------------------------
# Preparation
cat('Identifying Touchstone Perturbagens.\n')
tempEnvir = new.env()
load(subsetPath, envir = tempEnvir)
attach(tempEnvir)

lincsMetadata = lincs_l3Metadata[match(cmapMetadata$new_cval, lincs_l3Metadata$new_cval), ]
detach(tempEnvir)
rm(tempEnvir)

pertInfo = fread(infoPath, header = TRUE, sep = '\t')
touchstoneList = pertInfo[(pert_type == 'trt_cp' & is_touchstone == '1'), pert_iname]
touchstoneList = stri_trans_toupper(touchstoneList)
rm(pertInfo)

touchstoneFilter = (lincsMetadata$trt_iname %in% touchstoneList) & (lincsMetadata$source_dataset == 'GSE92742')

cat('Loading Pre-Calculated Results.\n')
tempPath = sprintf('%sspearman.dist.all.newdose.RDS.XZ', resultDir)
allResults = readRDS(tempPath)

tempPath = sprintf('%sspearman.dist.landmark.newdose.RDS.XZ', resultDir)
landmarkResults = readRDS(tempPath)

cat(sprintf('Number of Entries: %s\n', sum(touchstoneFilter)), file = outLogPath, sep = '', append = TRUE)

# ----- Generate Plots -----

# Spearman Histograms -----
tempPlotList = list()

for (resultType in c('allResults', 'landmarkResults')) {
	for (dataType in c('CLFC', 'CLMZS')) {
		tempDF = data.frame(x = .GlobalEnv[[resultType]][[dataType]][touchstoneFilter])
		tempBreaks = seq(from = -0.25, to = 0.75, by = 0.05)
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
tempPath = sprintf('%sFig_S5.pdf', outDir)
save_plot(tempPath, plot = tempCPlot, ncol = 2, nrow = 2)