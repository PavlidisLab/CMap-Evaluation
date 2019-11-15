# Generate Linear Model (In Main Manuscript)

# Load Libraries
require(NCustomLibs)
require(data.table)
require(lm.beta)

# Declaring Global Variables
resultDir = 'Data/Analysis_CMAP_LINCS/Results/IntraCmp/EXP/Packed/'

options(
	stringsAsFactors = FALSE,
	warn = 1
)

# Main Subroutine -------------------------------------------------------
# Preparation
cat('Loading Pre-Calculated Results.\n')
tempPath = sprintf('%scomplete.intra.dataset.lincs.l5.dt.RDS.XZ', resultDir)
completeDT = readRDS(tempPath)

# Model Fitting
setorder(completeDT, -DE.PWS.ALL)
completeDT = completeDT[, .SD[1], by = REF_CVAL]
completeDT[, DE.COUNT.ALL := log10(DE.COUNT.ALL)]

lmObject = completeDT[, lm(DE.PWS.ALL ~ DE.COUNT.ALL + WTH.MINMAX.ALL + BTW.MAX.ALL)]
std.lmObject = lm.beta(lmObject)

# Report on Linear Model Object
cat('Reporting Results In-console.\n')
summary(lmObject)
print(std.lmObject)