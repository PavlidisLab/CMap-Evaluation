#' XZ-supported (external XZ) implementation of data.table::fread
#'
#' @param file input file path
#' @param sep separator character
#' @param header has header?
#'
#' @return data.table object
#' @export
#'
#' @examples XZFread('FILE.xz')
XZFread = function(file, sep = '\t', header = TRUE) {
	# XZ-Supported data.table Fread
	stopifnot(is.character(file))
	require(data.table)

	# Pre-Decompress File
	tempOutDir = sprintf('%sXZFread_%s/', SysGetTempDir(useRAM = TRUE, useSSD = TRUE), SysGetSPID())
	dir.create(tempOutDir, showWarnings = FALSE, recursive = TRUE)

	secondStamp = format(Sys.time(), format = '%s')
	fileID = .FileNameScramble(file)
	tempFilePath = sprintf('%s%s_%s.TMP', tempOutDir, secondStamp, fileID)
	xzCommand = sprintf('xz -d -c %s > %s', file, tempFilePath)
	system(command = xzCommand, ignore.stdout = FALSE, ignore.stderr = TRUE)

	# Fast-Read
	tempDT = fread(file = tempFilePath, sep = sep, header = header, showProgress = FALSE)

	# Cleanup
	file.remove(tempFilePath)

	return(tempDT)
}

.FileNameScramble = function(charVec) {
	# Hidden Function For Scrambling File Names
	require(stringi)

	charVec = stri_trans_tolower(charVec)
	charVec = tail(unlist(strsplit(charVec, split = '/'), use.names = FALSE), n = 1)
	charVec = unlist(strsplit(charVec, split = ''), use.names = FALSE)
	charVec = paste0(na.omit(match(charVec, letters)), collapse = '')
	return(charVec)
}
