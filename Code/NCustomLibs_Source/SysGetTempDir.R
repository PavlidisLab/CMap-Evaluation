#' Generate temporary directory (preferably on RAM or SSD)
#'
#' @param useRAM use RAM disk?
#' @param useSSD use SSD?
#'
#' @return character value for temporary directory path
#' @export
#'
#' @examples x = SysGetTempDir(useRAM = FALSE, useSSD = TRUE)
SysGetTempDir = function(useRAM = TRUE, useSSD = TRUE) {
	# Generate Temporary Directory, preferably on RAM or SSD
	# Get Server Node
	currentNode = unlist(strsplit(Sys.info()['nodename'], '\\.'), use.names = FALSE)[1]

	# Check SSD
	if (useSSD) {
		cat('Using SSD Caching System.\n')
		tempDirPath = '/'
		return(tempDirPath)
	}

	# Check RAM
	if (useRAM) {
		cat('Using RAM Caching System.\n')
		tempDirPath = '/dev/shm/'
		return(tempDirPath)

	} else {
		# Revert to HDD
		cat('Using HDD Caching System.\n')
		tempDirPath = '/'
		return(tempDirPath)
	}
}
