#' Generate Temporary RDS (TRDS) path
#'
#' @param trdsDir temporary output directory (usually '/dev/shm/')
#' @param spid current server-process id
#' @param type operation type (e.g. correlation function)
#' @param idx batch index
#'
#' @return character vector of TRDS file path(s)
#' @export
#'
#' @examples .TRDSPathGen('/dev/shm/', 'SPID', 'NULL', 1)
.TRDSPathGen = function(trdsDir, spid, type, idx) {
	# TRDS Path Generator
	tempPrefix = sprintf('%s%s.%s.', trdsDir, spid, type)
	return(paste0(tempPrefix, idx, '.trds'))
}
