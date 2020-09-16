#' Parallel XZ-compression (external XZ) implementation of saveRDS
#'
#' @param obj R object to save
#' @param file path to RDS file
#' @param threads number of parallel threads
#' @param compression level of compression
#'
#' @return NULL
#' @export
#'
#' @examples XZSaveRDS(obj, 'PATH')
XZSaveRDS = function(obj, file, threads = 32, compression = 6) {
	# Parallel XZ-Compression SaveRDS
	stopifnot(is.character(file))

	QuietFileRemove(paths = file, check = TRUE)

	xzCommand = sprintf('xz -z -T %s -%s > %s', threads, compression, file)
	xzConnection = pipe(description = xzCommand, open = 'wb')
	saveRDS(object = obj, file = xzConnection)
	close(xzConnection)
}
