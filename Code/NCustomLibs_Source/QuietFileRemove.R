#' Quiet File Removal
#'
#' @param paths files/directories to remove
#' @param check enable file/directory existence check?
#'
#' @return Invisible logical (from file.remove)
#' @export
#'
#' @examples QuietFileRemove('REMOVE.THIS')
QuietFileRemove = function(paths, check = FALSE) {
	# Quiet File Removal, With Optional Presence Check
	stopifnot(is.character(paths))
	stopifnot(is.logical(check))

	if (check) {
		existsVector = file.exists(paths) | dir.exists(paths)
		invisible(suppressWarnings(file.remove(paths[existsVector])))
	} else {
		invisible(suppressWarnings(file.remove(paths)))
	}
}
