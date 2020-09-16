#' Quiet garbage collection
#'
#' @return Invisible null
#' @export
#'
#' @examples QuietGC()
QuietGC = function() {
	# Quiet Garbage Collection
	invisible(gc())
}
