#' Obtain server-process combined ID
#'
#' @return character vector of length 1
#' @export
#'
#' @examples SysGetSPID()
SysGetSPID = function() {
	# Obtain Server-Process Combined ID
	NODE = unlist(strsplit(Sys.info()['nodename'], '\\.'))[1]
	PID = Sys.getpid()
	tempString = sprintf('%s-%s', NODE, PID)
	return(tempString)
}
