#' Basic pairwise Spearman correlation
#'
#' @param x numeric vector
#' @param y numeric vector
#'
#' @return Spearman correlation value
#' @export
#'
#' @examples result = PairwiseSpearman(X, Y)
PairwiseSpearman = function(x, y) {
	# Basic Pairwise Spearman
	require(ccaPP)
	stopifnot(is.numeric(x))
	stopifnot(is.numeric(y))

	ok = complete.cases(x, y)
	return(corSpearman(x[ok], y[ok]))
}
