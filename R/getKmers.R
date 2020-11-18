#' Returns position weight matrices of all permutations of some NT length k.
#'
#' @param k bp length of each mer.
#' @return A list of logical matrices.
#' @seealso \code{\link{PWMatrixList}}
#' @export
getMatrices <- function(k,bases=c("A","C","G","T")){
	pwms <- expand.grid(as.data.frame(replicate(k,bases)))
	mat <- sapply(as.data.frame(t(pwms)), function(x) {
		res <- sapply(x,function(y) {
			y==bases
		})*1
		row.names(res) <- bases
		return(res)
	},simplify=F)
	names(mat) <- apply(pwms,1,paste,collapse="")
	return(mat)
}

#' Returns a PWMatrixList of all permutations of some NT length k.
#'
#' @param k bp length of each mer.
#' @return A PWMatrixList
#' @seealso \code{\link{PWMatrixList}}
#' @export
getPWMs <- function(k,bases=c("A","C","G","T")){
	mat <- getMatrices(k,bases)
	do.call(TFBSTools::PWMatrixList,mapply(
	       TFBSTools::PWMatrix,
	       names(mat),
	       profileMatrix=mat
	))
}

