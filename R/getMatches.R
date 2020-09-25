#' Searches a given set of genomic ranges for the motifs returned by mergeMotifs.
#' The background frequency is obtained from the total nucleotide frequencies in peaks.
#'
#' @param peaks A GenomicRanges or other subject input for matchMotifs.
#' @param genome A BSgenome object or other genome input for matchMotifs.
#' @param out Type of output to return.
#' @param motifs A PWMatrixList.
#' @param ... Additional arguments to \code{matchMotifs()}
#' @return A \code{matches} object.
#' @seealso \code{\link{mergeMotifs}}, \code{\link{motifmatchr::matchMotifs}}
#' @export
getMatches <- function(peaks,genome,out="matches",motifs=CrobustaMotifs,...){
	# get NT background freq from accessome
	bg <- letterFrequency(Views(genome,peaks),c("A","C","G","T"))
	bg <- apply(bg,2,sum)
	bg <- bg/sum(bg)

	#         motifs <- mergeMotifs()

	# match motifs to peaks
	matches <- matchMotifs(motifs,peaks,genome,bg=bg,out=out,...)
	row.names(matches) <- rowData(matches)$name
	rowData(matches)$width <- width(peaks)
	return(rmdup(motifs,matches,out))
}

#' Accepts a PWMatrixList and an output from motifmatchr::matchMotifs, and removes rows from the output corresponding to duplicate motif IDs.
#' Only the row with the greatest number of matches is left for each motif ID.
#' The background frequency is obtained from the total nucleotide frequencies in peaks.
#'
#' @param motifs The PWMatrixList used to compute matches.
#' @param matches An output from motifmatchr::matchMotifs.
#' @param out Type of output to return.
#' @return A \code{matches} object.
#' @seealso \code{\link{getMatches}}, \code{\link{motifmatchr::matchMotifs}}
#' @export
rmdup <- function(motifs,matches,metric='matches'){
	# split motif matches by TF
	sel <- split(names(motifs),ID(motifs))

	# count matches for each TF
	if(metric=="matches"){
		ct <- motifMatches(matches)
	}else if(metric=="scores"){
		ct <- motifScores(matches)
	}else return(matches)

	dupct <- lapply(sel,function(x) apply(ct[,x,drop=F],2,sum))
	# find motif for each TF with max number of matches
	sel <- sapply(dupct,function(x) names(x)[which.max(x)])

	bestmatch <- matches[,sel]
	return(bestmatch)
}

#' Accepts a vector of peak IDs, and performs a Poisson test for enrichmet of motifs in \code{peaks} versus a background of all elements in \code{matches}.
#'
#' @param peaks A character vector of peak IDs corresponding to the \code{row.names} of \code{matches}.
#' @param matches An output from motifmatchr::matchMotifs.
#' @param padj.method The method to be used by \code{p.adjust()}.
#' @param allow.duplicates Whether to remove duplicates from \code{peaks}.
#' @return A \code{data.frame} containing the test results.
#' @seealso \code{\link{getMatches}}, \code{\link{motifmatchr::matchMotifs}}
#' @export
matchPois <- function(peaks,matches,padj.method='fdr',allow.duplicates=F){
	if(!allow.duplicates) peaks <- unique(peaks)
	ct <- motifCounts(matches)
	bgct <- apply(ct,2,sum)
	testct <- apply(ct[peaks,],2,sum)
	mpkb <- bgct/sum(rowData(matches)$width)*1000
	testkbp <- sum(rowData(matches)[peaks,'width'])/1000
	res <- mapply(poisson.test,testct,testkbp,mpkb,SIMPLIFY=F)

	p <- sapply(res,'[[','p.value')
	padj <- p.adjust(p,method=padj.method)

	lor <- log2(testct/testkbp/mpkb)
	return(cbind(
		bgCount=bgct,
		testCount=testct,
		bgMotifPerKbp=mpkb,
		testMotifPerKbp=testct/testkbp,
		log2OR=lor,
		p=p,
		padj=padj
	))
}

#' Determines a matrix of potential regulatory TFs.
#'
#' @param geneToPeak A 2-column matrix with the first row corresponding to the \code{row.names} of \code{matches} and the second giving the Gene ID.
#' @param matches An output from \code{matchMotifs}.
#' @return A logical matrix with columns corresponding to the motifs in \code{matches} and rows corresponding to the genes in \code{geneToPeak}.
#' @seealso \code{\link{getOverlaps}}, \code{\link{matchMotifs}}
#' @export
regMat <- function(geneToPeak,matches){
	mat <- as.matrix(motifMatches(matches))
	geneID <- split(geneToPeak[,1],geneToPeak[,2])
	genemat <- lapply(geneID,function(x) mat[x,,drop=F])
	res <- t(sapply(genemat, function(x) as.numeric(apply(x,2,any))))
	colnames(res) <- colnames(matches)
	return(res)
}

