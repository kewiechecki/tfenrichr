% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getMatches.R
\name{rmdup}
\alias{rmdup}
\title{Accepts a PWMatrixList and an output from motifmatchr::matchMotifs, and removes rows from the output corresponding to duplicate motif IDs.
Only the row with the greatest number of matches is left for each motif ID.
The background frequency is obtained from the total nucleotide frequencies in peaks.}
\usage{
rmdup(motifs, matches, metric = "matches")
}
\arguments{
\item{motifs}{The PWMatrixList used to compute matches.}

\item{matches}{An output from motifmatchr::matchMotifs.}

\item{out}{Type of output to return.}
}
\value{
A \code{matches} object.
}
\description{
Accepts a PWMatrixList and an output from motifmatchr::matchMotifs, and removes rows from the output corresponding to duplicate motif IDs.
Only the row with the greatest number of matches is left for each motif ID.
The background frequency is obtained from the total nucleotide frequencies in peaks.
}
\seealso{
\code{\link{getMatches}}, \code{\link{motifmatchr::matchMotifs}}
}
