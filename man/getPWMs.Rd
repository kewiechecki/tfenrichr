% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getKmers.R
\name{getPWMs}
\alias{getPWMs}
\title{Returns a PWMatrixList of all permutations of some NT length k.}
\usage{
getPWMs(k, bases = c("A", "C", "G", "T"))
}
\arguments{
\item{k}{bp length of each mer.}
}
\value{
A PWMatrixList
}
\description{
Returns a PWMatrixList of all permutations of some NT length k.
}
\seealso{
\code{\link{PWMatrixList}}
}
