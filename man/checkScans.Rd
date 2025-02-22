% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/checkScans.R
\name{checkScans}
\alias{checkScans}
\title{Check if any samples are missing spectra}
\usage{
checkScans(spectra)
}
\arguments{
\item{spectra}{\code{Spectra} object}
}
\description{
Display error if they do.
Currently stops if a sample is detected that has less than 90\% of the mean
of spectra in all samples.
}
\author{
Pablo Vangeenderhuysen
}
