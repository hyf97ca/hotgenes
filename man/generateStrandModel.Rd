% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generateStrandModel.R
\name{generateStrandModel}
\alias{generateStrandModel}
\title{Create a Strand Model

generateStrandModel Creates a model of a DNA strand in the form of a matrix (data frames are too slow) for further processing/plotting.}
\usage{
generateStrandModel(startBase, endBase, fc, chr, strand, scaling = 1000,
  updateProgressBar = NULL)
}
\arguments{
\item{startBase}{The starting base pair number}

\item{endBase}{The ending base pair number}

\item{fc}{The featurecounts object from rsubread after counting your alignment}

\item{chr}{Chromosome; Should match Chr in fc}

\item{strand}{+ or - strand to generate models for}

\item{scaling}{scaling factor in how many base pairs should be represented per index (ex. 1000) to
reduce computational complexity when simulating physics/plotting}

\item{updateProgressBar}{a function for linking with shiny's progress bar. This will be called if not NULL
when progress bar is to be updated (every 1000 processed features)}
}
\value{
A matrix containing the strand model.
}
\description{
Create a Strand Model

generateStrandModel Creates a model of a DNA strand in the form of a matrix (data frames are too slow) for further processing/plotting.
}
\examples{
\dontrun{
   #featureCounts() is from RSubread package; only run this if that is installed
   fc <- featureCounts(files=bamFile, annot.ext=annotation, useMetaFeatures=FALSE)
   #replace fc with the one loaded from featureCounts if you are using your own data
   sm <- generateStrandModel(startBase=1, endBase=195471971,
    fc=musCh1fc, chr="chr1", strand="-", scaling=100000)
   sm1 <- generateStrandModel(1, 195471971, musCh1fc, "chr1", "+")
   sm2 <- generateStrandModel(3000000, 3217000, musCh1fc, "chr1", "-")
}

}
