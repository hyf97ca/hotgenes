% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotHeatedMap.R
\name{plotHeatedMap}
\alias{plotHeatedMap}
\title{\code{plotHeatedMap}}
\usage{
plotHeatedMap(expression, locations, palette = "Plasma")
}
\arguments{
\item{expression}{Input matrix for adjusted RNA-seq values. the dataframe should be organized such that RNA-seq values run horizontally.\cr
This means that after loading column vectors the matrix should be transposed with t(). generateStrandModel already does this. In particular if you get "figure margins too large"
you are probably not transposing either this or the other parameter correctly.}

\item{locations}{Input sequence/matrix for the base pair number cooresponding with expression. This should have the same number of columns (and rows if matrix) as expression.}

\item{palette}{R color palette string, see hcl.colors. it is recommended to select a palette from \code{hcl.pals(type="sequential")}}
}
\value{
No return value, but tries to save par() and scipen(), create a plot with custom par and scipen() settings and restore old settings.
}
\description{
Creates a plot of a heated map of the specified strand model and location model. Plot is adjusted to roughly look like a DNA strand/track.
}
\examples{
\dontrun{
  sm <- generateStrandModel(1, 195471971, musCh1fc, "chr1", "-", 100000)
  x <- generateLocationModel(1, 195471971, 100000)
  sm <- simulateHeatSpread(sm, 0.001, 1000)
  plotHeatedMap(sm, x)
  plotHeatedMap(expression=sm, locations=x, palette="Heat")
}

}
