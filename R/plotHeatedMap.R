# plotHeatedMap.R

#' \code{plotHeatedMap}
#'
#' Creates a plot of a heated map of the specified strand model and location model. Plot is adjusted to roughly look like a DNA strand/track.
#'
#' @name plotHeatedMap
#' @param expression
#' Input matrix for adjusted RNA-seq values. the dataframe should be organized such that RNA-seq values run horizontally.\cr
#' This means that after loading column vectors the matrix should be transposed with t(). generateStrandModel already does this. In particular if you get "figure margins too large"
#' you are probably not transposing either this or the other parameter correctly.
#' @param locations
#' Input sequence/matrix for the base pair number cooresponding with expression. This should have the same number of columns (and rows if matrix) as expression.
#' @param palette
#' R color palette string, see hcl.colors. it is recommended to select a palette from \code{hcl.pals(type="sequential")}
#'
#' @return No return value, but tries to save par() and scipen(), create a plot with custom par and scipen() settings and restore old settings.
#'
#' @import graphics
#' @importFrom grDevices hcl.colors
#' @export
#'
#' @examples
#' \dontrun{
#'   sm <- generateStrandModel(1, 195471971, musCh1fc, "chr1", "-", 100000)
#'   x <- generateLocationModel(1, 195471971, 100000)
#'   sm <- simulateHeatSpread(sm, 0.001, 1000)
#'   plotHeatedMap(sm, x)
#'   plotHeatedMap(expression=sm, locations=x, palette="Heat")
#' }
#'

plotHeatedMap <- function(expression, locations, palette="Plasma")#c(2,0.5,0.4,2)
{
  #only returns parameters that can be set by subsequent par()
  old.par <- par(no.readonly = TRUE)
  old.scipen <- getOption("scipen")
  options(scipen=10)
  par(mfrow=c(max(12, nrow(expression)),1), mar=c(2,1,0.4,2))
  for (columnIndex in 1 : (nrow(expression)))
  {
    Z = as.matrix(expression[columnIndex,])
    if (is.matrix(locations) | is.data.frame(locations)){
      X = as.matrix(locations[columnIndex,])
    }
    else
    {
      X = as.matrix(locations)
    }
    head(Z)
    head(X)
    graphics::image(x=X, y=1, z = log(Z), yaxt = "n",
                    xlab = "bp", ylab = "", col = grDevices::hcl.colors(1024, palette=palette))
  }
  #restore old parameters. pin is supposed to be settable by par
  par(old.par)
  options(scipen=old.scipen)
}
