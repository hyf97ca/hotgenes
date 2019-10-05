# plotHeatedMap.R

#' \code{plotHeatedMap}
#'
#' \code{plotHeatedMap} is used to create a plot of a heated map of the specified matrix.
#'
#' @param expression
#' Input matrix for adjusted RNA-seq values. the dataframe should be organized such that RNA-seq values run horizontally;
#'  this means that after loading column vectors the matrix should be transposed with t(). generateStrandModel already does this. In particular if you get "figure margins too large"
#'  you are probably not transposing either this or the other parameter correctly.
#' @param locations
#' Input sequence/matrix for the base pair number cooresponding with expression. This should have the same number of columns (and rows if matrix) as expression.
#' @param palette
#' R color palette string, see hcl.colors. it is recommended to select a palette from \code{hcl.pals(type="sequential")}
#'
#' @export
#'
# @examples
# \dontrun
# {
#
# }
#

plotHeatedMap <- function(expression, locations, palette="Plasma")#c(2,0.5,0.4,2)
{
  old.par <- par(no.readonly = TRUE)
  old.scipen <- getOption("scipen")
  options(scipen=10)
  par(mfrow=c(max(12, nrow(expression)),1), mar=c(2,1,0.4,2))
  for (columnIndex in 1:(nrow(expression)))
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
    image(x=X, y=1, z = log(Z), yaxt = "n", xlab = "bp", ylab = "", col = hcl.colors(1024, palette=palette))
  }
  par(old.par)
  options(scipen=old.scipen)
}
