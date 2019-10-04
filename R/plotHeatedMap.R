# plotHeatedMap.R

#' \code{plotHeatedMap}
#'
#' \code{plotHeatedMap} is used to create a plot of a heated map of the specified data frame.
#'
#' @param expression
#' Input dataframe for adjusted RNA-seq values. the dataframe should be organized such that RNA-seq values run horizontally;
#'  this means that after loading column vectors the data frame should be transposed with t().
#' @param locations
#' Input dataframe for the base pair number cooresponding with expression. This should have the same dimensions and organization as expression.
#'
#' @export
#'
# @examples
# \dontrun
# {
#
# }
#

plotHeatedMap <- function(expression, locations)#c(2,0.5,0.4,2)
{
  par(mfrow=c(max(12, nrow(expression)),1), mar=c(2,0.5,0.4,2))
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
    image(x=X, y=1, z = log(Z), yaxt = "n", xlab = "bp", ylab = "", col = heat.colors(1024))
  }
}
