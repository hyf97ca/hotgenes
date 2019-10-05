#' \code{generateLocationModel}
#'
#' @param startBase The starting base pair number
#' @param endBase The ending base pair number
#' @param scaling scaling factor in how many base pairs should be represented per index (ex. 1000) to reduce computational complexity when simulating physics/plotting;
#' Should match value used in generateStrandModel
#'
#' @return A sequence containing the location model.
#' @export
#'
# @examples
#\dontrun
#{
#    generateLocationModel(1, 195471971)
#}
#

generateLocationModel <- function(startBase, endBase, scaling=1000)
{
  return(t(as.matrix(seq(from=startBase, to=endBase + scaling + 1, by=scaling))))
  #return(format(t(as.matrix(seq(from=startBase, to=endBase + scaling, by=scaling))), scientific = FALSE))
}
