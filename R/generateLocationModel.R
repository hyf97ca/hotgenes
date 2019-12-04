#generateLocationModel.R

#' Create a 'Location Model'
#'
#' Convenience function for generating the locations on the DNA strand where exons are so that the x scale on the plot roughly makes sense.
#' Yes, this is a wrapper for a regular sequence generating function, but it takes similar inputs as generateStrandModel and more importantly adds context to the
#' parameter list.
#'
#' @param startBase The starting base pair number
#' @param endBase The ending base pair number
#' @param scaling scaling factor in how many base pairs should be represented per index (ex. 1000) to reduce computational complexity when simulating physics/plotting;
#' Should match value used in generateStrandModel
#'
#' @return A sequence containing the location model.
#' @export
#'
#' @examples
#'\dontrun{
#'    generateLocationModel(startBase=1, endBase=195471971, scaling=100000)
#'}
#

generateLocationModel <- function(startBase, endBase, scaling=1000)
{
  locationModel <- t(as.matrix(seq(from=startBase, to=endBase + scaling - 1, by=scaling)))
  return(locationModel)
  #return(format(t(as.matrix(seq(from=startBase, to=endBase + scaling, by=scaling))), scientific = FALSE))
}
