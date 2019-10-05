#' \code{simulateHeatSpread}
#'
#' Simulates heat spread through conductance on (a) strand model(s). Kind of looks like anti-aliasing. Makes smaller expressed exons more visible.
#'
#' @name simulateHeatSpread
#' @param .strandModel Strand model from generateStrandModel
#' @param conductivity numeric between 0-1: 1 implies perfect equalization of heat between two adjacent cells in 1 iteration, 0 implies no heat transfer
#' @param iterations integer number of cycles to simulate: should be 1 or greater
#'
#' @return simulated strand model
#' @import utils
#' @export
#'
#' @examples
#' \dontrun{
#' sm <- generateStrandModels(1, 195471971, musCh1fc, "chr1", "-", 100000)
#' smh <- simulateHeatSpread(sm, 0.001, 1000)
#' smh1 <- simulateHeatSpread(sm, 1, 100)
#' }
#'

simulateHeatSpread <- function(.strandModel, conductivity, iterations)
{
  message("simulating heat spread with conductivity ", conductivity, " for ", iterations, " iterations...")
  prog <- txtProgressBar(min = 0, max = iterations, style=3)
  for (it in 1:iterations)
  {
    setTxtProgressBar(prog, it)
    for (columnIndex in 1:(nrow(.strandModel)))
    {
      strand <- .strandModel[columnIndex,]

      currentMaxIndex <- which.max(strand)
      while(strand[currentMaxIndex] > 0.1)
      {
        #do physics on two sides of actual strand then negate
        if (currentMaxIndex > 1)
          if(.strandModel[columnIndex, currentMaxIndex] > .strandModel[columnIndex, currentMaxIndex - 1])
          {
            transfer = (.strandModel[columnIndex, currentMaxIndex] - .strandModel[columnIndex, currentMaxIndex - 1])/2 * conductivity
            .strandModel[columnIndex, currentMaxIndex] = .strandModel[columnIndex, currentMaxIndex] - transfer
            .strandModel[columnIndex, currentMaxIndex - 1] = .strandModel[columnIndex, currentMaxIndex - 1] + transfer
          }
        if (currentMaxIndex < length(strand))
          if (.strandModel[columnIndex, currentMaxIndex] > .strandModel[columnIndex, currentMaxIndex + 1])
          {
            transfer = (.strandModel[columnIndex, currentMaxIndex] - .strandModel[columnIndex, currentMaxIndex + 1])/2 * conductivity
            .strandModel[columnIndex, currentMaxIndex] = .strandModel[columnIndex, currentMaxIndex] - transfer
            .strandModel[columnIndex, currentMaxIndex + 1] = .strandModel[columnIndex, currentMaxIndex + 1] + transfer
          }
        #negate so next max is chosen
        strand[currentMaxIndex] <- -strand[currentMaxIndex]
        currentMaxIndex <- which.max(strand)
      }
    }
  }
  return(.strandModel)
}
