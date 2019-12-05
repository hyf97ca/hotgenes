#generateStrandModel.R

#' @title Create a Strand Model
#'
#' generateStrandModel Creates a model of a DNA strand in the form of a matrix (data frames are too slow) for further processing/plotting.
#'
#' @name generateStrandModel
#' @param startBase The starting base pair number
#' @param endBase The ending base pair number
#' @param fc The featurecounts object from rsubread after counting your alignment
#' @param chr Chromosome; Should match Chr in fc
#' @param strand \+ or \- strand to generate models for
#' @param scaling scaling factor in how many base pairs should be represented per index (ex. 1000) to
#'                reduce computational complexity when simulating physics/plotting
#' @param updateProgressBar a function for linking with shiny's progress bar. This will be called if not NULL
#'                          when progress bar is to be updated (every 1000 processed features)
#'
#' @return A matrix containing the strand model.
#' @import utils
#' @export
#'
#' @examples
#' \dontrun{
#'    #featureCounts() is from RSubread package; only run this if that is installed
#'    fc <- featureCounts(files=bamFile, annot.ext=annotation, useMetaFeatures=FALSE)
#'    #replace fc with the one loaded from featureCounts if you are using your own data
#'    sm <- generateStrandModel(startBase=1, endBase=195471971,
#'     fc=musCh1fc, chr="chr1", strand="-", scaling=100000)
#'    sm1 <- generateStrandModel(1, 195471971, musCh1fc, "chr1", "+")
#'    sm2 <- generateStrandModel(3000000, 3217000, musCh1fc, "chr1", "-")
#' }
#'

generateStrandModel <- function(startBase, endBase, fc, chr, strand, scaling=1000, updateProgressBar=NULL)
{
  if (!(startBase < endBase & !is.null(fc) & !is.null(chr) & !is.null(strand)))
    stop("Invalid parameters")

  fcLength = length(fc[["counts"]][,1])
  numDatasets = length(fc[["counts"]][1,])

  models <- list()

  for (i in 1:numDatasets)
  {
    if (!shiny::isRunning())
    {
      message("Starting processing of chromsome ", chr, " strand ", strand, " in ", fc[["targets"]][i])
    }
    genome <- vector(mode = "numeric", length = ceiling((endBase - startBase)/scaling) + 1)

    #transcripts per kilobase million sum value
    # TPM, the idea, was originally presented in the following paper:
    # Wagner GP, Kin K, Lynch VJ. Measurement of mRNA abundance using RNA-seq data: RPKM measure
    # is inconsistent among samples. Theory in biosciences. 2012 Dec 1;131(4):281-5.
    tpmSum <- 0
    if (!shiny::isRunning())
    {
      prog <- txtProgressBar(min = 0, max = fcLength, style=3)
    }
    for (j in 1:fcLength)
    {
      entry <- fc[["annotation"]][j,]
      if (!shiny::isRunning())
      {
        setTxtProgressBar(prog, j)
      }
      else if (!is.null(updateProgressBar) & j %% 1000 == 0)
      {
        updateProgressBar()
      }
      #reads per kilobase
      rpk <- fc[["counts"]][j,i] / ((entry[["Length"]])/1000)
      tpmSum <- tpmSum + rpk

      if (entry[["Strand"]] == strand
          & entry[["Chr"]] == chr
          & entry[["Start"]] >= startBase
          & entry[["End"]] <= endBase)
      {
        #Entry within scope of view

        startInVector <- floor((entry[["Start"]] - startBase)/scaling)
        endInVector <- ceiling((entry[["End"]] - startBase)/scaling)
        if (startInVector == 0)
          startInVector <- 1
        for (k in startInVector : endInVector)
        {
          #case where the feature only covers end half of a cell (ie. starting end)
          if (k == startInVector & (entry[["Start"]] - startBase)/scaling > k & entry[["Length"]] > scaling)
          {
            genome[k] <- genome[k] + rpk * ((entry[["Start"]] - startBase) - startInVector*scaling)/(scaling*(endInVector-startInVector))
          }
          #case where feature only covers start half of a cell (ie. ending end)
          else if (k == endInVector & (entry[["End"]] - startBase)/scaling < k & entry[["Length"]] > scaling)
          {
            genome[k] <- genome[k] + rpk * (endInVector*scaling - (entry[["End"]] - startBase))/(scaling*(endInVector-startInVector))
          }
          #general case; feature covers entire cell or is entirely contained
          else
          {
            genome[k] <- genome[k] + rpk/(endInVector-startInVector)
          }
        }
      }

    }

    # normalize dataset to TPM + 0.1
    tpmDivFactor <- tpmSum / 1000000

    #don't divide by 0
    if (tpmDivFactor == 0)
    {
      tpmDivFactor <- 1
    }

    genome <- (genome/tpmDivFactor) + rep_len(x=0.1, length.out = length(genome))

    #save into models
    if (!shiny::isRunning())
    {
      message("\n", fc[["targets"]][i], i, " done")
    }
    models[[i]] <- genome

  }
  t(as.data.frame(models))
}

#' Rebuild a Strand Model
#'
#' This function is cheaper than rebuilding a brand new strand model from scratch. These savings are realized by using parts of the old strand model.
#'
#' @param strandModel a strand model that is a superset of the new strand model with tighter scaling
#'
#' @param newStartBase an integer depicting the desired start base pair number.
#' @param newEndBase an integer depicting the desired end base pair number.
#' @param newScaling scaling factor in how many base pairs should be represented per index. Should be larger (eg. 1000000) than old scaling factor.
#'                   IMPORTANT: newScaling is assumed to be divisible by scaling.
#'
#' @param startBase an integer depicting the old strand model's starting base pair number. Should be the beginning of the chromosome to be rendered (eg. 1).
#' @param endBase an integer depicting the old strand model's starting base pair number. Should be the end of the chromosome to be rendered (eg. 3000000).
#' @param scaling scaling factor in how many base pairs should be represented per index in the existing strand model. Should be smaller than any rebuilt choice (eg. 1)
#'
#' @return a strand model with the new start-end base and/or more relaxed scaling
#' @export
#'
#' @examples
#' \dontrun{
#'    sm <- generateStrandModel(startBase=1, endBase=195471971,
#'     fc=musCh1fc, chr="chr1", strand="-", scaling=1)
#'    sm1 <- rebuildStrandModel(strandModel=sm, newStartBase=1000,
#'     newEndBase=2000, newScaling=100000, startBase=1, endBase=195471971, scaling=1)
#' }
#'
rebuildStrandModel <- function(strandModel, newStartBase, newEndBase, newScaling, startBase, endBase, scaling)
{
  if (!(newScaling %% scaling == 0 & startBase < endBase & newStartBase < newEndBase & startBase <= newStartBase & endBase >= newEndBase))
    stop("Invalid parameters")

  rebuilt <- list()

  #the start of the section to be read in strandModel
  startIndex <- ceiling((newStartBase - startBase + 1)/scaling)
  #when to stop reading strandModel
  endIndex <- ceiling((newEndBase - startBase + 1)/scaling)
  #how many cells to combine into the new cell.
  scale <- newScaling / scaling

  for (columnIndex in 1:(nrow(strandModel)))
  {
    strand <- strandModel[columnIndex,]
    newStrand <- vector(mode = "numeric", length = ceiling((newEndBase - newStartBase)/newScaling) + 1)
    strandCounter <- 1
    for (rowIndex in startIndex:endIndex)
    {
      newStrandIndex <- ceiling(strandCounter/scale)
      #map old strand to new strand, remove the 0.1 from old strand that we added in generateStrandModel
      newStrand[newStrandIndex] <- newStrand[newStrandIndex] + strand[rowIndex] - 0.1
      strandCounter <- strandCounter + 1
    }
    #Re-add 0.1 to every index
    newStrand <- newStrand + rep_len(x=0.1, length.out = length(newStrand))
    rebuilt[[columnIndex]] <- newStrand
  }
  t(as.data.frame(rebuilt))
}

