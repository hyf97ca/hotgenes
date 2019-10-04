#' \code{generateStrandModels}
#'
#' @param startBase The starting base pair number
#' @param endBase The ending base pair number
#' @param fc The featurecounts object from rsubread after counting your alignment
#' @param chr Chromosome; Should match Chr in fc
#' @param strand \+ or \- strand to generate models for
#' @param scaling scaling factor in how many base pairs should be represented per index (ex. 1000) to reduce computational complexity when simulating physics/plotting
#'
#' @return A matrix containing the strand model.
#' @export
#'
# @examples
# \dontrun
# {
#     fc <- featureCounts(files=bamFile, annot.ext=annotation, useMetaFeatures=FALSE)
#     sm <- generateStrandModels(1, 195471971, fc, "Ch1", "-")
# }
#
generateStrandModels <- function(startBase, endBase, fc, chr, strand, scaling=1000)
{
  fcLength = length(fc[["counts"]][,1])
  numDatasets = length(fc[["counts"]][1,])

  models <- list()

  for (i in 1:numDatasets)
  {
    message("Starting processing of chromsome ", chr, " strand ", strand, " in ", fc[["targets"]][i])
    genome <- vector(mode = "numeric", length = ceiling((endBase - startBase)/scaling) + 1)
    #transcripts per kilobase million sum value
    tpmSum <- 0

    prog <- txtProgressBar(min = 0, max = fcLength, style=3)
    for (j in 1:fcLength)
    {
      entry <- fc[["annotation"]][j,]
      setTxtProgressBar(prog, j)
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
        for (k in startInVector:endInVector)
        {
          if (k == startInVector & (entry[["Start"]] - startBase)/scaling > k)
          {
            genome[k] <- genome[k] + rpk * ((entry[["Start"]] - startBase) - startInVector*scaling)/scaling
          }
          else if (k == endInVector & (entry[["End"]] - startBase)/scaling < k)
          {
            genome[k] <- genome[k] + rpk * (endInVector*scaling - (entry[["End"]] - startBase))/scaling
          }
          else
          {
            genome[k] <- genome[k] + rpk
          }
        }
      }

    }

    #normalize dataset to TPM + 0.1
    tpmDivFactor <- tpmSum / 1000000

    #don't divide by 0
    if (tpmDivFactor == 0)
    {
      tpmDivFactor <- 1
    }

    genome <- (genome/tpmDivFactor) + rep_len(x=0.1, length.out=length(genome))

    #save into models
    #names(genome) <- fc[["targets"]][i]
    message("\n", fc[["targets"]][i], i, " done")
    models[[i]] <- genome
    #names(models[[i]]) <- fc[["targets"]][i]

  }
  return(t(as.data.frame(models)))
}


#' @rdname generateStrandModels
#' @export
# @examples
generateStrandModel <- generateStrandModels
