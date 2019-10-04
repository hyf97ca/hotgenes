#' \code{generateStrandModels}
#'
#' @param startBase The starting base pair number
#'
#' @param endBase The ending base pair number
#' @param fc The featurecounts object from rsubread after counting your alignment
#' @param strand \+ or \- strand to generate models for
#'
#' @return A data frame containing the strand model.
#' @export
#'
# @examples
# \dontrun
# {
#     fc <- featureCounts(files=bamFile, annot.ext=annotation, useMetaFeatures=FALSE)
# }
#
generateStrandModels <- function(startBase, endBase, fc, strand)
{
  fcLength = length(fc[["counts"]][,1])
  numDatasets = length(fc[["counts"]][1,])

  models <- list()

  for (i in 1:numDatasets)
  {
    #print(c("Starting ", fc[["targets"]][i]))
    genome <- vector(mode = "numeric", length = (endBase - startBase + 1))
    #transcripts per kilobase million sum value
    tpmSum <- 0
    for (j in 1:fcLength)
    {
      entry <- fc[["annotation"]][j,]
      if (entry[["Strand"]] == strand
       & entry[["Start"]] >= startBase
       & entry[["End"]] <= endBase)
      {
        #Entry within scope of view

        #reads per kilobase
        rpk <- fc[["counts"]][j,i] / ((entry[["Length"]])/1000)
        tpmSum <- tpmSum + rpk
        for (k in entry[["Start"]]:entry[["End"]])
        {
          genome[k - startBase + 1] <- genome[k - startBase + 1] + rpk
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
    print(c(fc[["targets"]][i], i, "done"))
    models[[i]] <- genome
    #names(models[[i]]) <- fc[["targets"]][i]
  }
  return(as.data.frame(models, col.names = fc[["targets"]]))
}
