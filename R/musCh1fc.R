#' SRX681996: GSM1480302: Basal 2 dL #2; Mus musculus; RNA-Seq
#'
#' A sample dataset processed following instructions from \url{https://bioinformatics-core-shared-training.github.io/RNAseq-R/align-and-count.nb.html};
#' In the provided instructions truncated data was used; however this dataset uses the full data of the last sequencing run, although the index was constructed using only chr1.
#'
#'Phipson, B., Doyle, M., & Dashnow, H. (2017). RNA-seq analysis in R. Retrieved from https://bioinformatics-core-shared-training.github.io/RNAseq-R/align-and-count.nb.html
#'
#' @format A list of 4 objects
#' \describe{
#'   \item{counts}{count per exon}
#'   \item{annotation}{GeneID, Chr, Start, End, Strand, Length: annotation of exon, index is shared with that in counts}
#'   \item{targets}{dataset listing. Only 1 is used in this featureCounts output: SRR1552455.1.fastq.subread.BAM which delineates the pipeline it went through previously:
#'    fastq -> rsubread }
#'   \item{stat}{Statistics about this featureCounts run. Not used in this package, but is provided for completion.}
#' }
#' @source \url{https://www.ncbi.nlm.nih.gov/sra?term=SRX681996}
#' Fu, N. Y., Rios, A. C., Pal, B., Soetanto, R., Lun, A. T. L., Liu, K., . . . Visvader, J. E. (2015). EGF-mediated
#' induction of Mcl-1 at the switch to lactation is essential for alveolar cell survival. Nature Cell Biology, 17(4),
#' 365–375. doi: 10.1038/ncb3117
#'
"musCh1fc"
