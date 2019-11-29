# hotgenes

## Description

The goal of hotgenes is to take as input RNA-seq data (aligned with Rsubread and counted using featureCounts) and create a plot of the genome of a selected area of the data. The visualization is adjusted using a simulation of heat conduction to fill in the gaps between exons to give a rough idea where the 'hotspots' of expression are in the selected area.

## Installation

You can install the latest version of hotgenes from github using the following code:

``` r
library(devtools)
install_github("hyf97ca/hotgenes",  build_vignettes = TRUE)
library(hotgenes)
```

## Overview
![](./inst/extdata/Pitch.PNG)

  This package contains three functions that form the end of a pipeline for displaying 

![](./inst/extdata/musCh1.PNG)

## Contributions

The author of the package is Yi Fei Huang. All code in this package was written from scratch by Yi Fei Huang. The TPM calculation used to normalize the counts was originally proposed by Wagner GP, Kin K, Lynch VJ. The package assumes preprocessing of data by Rsubread and its featurecounts subcomponent. It directly takes in a featurecounts object provided by Rsubread as part of the input. Inspiration for package structure and comments was drawn from Rppt by Boris Steipe. The example data was constructed using instructions from Phipson, B., Doyle, M., & Dashnow, H. using raw data from Fu, N. Y., Rios, A. C., Pal, B., Soetanto, R. et al.

## References

 Wagner GP, Kin K, Lynch VJ. Measurement of mRNA abundance using RNA-seq data: RPKM measure
    # is inconsistent among samples. Theory in biosciences. 2012 Dec 1;131(4):281-5.

Fu, N. Y., Rios, A. C., Pal, B., Soetanto, R., Lun, A. T. L., Liu, K., . . . Visvader, J. E. (2015). EGF-mediated
 induction of Mcl-1 at the switch to lactation is essential for alveolar cell survival. Nature Cell Biology, 17(4),
 365–375. doi: 10.1038/ncb3117

## Example
(run rsubread align and featureCounts on your sequenced RNA to get a featureCounts data structure)
``` r
library(hotgenes)
data(musCh1fc)
sm <- generateStrandModels(1, 195471971, musCh1fc, "chr1", "-", 100000)
x <- generateLocationModel(1, 195471971, 100000)
sm <- simulateHeatSpread(sm, 0.001, 1000)
plotHeatedMap(sm, x)
```
