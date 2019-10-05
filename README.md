
# hotgenes

<!-- badges: start -->
<!-- badges: end -->
![](./inst/extdata/musCh1.PNG)
The goal of hotgenes is to take as input RNA-seq data (aligned with Rsubread and counted using featureCounts) and create a plot of the genome of a selected area of the data. The visualization is adjusted using a simulation of heat conduction to fill in the gaps between exons to give a rough idea where the 'hotspots' of expression are in the selected area.

## Installation

You can install the latest version of hotgenes from github using the following code:

``` r
library(devtools)
installgithub("hyf97ca/hotgenes")
library(hotgenes)
```

## Pitch
![](./inst/extdata/Pitch.PNG)

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
