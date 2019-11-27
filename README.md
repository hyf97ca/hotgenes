# hotgenes

<!-- badges: start -->
<!-- badges: end -->
![](./inst/extdata/musCh1.PNG)<br/>

## Installation

You can install the latest version of hotgenes from github using the following code:

``` r
library(devtools)
install_github("hyf97ca/hotgenes")
library(hotgenes)
```

## Overview
![](./inst/extdata/Pitch.PNG)

The goal of hotgenes is to take as input RNA-seq data (aligned with Rsubread and counted using featureCounts) and create a plot of the genome of a selected area of the data. The visualization is adjusted using a simulation of heat conduction to fill in the gaps between exons to give a rough idea where the 'hotspots' of expression are in the selected area.

## Contributions

All code in this package was written by Yi Fei Huang. The package assumes preprocessing of data by Rsubread and its featurecounts component and directly takes in a featurecounts object as provided by Rsubread as part of the input. Inspiration for package structure and comments was drawn from Rppt by Boris Steipe.

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
