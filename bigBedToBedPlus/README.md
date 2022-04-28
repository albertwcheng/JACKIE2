# bigBedToBedPlus

bigBedToBedPlus is a modified version of Kent's bigBedToBed which allows the extraction of multiple regions simultaneously, as well as allowing for filtering of integer-based fields of custom bigBed format. Also allows for the selection of a specified number of best bed items according to a series of column sorting criteria.

For compilation, skip to [Building bigBedToBedPlus](#building-bigBedToBedPlus)

Usage:

The genomic regions of interest are supplied in a txt file and supplied as an ```-filter filter.txt``` option to the program
The filter.txt file should contains each line chr <tab> start <tab> end, e.g.,

```
chr1	1	500000
chr2	1	1000000

```

# Building bigBedToBedPlus
To compile this, place this folder under kent source tree under folder kent/src/utils/ (i.e., make this folder kent/src/utils/bigBedToBedPlus) , then run ```make``` on kent/src level first, then go into bigBedToBedPlus and run ```make```.

