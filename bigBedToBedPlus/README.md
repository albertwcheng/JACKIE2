# bigBedToBedPlus

bigBedToBedPlus is a modified version of Kent's bigBedToBed which allows the extraction of multiple regions simultaneously, as well as allowing for filtering of integer-based fields of custom bigBed format. Also allows for the selection of a specified number of best bed items according to a series of column sorting criteria.

For compilation, skip to [Building bigBedToBedPlus](#building-bigBedToBedPlus)

Usage: bedBedToBedPlus [-filter <filterFile>] input.bb 
For convenient piping, I have decided to change output to STDOUT instead of writing to a file in bigBedToBed. To write to file, use redirection (i.e., bedBedToBedPlus [-filter <filterFile>] input.bb  > output.bed )


# Basic operations with examples
The genomic regions of interest are supplied in a txt file and supplied as an ```-filter filter.txt``` option to the program

Our example bigbed https://albertcheng.info/jackie_downloads/hg38PAM.1copy.offSiteCounts.wGCT.bb has the following structure (in autosql)
```
table hg38OneCopyPAMOffSiteCounts
"One-copy CRISPR SpCas9 NGG sites and off-target counts up to 3 mismatches"
(
string  chrom;		"Reference sequence chromosome or scaffold"
uint    chromStart;	"Start position of feature on chromosome"
uint    chromEnd;	"End position of feature on chromosome"
string  name;		"Name of gene"
uint    score;		"Score"
char[1] strand;		"+ or - for strand"
uint    thickStart;	"Coding region start"
uint    thickEnd;	"Coding region end"
uint  	reserved;	"Color"
uint    n0mismatches;  "Number of exact matched sites"
uint    n1mismatches;  "Number of 1-mismatch sites"
uint    n2mismatches;  "Number of 2-mismatch sites"
uint    n3mismatches;  "Number of 3-mismatch sites"
uint    totalOffSites; "Total number of 1,2,3-mismatch sites"
string  offSiteCounts; "String representation of offsite counts separated by slashes"
string  spacerSeq;  "Spacer sequence of gRNA"
uint    percentGC;  "Percent GC of spacer sequence"
uint    longestTandemT; "Longest run of T"
)
```
or in a table

Column Number | Name | Description
--- | --- | --- 
1 | chrom | Reference sequence chromosome or scaffold
2 | chromStart | Start position of feature on chromosome

The filter.txt file should contains each line chr <tab> start <tab> end, e.g.,
```
chr1	1	500000
chr2	1	1000000
```
Example: To run with filter file:
```
bedBedToBedPlus -filter filter.txt https://albertcheng.info/jackie_downloads/hg38PAM.1copy.offSiteCounts.wGCT.bb 
```
The filter.txt file can also contain the input.bb filename by adding :input.bb, and in command line only run with bedBedToBedPlus -filter filter.txt 
```
:https://albertcheng.info/jackie_downloads/hg38PAM.1copy.offSiteCounts.wGCT.bb
chr1	1	500000
chr2	1	1000000
```
# Filtering
The filter.txt can contain $column_number <tab> min <tab> max to filter items with column integer values within min and max inclusive. tab delimited fields 4+ are ignored and can be used for commenting, #COMMENTS. 
```
:https://albertcheng.info/jackie_downloads/hg38PAM.1copy.offSiteCounts.wGCT.bb
chr1	1	500000
chr2	1	1000000
$17	40	60	#filter for column 17 in range [40,60]
$14	0	4	#filter for column 14 in range [0,4]
$18	0	5	#filter for column 18 in range [0,5]
```

# Selecting best items
The filter.txt can direct the program to only print out the best N items according to a sort lists
```
:https://albertcheng.info/jackie_downloads/hg38PAM.1copy.offSiteCounts.wGCT.bb
chr1	1	500000
chr2	1	1000000
$17	40	60	#filter for column 17 in range [40,60]
$14	0	4	#filter for column 14 in range [0,4]
$18	0	5	#filter for column 18 in range [0,5]
!BEST	4	#print best 4 items from each region, according to the list of sorting belong
!MIN	14	#minimize column 14 first (ascending)
!MAX	17	#then maximize column 17 (descending)
```


# Building bigBedToBedPlus
To compile this, place this folder under kent source tree under folder kent/src/utils/ (i.e., make this folder kent/src/utils/bigBedToBedPlus) , then run ```make``` on kent/src level first, then go into bigBedToBedPlus and run ```make```.

