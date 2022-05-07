table OneCopyOffSiteCounts
"One-copy sites and off-target counts up to 3 mismatches"
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
