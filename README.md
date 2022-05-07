# JACKIE2
```
     ██╗ █████╗  ██████╗██╗  ██╗██╗███████╗
     ██║██╔══██╗██╔════╝██║ ██╔╝██║██╔════╝
     ██║███████║██║     █████╔╝ ██║█████╗  
██   ██║██╔══██║██║     ██╔═██╗ ██║██╔══╝  
╚█████╔╝██║  ██║╚██████╗██║  ██╗██║███████╗
 ╚════╝ ╚═╝  ╚═╝ ╚═════╝╚═╝  ╚═╝╚═╝╚══════╝
```

JACKIE (Jackie and Albert's CRISPR K-mer Instances Enumerator) [yes, a recursive acronym!], is a software pipeline written mainly in C++ with accessory scripts written in bash and python languges, that allows enumeration of all kmer (e.g., SpCas9 binding sites, ZFPs) in a genome and output their sequences, copy numbers and locations. We have demonstrated its application to design sgRNA with clustered repetitive binding sites for imaging genomic region. Please see the [preprint](https://doi.org/10.1101/2020.02.27.968933) for more details. JACKIE2 is a new version that incorporates the fast computation of off-target sequence neighbors, allowing for example >50 million sequences to be checked for off-target sequence neighbors up to 3 mismatches in less than 12 hours.

## Precomputed CRISPR sites (JACKIEdb)
Precomputed CRISPR sites (JACKIEdb) available for hg38 and mm10, available at http://cheng.bio/JACKIE
Users interested in starting from precomputed sites to further filter or identify sites within regions, etc, can jump to the [later sections](README.md#filtering-examples)

## Install binaries

Binaries available for:
Architecture | Folder
--- | ---
M1 Mac | [Binaries/arm64-apple-darwin20](Binaries/arm64-apple-darwin20)
Intel Mac x86_64 | [Binaries/x86_64-Mac](Binaries/x86_64-Mac)
Linux x86_64 | [Binaries/x86_64-linux](Binaries/x86_64-linux)

Installation of binaries
with root privilege

```
cp Binaries/arm64-apple-darwin20/* /usr/bin/ #for M1 Mac, change according to above table
```

without root privilege
```
mkdir ~/bin/
cp Binaries/arm64-apple-darwin20/* ~/bin/ #for M1 Mac, change according to above table
```
Add JACKIE to path. In your ~/.bashrc file, add a line:
```
export PATH=~/bin/:${PATH}
```

## Build from source

With root privilege:

```
./configure
make
make install
```

Install to specific path:

```
./configure --prefix=/path/to/install/
make
make install
```
Build JACKIE.queryDB, see [JACKIE.queryDB](JACKIE.queryDB)

Add JACKIE to path. In your ~/.bashrc file, add a line:
```
export PATH=/path/to/install/:${PATH}
```
## Example run
Download genome fasta files and produce a merged files for "non-random" chromosomes
```
genome=<fill in your genome> #e.g., hg38
genomesRoot=<fill in your genomes data root path> 
pathToGenome=$genomesRoot/$genome
genomeFasta=$pathToGenome/$genome.nr.fa

cd $pathToGenome
wget --timestamping "ftp://hgdownload.cse.ucsc.edu/goldenPath/$genome/chromosomes/*"
gunzip *.gz
mkdir random
mv *random*.fa random/
mv chrUn*.fa random/
mkdir nr
mv *.fa nr
cat nr/*.fa > $genome.nr.fa

```
Run JACKIE2 using SLURM cluster

Parameters. Change these to fit your case.
```
#cluster script header, modify to fit your cluster. This is for Slurm. 256G memory is needed for encodeSeqSpace.
CLUSTER_SCRIPT_HEADER='#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem 256G
#SBATCH -t 0-18:00:00
#SBATCH -p compute
#SBATCH -q batch'

JACKIE_DIR=/path/to/JACKIE2/ #If JACKIE2 is added to $PATH, then you don't need to add ${JACKIE_DIR}/ in the codes below.
GENOME_DIR=/path/to/genome/
GENOME=hg38 #genome
kmer=20 #kmer
mm=3 #mismatches
pattern=XXXXXXXXXXXXXXXXXXXXNGG #pattern, can be XXXXXXXXXXXXXXXXXXXX for unrestricted 20mer
shortPattern=20NGG #shortened name for pattern
```
Create bin files XXXXXX.bin
```
mkdir ${GENOME_DIR}/${GENOME}/pamFold-${shortPattern}
cd ${GENOME_DIR}/${GENOME}/pamFold-${shortPattern}
for N in A C G T; do
echo "${CLUSTER_SCRIPT_HEADER}" > ${GENOME_DIR}/${GENOME}/pamFold-${shortPattern}/$N.b4.slurmjob.sh
echo "#SBATCH -e ${GENOME_DIR}/${GENOME}/pamFold-${shortPattern}/$N.b4.slurmjob.stderr" >> ${GENOME_DIR}/${GENOME}/pamFold-${shortPattern}/$N.b4.slurmjob.sh
echo "#SBATCH -o ${GENOME_DIR}/${GENOME}/pamFold-${shortPattern}/$N.b4.slurmjob.stdout" >> ${GENOME_DIR}/${GENOME}/pamFold-${shortPattern}/$N.b4.slurmjob.sh
echo "date; ${JACKIE_DIR}/JACKIE.bin ${GENOME_DIR}/${GENOME}/${GENOME}.nr.fa ${GENOME_DIR}/${GENOME}/pamFold-${shortPattern}/ .bin 6 ${GENOME_DIR}/${GENOME}/pamFold-${shortPattern}/$N.ref.txt $N n ${pattern}; date" >> ${GENOME_DIR}/${GENOME}/pamFold-${shortPattern}/$N.b4.slurmjob.sh
sbatch ${GENOME_DIR}/${GENOME}/pamFold-${shortPattern}/$N.b4.slurmjob.sh
done
```
Fold bin files to XXXXXX.bed, 16 parallel jobs, each focusing on 2bp prefix (4^2)
```
for prefix in AA AC AT AG CA CC CT CG TA TC TT TG GA GC GT GG; do
echo "${CLUSTER_SCRIPT_HEADER}" > ${GENOME_DIR}/${GENOME}/pamFold-${shortPattern}/$prefix.outbed.sh
echo "#SBATCH -e ${GENOME_DIR}/${GENOME}/pamFold-${shortPattern}/$prefix.outbed.stderr" >> ${GENOME_DIR}/${GENOME}/pamFold-${shortPattern}/$prefix.outbed.sh
echo "#SBATCH -o ${GENOME_DIR}/${GENOME}/pamFold-${shortPattern}/$prefix.outbed.stdout" >> ${GENOME_DIR}/${GENOME}/pamFold-${shortPattern}/$prefix.outbed.sh
echo "bash ${JACKIE_DIR}/outbedForPrefixJob.sh ${GENOME_DIR}/${GENOME}/pamFold-${shortPattern}/ $prefix" >> ${GENOME_DIR}/${GENOME}/pamFold-${shortPattern}/$prefix.outbed.sh
sbatch ${GENOME_DIR}/${GENOME}/pamFold-${shortPattern}/$prefix.outbed.sh
done
```
combine XXXXXX.bed to one BED for the whole genome
```
echo "${CLUSTER_SCRIPT_HEADER}" > ${GENOME_DIR}/${GENOME}/pamFold-${shortPattern}/combineBed.sh
echo "#SBATCH -e ${GENOME_DIR}/${GENOME}/pamFold-${shortPattern}/combineBed.stderr" >> ${GENOME_DIR}/${GENOME}/pamFold-${shortPattern}/combineBed.sh
echo "#SBATCH -o ${GENOME_DIR}/${GENOME}/pamFold-${shortPattern}/combineBed.stdout" >> ${GENOME_DIR}/${GENOME}/pamFold-${shortPattern}/combineBed.sh
echo "cat ${GENOME_DIR}/${GENOME}/pamFold-${shortPattern}/*.bed > ${GENOME_DIR}/${GENOME}/pamFold-${shortPattern}/${GENOME}.${shortPattern}.BED" >> ${GENOME_DIR}/${GENOME}/pamFold-${shortPattern}/combineBed.sh
sbatch ${GENOME_DIR}/${GENOME}/pamFold-${shortPattern}/combineBed.sh
```

Optional step - collapse CRISPR sites from same chromosome into an extended bed format, with each line containing sites of the same sequence:
```
#collapse sgRNA binding locations with same sequecnes into an extended bed format (exclusively within a chromosome)
chainSameChromSites.py ${GENOME_DIR}/${GENOME}/pamFold-${shortPattern}/${GENOME}.${shortPattern}.BED 0 2 > ${GENOME_DIR}/${GENOME}/pamFold-${shortPattern}/${GENOME}.${shortPattern}.sameChr.tx.bed
```

Optional step - sort extended bed file:
```
sort -k1,1 -k2,2n ${GENOME_DIR}/${GENOME}/pamFold-${shortPattern}/${GENOME}.${shortPattern}.sameChr.tx.bed > ${GENOME_DIR}/${GENOME}/pamFold-${shortPattern}/${GENOME}.${shortPattern}.sameChr.tx.sorted.bed

#some CRISPR sites are overlapping, and will crash browser visualization. Remove entries with overlapping sgRNA sites

```

Optional: If you want to compute off-target profiles.
Encode k-mer sequence NGG-subspace of the genome using `JACKIE.encodeSeqSpaceNGG`. If no NGG restriction, use `JACKIE.encodeSeqSpace`.
```
echo "${CLUSTER_SCRIPT_HEADER}" > ${GENOME_DIR}/${GENOME}/encodeSeqSpaceNGG.$kmer.slurmjob.sh
echo "#SBATCH -e ${GENOME_DIR}/${GENOME}/encodeSeqSpaceNGG.$kmer.slurmjob.stderr" >> ${GENOME_DIR}/${GENOME}/encodeSeqSpaceNGG.$kmer.slurmjob.sh
echo "#SBATCH -o ${GENOME_DIR}/${GENOME}/encodeSeqSpaceNGG.$kmer.slurmjob.stdout" >> ${GENOME_DIR}/${GENOME}/encodeSeqSpaceNGG.$kmer.slurmjob.sh
echo "date; ${JACKIE_DIR}/JACKIE.encodeSeqSpaceNGG ${GENOME_DIR}/${GENOME}/$GENOME.$kmer.NGG.seqbits.gz $kmer ${GENOME_DIR}/${GENOME}/nr/*.fa; date" >> ${GENOME_DIR}/${GENOME}/encodeSeqSpaceNGG.$kmer.slurmjob.sh
sbatch ${GENOME_DIR}/${GENOME}/encodeSeqSpaceNGG.$kmer.slurmjob.sh
```

JACKIE.encodeSeqSpace and JACKIE.encodeSeqSpaceNGG with kmer=20 will require 128GB memory (RAM virtual). If such memory is not available. 

A newer "divide-and-conquer" approach is available via JACKIE.encodeSeqSpace.prefixed
For example, the sequence space is divided into AA,AC,AT subspaces. This will require only 8GB memory.

```
PAM=NGG #PAM=- for no-PAM restrictions, i.e., all possible k-mers

for PREFIX in AA AC AT AG CA CC CT CG TA TC TT TG GA GC GT GG; do
date; date +%s;
echo "working on $PREFIX"
echo "date; date +%s; ${JACKIE_DIR}/JACKIE.encodeSeqSpace.prefixed ${GENOME_DIR}/${GENOME}/$GENOME.$kmer.$PREFIX.$PAM.seqbits.gz $kmer $PREFIX $PAM ${GENOME_DIR}/${GENOME}/nr/*.fa; date; date +%s;" > ${GENOME_DIR}/${GENOME}/encodeSeqSpace.prefixed.$kmer.$PREFIX.$PAM.slurmjob.sh
bash ${GENOME_DIR}/${GENOME}/encodeSeqSpace.prefixed.$kmer.$PREFIX.$PAM.slurmjob.sh > ${GENOME_DIR}/${GENOME}/encodeSeqSpace.prefixed.$kmer.$PREFIX.$PAM.slurmjob.stdout 2> ${GENOME_DIR}/${GENOME}/encodeSeqSpace.prefixed.$kmer.$PREFIX.$PAM.slurmjob.stderr
done;
date; date +%s;
```

Similarly, a "divide-and-conquer" approach for counting sequence neighbors:
```

gcRange=0.4,0.6
cpRange=1,1

PAM=NGG
inputBed=${GENOME_DIR}/${GENOME}/pamFold-${shortPattern}/${GENOME}.${shortPattern}.cpRange${cpRange}.GC${gcRange}.BED      #change to your input bed file name
outputBed=${GENOME_DIR}/${GENOME}/pamFold-${shortPattern}/${GENOME}.${shortPattern}.cpRange${cpRange}.GC${gcRange}.offProfile.BED  #change to your output bed file name
sequenceField="4,/,2" #fourth column => split by "/" => second element


echo "date; date +%s; ${JACKIE_DIR}/JACKIE.countSeqNeighbors.pmulti ${GENOME_DIR}/${GENOME}/$GENOME.$kmer.[AA,AC,AG,AT,CA,CC,CG,CT,GA,GC,GG,GT,TA,TC,TG,TT].$PAM.seqbits.gz $kmer $mm $inputBed $sequenceField > $outputBed ; date; date +%s;"  > ${GENOME_DIR}/${GENOME}/countSeqNeighbors.pmulti.$kmer.$PAM.slurmjob.sh

bash ${GENOME_DIR}/${GENOME}/countSeqNeighbors.pmulti.$kmer.$PAM.slurmjob.sh

```

To allow for exact number of target sites to be identified, use JACKIE.encodeSeqCountDatabase to encode a SeqCountDatabase.

```
numBitsPerSeq=3 #3-bits hold up to 2^3-2=6 in the bit array, and put the "overflow" count to map

for PREFIX in AA AC AT AG CA CC CT CG TA TC TT TG GA GC GT GG ; do
echo "working on $PREFIX"
date; date +%s;
JACKIE.encodeSeqCountDatabase ${GENOME_DIR}/${GENOME}/$GENOME.$kmer.$PREFIX.$PAM.${numBitsPerSeq}bits.seqbits.gz $kmer $PREFIX $PAM ${numBitsPerSeq} ${GENOME_DIR}/${GENOME}/nr/*.fa
date; date +%s;
done;

```

Use JACKIE.countOffSites to count number of off-targets per query

```

numBitsPerSeq=3 

inputBed=${GENOME_DIR}/${GENOME}/pamFold-${shortPattern}/${GENOME}.${shortPattern}.cpRange${cpRange}.GC${gcRange}.top1000000.BED  #change to your input bed file name
outputBed=${GENOME_DIR}/${GENOME}/pamFold-${shortPattern}/${GENOME}.${shortPattern}.cpRange${cpRange}.GC${gcRange}.top1000000.offProfile.pmulti.BED  #change to your output bed file name
sequenceField="4,/,2" #fourth column => split by "/" => second element"


JACKIE.countOffSites ${GENOME_DIR}/${GENOME}/$GENOME.$kmer.[AA,AC,AG,AT,CA,CC,CG,CT,GA,GC,GG,GT,TA,TC,TG,TT].$PAM.${numBitsPerSeq}bits.seqbits.gz $kmer $mm $inputBed $sequenceField > $outputBed
```

Make bigBed-formatted JACKIEdb according to the following AutoSql structures
For one-copy sites

```
table OneCopyPAMOffSiteCounts
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
We need to duplicate some fields (thickStart, thickEnd), generate the nXmismatches, totalOffsites, percentGC longestTandemT fields etc
The bed file we have now is structured.

Column | Example field
--- | ---
Column 1 | chr1
Column 2 | 10607
Column 3 | 10627
Column 4 | 8073627651782476488.1/ACGCAACTCCGCCGTTGCAA
Column 5 | 1
Column 6 | +
Column 7 | 1/2/21/26

```
SevenColumnBedFile= #fill input 7-column bed file
EighteenColumnBedFile= #fill output 18-column bed file

#convert 7-col to 18-col
awk -v FS="\t" -v OFS="\t" awk -v FS="\t" -v OFS="\t" $SevenColumnBedFile ' {split($4,a,"/"); seq=a[2]; offString=$(7); split(offString,m,"/"); $(7)=$(2); $(8)=$(3); $(9)="0,0,0"; $(10)=m[1]; $(11)=m[2]; $(12)=m[3]; $(13)=m[4]; $(14)=m[2]+m[3]+m[4]; $(15)=offString; $(16)=seq; TT=0; GC=0; maxT=0; for(i=1;i<=length(seq);i++){thisBase=substr(seq,i,1); if(thisBase=="T"){TT++;if(TT>maxT){maxT=TT;}}else{TT=0;} if(thisBase=="G" || thisBase=="C"){GC++;}}$(17)=GC/length(seq)*100; $(18)=maxT; print;}' > $EighteenColumnBedFile

#sort by coordinates in order for bedToBigBed to work
sort -k1,1 -k2,2n $EighteenColumnBedFile > ${EighteenColumnBedFile/.bed/}.sorted.bed

#convert to bigBed
bedToBigBed -as=hg38PAM.1copy.offSiteCounts.as -type=bed9+9 ${EighteenColumnBedFile/.bed/}.sorted.bed hg38.chrom.sizes ${EighteenColumnBedFile/.bed/}.sorted.bb
```

and for 2+ copies clusters
```
table hg38TwoPlusCopyPAMOffSiteCounts
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
int blockCount; "Number of blocks"
int[blockCount] blockSizes; "Comma separated list of block sizes"
int[blockCount] chromStarts; "Start positions relative to chromStart"
uint    clusterSize;    "Size of gRNA cluster"
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

## For filtering, selection of best N gRNAs within specified genomic regions with JACKIE.queryDB, see [JACKIE.queryDB](JACKIE.queryDB)



