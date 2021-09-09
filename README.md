# JACKIE

JACKIE (Jackie and Albert's CRISPR K-mer Instances Enumerator) [yes, a recursive acronym!], is a software pipeline written mainly in C++ with accessory scripts written in bash and python languges, that allows enumeration of all potential SpCas9 binding sites in a genome and output their sequences, copy numbers and locations. We have demonstrated its application to design sgRNA with clustered repetitive binding sites for imaging genomic region. Please see the [preprint](https://doi.org/10.1101/2020.02.27.968933) for more details.

## Precomputed CRISPR sites (JACKIEdb)
Precomputed CRISPR sites (JACKIEdb) available for hg38 and mm10, available at http://crispr.software/JACKIE
Users interested in strating from precomputed sites to further filter or identify sites within regions, etc, can jump to the [later sections](README.md#filtering-examples)

## Installation

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

Run first step of JACKIE, assuming your cluster uses `qsub`:

```
genome=<fill in your genome> #e.g., hg38
genomesRoot=<fill in your genomes data root path> 
pathToGenome=$genomesRoot/$genome
genomeFasta=$pathToGenome/$genome.nr.fa

jackieDB=$pathToGenome/jackieDB/
mkdir $jackieDB

#generate binary represetation of sgRNA binding locations

for N in A C G T; do
echo "date; JACKIE -b2 $genomeFasta $jackieDB .bin 6 $jackieDB/$N.ref.txt $N n; date" | qsub -l walltime=48:00:00
done

```
Please make sure all four jobs for step one have completed successfully, e.g., by using `qstat` and checking output files.
Second step, load binary files (*.bin), sort by sequence, output bed files:
```
#output bed file from binary files.
for prefix in AA AC AT AG CA CC CT CG TA TC TT TG GA GC GT GG; do
echo "date; outbedForPrefixJob.sh $jackieDB $prefix 1 0; date" | qsub -l walltime=24:00:00 -e `pwd`/$prefix.stderr.txt -o `pwd`/$prefix.stdout.txt	
done
```
Please make sure all jobs in second step have completed successfully.
Third step, merge all bed files into one:
```
#concatenate all bed files into one
echo "cat $jackieDB/*.bed > $jackieDB/${genome}PAM.BED" | qsub -l walltime=24:00:00
```

Optional step - collapse CRISPR sites from same chromosome into an extended bed format, with each line containing sites of the same sequence:
```
#collapse sgRNA binding locations with same sequecnes into an extended bed format
echo "chainExonBedsToTranscriptBed.py $jackieDB/${genome}PAM.BED 0 > $jackieDB/${genome}PAM.sameChr.tx.bed" | qsub -l walltime=24:00:00
```
Optional step - sort extended bed file:
```
#sort extended bed file
echo "sort -k1,1 -k2,2n $jackieDB/${genome}PAM.sameChr.tx.bed > $jackieDB/${genome}PAM.sameChr.tx.sorted.bed" | qsub -l walltime=24:00:00

#some CRISPR sites are overlapping, and will crash browser visualization. Remove entries with overlapping sgRNA sites
echo "removeIllegalBlockEntries.py $jackieDB/${genome}PAM.sameChr.tx.sorted.bed $jackieDB/${genome}PAM.sameChr.tx.sorted.legal.bed $jackieDB/${genome}PAM.sameChr.tx.sorted.illegal.bed" | qsub -l walltime=24:00:00
```

## Filtering examples

Select clustered sgRNA with (minBS)5 to (maxBS)8 binding sites and within (minDist)5kb to (maxDist)10kb distance
Precomputed [hg38 same-chromosome sites](http://albertcheng.info/jackie_downloads/hg38PAM.sameChr.tx.sorted.legal.bed.gz) and [mm10 same-chromosome sites](http://albertcheng.info/jackie_downloads/mm10PAM.sameChr.tx.sorted.legal.bed.gz)
```
#select clustered sgRNA with (minBS)5 to (maxBS)8 binding sites and within (minDist)5kb to (maxDist)10kb distance
minBS=5
maxBS=8
minDist=5000
maxDist=10000
genome=hg38 #human genome example
awk -v FS="\t" -v OFS="\t" -v minBS=$minBS -v maxBS=$maxBS -v minDist=$minDist -v maxDist=$maxDist '($3-$2>=minDist && $3-$2<=maxDist && $5>=minBS && $5<=maxBS)' $genome.sameChr.tx.sorted.legal.bed > $genome.sameChr.tx.sorted.legal.bed
```
Select unique sgRNA sites
```
#select unique sgRNA sites
awk -v FS="\t" -v OFS="\t" '($5==1)' $jackieDB/${genome}PAM.BED > $jackieDB/${genome}PAM.1copy.BED
```
Precomputed [hg38 1-copy sites](http://albertcheng.info/jackie_downloads/hg38PAM.1copy.BED.gz) and [mm10 1-copy sites](http://albertcheng.info/jackie_downloads/mm10PAM.1copy.BED.gz)

## Selecting (clustered or single) CRISPR sites overlapping regions of interest
Put regions of interest into a bed file, say,  `selection.bed`:
```
chr13	75845001	75850000	Region_A
chr13	76493001	76498000	Region_B
```
Overlap selection.bed with one copy sites:
```
fastjoinBedByOverlap.py selection.bed  $jackieDB/hg38PAM.1copy.BED > selection.overlap.hg38PAM.1copy.BED
```
## Run Cas-OFFinder 
Requires offline version of Cas-OFFinder at http://www.rgenome.net/cas-offinder/portable)
Also, cas-offinder should be in `$PATH`

For example, from selection.overlap.hg38PAM.1copy.BED above:
Let say, up to 3 mismatches. Sequence is encoded in the itemName (on the 8th column, second component of a split with "/"), so seqColExtract=8,/,2
```

genome=<fill in your genome> #e.g., hg38
genomesRoot=<fill in your genomes data root path> 
pathToGenome=$genomesRoot/$genome

BEDFile=selection.overlap.hg38PAM.1copy.BED
cas_outDir=/path/to/IOForCasOffFinder
maxNmismatches=3
seqColExtract=8,/,2
runCasOFFinderOnSequences.py $BEDFile $seqColExtract $maxNmismatches $pathToGenome $cas_outDir > $BEDFile.cas_off.txt
```

<!--


#run CasOffFinder (requires offline version of Cas-OFFinder at http://www.rgenome.net/cas-offinder/portable)
#export PATH=/usr/bin/:${PATH}
#export PATH=~/Dropbox/unixEnv/scripts:${PATH}
runCasOFFinderOnSequences.py <file> 8,/,2 3 ~/Dropbox/unixEnv/genomes/hg38/ casOffinder_outputDir > ???

#runCasOFFinderOnSequences.py newSelectionLoop.overlap.hg38PAM.sameChr.tx.sorted.legal.1copy.GC40to60.no5T.noLowercase.2.bed 17 3 ~/Dropbox/unixEnv/genomes/hg38/ newSelectionLoop.overlap.hg38PAM_off > newSelectionLoop.overlap.hg38PAM.sameChr.tx.sorted.legal.1copy.GC40to60.no5T.noLowercase.2.casOffinder.txt
```
-->
