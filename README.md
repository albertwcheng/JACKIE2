# JACKIE2

JACKIE (Jackie and Albert's CRISPR K-mer Instances Enumerator) [yes, a recursive acronym!], is a software pipeline written mainly in C++ with accessory scripts written in bash and python languges, that allows enumeration of all kmer (e.g., SpCas9 binding sites, ZFPs) in a genome and output their sequences, copy numbers and locations. We have demonstrated its application to design sgRNA with clustered repetitive binding sites for imaging genomic region. Please see the [preprint](https://doi.org/10.1101/2020.02.27.968933) for more details. JACKIE2 is a new version that incorporates the fast computation of off-target sequence neighbors, allowing for example >50 million sequences to be checked for off-target sequence neighbors up to 3 mismatches in less than 12 hours.

## Precomputed CRISPR sites (JACKIEdb)
Precomputed CRISPR sites (JACKIEdb) available for hg38 and mm10, available at http://cheng.bio/JACKIE
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
#collapse sgRNA binding locations with same sequecnes into an extended bed format
chainExonBedsToTranscriptBed.py ${GENOME_DIR}/${GENOME}/pamFold-${shortPattern}/${GENOME}.${shortPattern}.BED 0 > ${GENOME_DIR}/${GENOME}/pamFold-${shortPattern}/${GENOME}.${shortPattern}.sameChr.tx.bed
```

Optional step - sort extended bed file:
```
sort -k1,1 -k2,2n ${GENOME_DIR}/${GENOME}/pamFold-${shortPattern}/${GENOME}.${shortPattern}.sameChr.tx.bed > ${GENOME_DIR}/${GENOME}/pamFold-${shortPattern}/${GENOME}.${shortPattern}.sameChr.tx.sorted.bed

#some CRISPR sites are overlapping, and will crash browser visualization. Remove entries with overlapping sgRNA sites

removeIllegalBlockEntries.py ${GENOME_DIR}/${GENOME}/pamFold-${shortPattern}/${GENOME}.${shortPattern}.sameChr.tx.sorted.bed ${GENOME_DIR}/${GENOME}/pamFold-${shortPattern}/${GENOME}.${shortPattern}.sameChr.tx.sorted.legal.bed ${GENOME_DIR}/${GENOME}/pamFold-${shortPattern}/${GENOME}.${shortPattern}.sameChr.tx.sorted.illegal.bed"

```

Optional: filter for gc range 40%-60%, copy number = 1, max tandem run of T = 4
```
gcRange=0.4,0.6
cpRange=1,1
maxTandemT=4
echo "${CLUSTER_SCRIPT_HEADER}" > ${GENOME_DIR}/${GENOME}/pamFold-${shortPattern}/filterGC.slurmjob.sh
echo "#SBATCH -e ${GENOME_DIR}/${GENOME}/pamFold-${shortPattern}/filterGC.slurmjob.stderr" >> ${GENOME_DIR}/${GENOME}/pamFold-${shortPattern}/filterGC.slurmjob.sh
echo "#SBATCH -o ${GENOME_DIR}/${GENOME}/pamFold-${shortPattern}/filterGC.slurmjob.stdout" >> ${GENOME_DIR}/${GENOME}/pamFold-${shortPattern}/filterGC.slurmjob.sh
echo "python ${JACKIE_DIR}/filterPAMFoldGC.py --gc-range $gcRange --cp-range $cpRange --max-tandem-t $maxTandemT ${GENOME_DIR}/${GENOME}/pamFold-${shortPattern}/${GENOME}.${shortPattern}.BED > ${GENOME_DIR}/${GENOME}/pamFold-${shortPattern}/${GENOME}.${shortPattern}.cpRange${cpRange}.GC${gcRange}.BED" >> ${GENOME_DIR}/${GENOME}/pamFold-${shortPattern}/filterGC.slurmjob.sh
sbatch ${GENOME_DIR}/${GENOME}/pamFold-${shortPattern}/filterGC.slurmjob.sh
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




Optional example:
Generate 3-mismatch off-target profiles for sgRNA with gc range 40%-60%, copy number = 1, max tandem run of T = 4. piped into awk script to put result string into name of the bed file to preserve bed formatting. remove NucKey.
```
echo "${CLUSTER_SCRIPT_HEADER}" > ${GENOME_DIR}/${GENOME}/countSeqNeighbors.$kmer.slurmjob.sh
echo "#SBATCH -e ${GENOME_DIR}/${GENOME}/countSeqNeighbors.$kmer.slurmjob.stderr" >> ${GENOME_DIR}/${GENOME}/countSeqNeighbors.$kmer.slurmjob.sh
echo "#SBATCH -o ${GENOME_DIR}/${GENOME}/countSeqNeighbors.$kmer.slurmjob.stdout" >> ${GENOME_DIR}/${GENOME}/countSeqNeighbors.$kmer.slurmjob.sh
echo "date; ${JACKIE_DIR}/JACKIE.countSeqNeighbors ${GENOME_DIR}/${GENOME}/$GENOME.$kmer.NGG.seqbits.gz $kmer $mm ${GENOME_DIR}/${GENOME}/pamFold-${shortPattern}/${GENOME}.${shortPattern}.cpRange${cpRange}.GC${gcRange}.BED 4,/,2 | awk -v FS=\"\\t\" -v OFS=\"\\t\" '{split(\$4,a,\"/\"); \$4=a[2] \"/\" \$7; print \$1,\$2,\$3,\$4,\$5,\$6}' > ${GENOME_DIR}/${GENOME}/pamFold-${shortPattern}/${GENOME}.${shortPattern}.cpRange${cpRange}.GC${gcRange}.offProfile.BED ; date"  >> ${GENOME_DIR}/${GENOME}/countSeqNeighbors.$kmer.slurmjob.sh
sbatch ${GENOME_DIR}/${GENOME}/countSeqNeighbors.$kmer.slurmjob.sh
```


JACKIE.encodeSeqSpace and JACKIE.encodeSeqSpaceNGG with kmer=20 will require 128GB memory (RAM virtual). If such memory is not available. A newer "divide-and-conquer" approach is available via JACKIE.encodeSeqSpace.prefixed
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

PAM=NGG #PAM=- for no-PAM restrictions, i.e., all possible k-mers
inputBed=${GENOME_DIR}/${GENOME}/pamFold-${shortPattern}/${GENOME}.${shortPattern}.cpRange${cpRange}.GC${gcRange}.BED      #change to your input bed file name
outputBed=${GENOME_DIR}/${GENOME}/pamFold-${shortPattern}/${GENOME}.${shortPattern}.cpRange${cpRange}.GC${gcRange}.offProfile.BED  #change to your output bed file name
sequenceField="4,/,2" #fourth column => split by "/" => second element

date; date +%s;
PREFIX=AC
echo "working on $PREFIX"
echo "date; date +%s; ${JACKIE_DIR}/JACKIE.countSeqNeighbors.prefixed ${GENOME_DIR}/${GENOME}/$GENOME.$kmer.$PREFIX.$PAM.seqbits.gz $kmer $mm $inputBed $PREFIX 0 $sequenceField > tmp.00 ; date; date +%s;"  > ${GENOME_DIR}/${GENOME}/countSeqNeighbors.prefxied.$kmer.$PREFIX.$PAM.slurmjob.sh

bash ${GENOME_DIR}/${GENOME}/countSeqNeighbors.prefxied.$kmer.$PREFIX.$PAM.slurmjob.sh > ${GENOME_DIR}/${GENOME}/countSeqNeighbors.prefxied.$kmer.$PREFIX.$PAM.slurmjob.stdout 2> ${GENOME_DIR}/${GENOME}/countSeqNeighbors.prefxied.$kmer.$PREFIX.$PAM.slurmjob.stderr

mv tmp.00 tmp.01 

for PREFIX in AC AT AG CA CC CT CG TA TC TT TG GA GC GT GG; do
date; date +%s;
echo "working on $PREFIX"

echo "date; date +%s; ${JACKIE_DIR}/JACKIE.countSeqNeighbors.prefixed ${GENOME_DIR}/${GENOME}/$GENOME.$kmer.$PREFIX.$PAM.seqbits.gz $kmer $mm tmp.01 $PREFIX 1 $sequenceField > tmp.00 ; date; date +%s;"  > ${GENOME_DIR}/${GENOME}/countSeqNeighbors.prefxied.$kmer.$PREFIX.$PAM.slurmjob.sh

bash ${GENOME_DIR}/${GENOME}/countSeqNeighbors.prefxied.$kmer.$PREFIX.$PAM.slurmjob.sh > ${GENOME_DIR}/${GENOME}/countSeqNeighbors.prefxied.$kmer.$PREFIX.$PAM.slurmjob.stdout 2> ${GENOME_DIR}/${GENOME}/countSeqNeighbors.prefxied.$kmer.$PREFIX.$PAM.slurmjob.stderr

mv tmp.00 tmp.01 

done;

mv tmp.01 $outputBed

date; date +%s;



```




Optional example:
filter for 1/0/0/0 (sgRNAs with no off-target matches up to 3-mismatches)
```
echo "${CLUSTER_SCRIPT_HEADER}" > ${GENOME_DIR}/${GENOME}/filterOffProfile.slurmjob.sh
echo "#SBATCH -e ${GENOME_DIR}/${GENOME}/filterOffProfile.slurmjob.stderr" >> ${GENOME_DIR}/${GENOME}/filterOffProfile.slurmjob.sh
echo "#SBATCH -o ${GENOME_DIR}/${GENOME}/filterOffProfile.slurmjob.stdout" >> ${GENOME_DIR}/${GENOME}/filterOffProfile.slurmjob.sh
echo "date; grep 1/0/0/0 ${GENOME_DIR}/${GENOME}/pamFold-${shortPattern}/${GENOME}.${shortPattern}.cpRange${cpRange}.GC${gcRange}.offProfile.BED > ${GENOME_DIR}/${GENOME}/pamFold-${shortPattern}/${GENOME}.${shortPattern}.cpRange${cpRange}.GC${gcRange}.offp1000.BED; date"  >> ${GENOME_DIR}/${GENOME}/filterOffProfile.slurmjob.sh
sbatch ${GENOME_DIR}/${GENOME}/filterOffProfile.slurmjob.sh
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
Precomputed [hg38 1-copy sites](https://albertcheng.info/jackie_downloads/hg38PAM.1copy.BED.gz) and [mm10 1-copy sites](https://albertcheng.info/jackie_downloads/mm10PAM.1copy.BED.gz)

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

If you want to run Cas-OFFinder, follow these steps. If you are only interested in finding gRNAs without matches up to a certain mismatches, the `JACKIE.countSeqNeighbors` approach mentioned above is way faster when looking at millions of gRNAs.

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
