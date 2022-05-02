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
Users interested in strating from precomputed sites to further filter or identify sites within regions, etc, can jump to the [later sections](README.md#filtering-examples)

## Binaries

Binaries available for:
Architecture | Folder
--- | ---
Intel Mac | [Binaries/x86_64-Mac](Binaries/x86_64-Mac)
M1 Mac | [Binaries/arm64-apple-darwin20](Binaries/arm64-apple-darwin20)
Linux x86_64 | [Binaries/x86_64-linux](Binaries/x86_64-linux)

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

## For filtering, selection of best N gRNAs within specified genomic regions with JACKIE.queryDB, see [JACKIE.queryDB](JACKIE.queryDB)



