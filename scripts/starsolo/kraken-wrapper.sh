#!/bin/bash
set -eo pipefail

#run with three positional arguments, like so:
SAMPLE=$1
REFERENCE=$2
CHEMISTRY=$3

#translate chemistry to actual starsolo parameters
if [[ $CHEMISTRY == "v2" ]]
then
        UMI=10
        STRAND=Forward
elif [[ $CHEMISTRY == "v3" ]]
then
        UMI=12
        STRAND=Forward
elif [[ $CHEMISTRY == "v2-5" ]]
then
        UMI=10
        STRAND=Reverse
        #set this to v2 for barcode purposes
        CHEMISTRY=v2
fi

#set up nested folder architecture to match /archive/HCA irods design
mkdir $SAMPLE && cd $SAMPLE && mkdir starsolo && cd starsolo

#start by getting the FASTQs
#this might be a CITE ID. split on "-" and take the first part
#if not CITE, this operation will change nothing
GEX=`echo $SAMPLE | sed 's/-.*//'`
bash /mnt/mapcloud/scripts/10x/utils/sample_fastq.sh $GEX
#collapse to single R1/R2. I1 not even ingested
mkdir fastq
cat *R1_001.fastq.gz > fastq/R1.fastq.gz
cat *R2_001.fastq.gz > fastq/R2.fastq.gz
rm *.fastq.gz

#actually run starsolo
~/STAR-2.7.3a/bin/Linux_x86_64/STAR --runThreadN `grep -c ^processor /proc/cpuinfo` \
	--genomeDir ~/STAR-2.7.3a/$REFERENCE \
	--readFilesIn fastq/R2.fastq.gz fastq/R1.fastq.gz \
	--readFilesCommand zcat \
	--outSAMtype BAM SortedByCoordinate \
	--outSAMmultNmax 1 \
	--outMultimapperOrder Random \
	--runRNGseed 1 \
	--outSAMattributes NH HI AS nM CB UB GX GN \
	--soloType CB_UMI_Simple \
	--soloCBwhitelist ~/STAR-2.7.3a/whitelists/$CHEMISTRY.txt \
	--soloBarcodeReadLength 0 \
	--soloUMIlen $UMI \
	--soloStrand $STRAND \
	--soloUMIfiltering MultiGeneUMI \
	--soloCellFilter CellRanger2.2 3000 0.99 10 \
	--soloFeatures Gene GeneFull Velocyto \
	--soloOutFileNames logs/ genes.tsv barcodes.tsv matrix.mtx \
	--outSAMunmapped Within

#okay, that's as much as we need from this. delete FASTQs for space conservation
rm -r fastq && cd ..

#kraken time!

#script: https://github.com/julianeweller/mg2sc
#DB: https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20201202.tar.gz
#hierarchy collapsing of output: /lustre/scratch117/cellgen/team205/kp9/jupyter/jw38/metagenomics-covid-check/

python3 /mnt/mg2sc/src/scMeG-kraken.py --input starsolo/Aligned.sortedByCoord.out.bam \
	--outdir kraken \
	--DBpath /mnt/k2_standard_20201202 \
	--threads `grep -c ^processor /proc/cpuinfo` \
	--prefix metagenomics \
	--kraken /mnt/kraken2/kraken2
#kick out the starsolo mapping and the kraken fastqs
rm kraken/metagenomics_unmapped.fq && rm -r starsolo