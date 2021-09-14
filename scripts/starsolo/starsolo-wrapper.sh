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
bash /mnt/mapcloud/scripts/10x/utils/sample_fastq.sh $SAMPLE
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
	--soloOutFileNames logs/ genes.tsv barcodes.tsv matrix.mtx

#starsolo's take on velocyto is a three value column mtx, with spliced, unspliced, ambiguous
#in order to let downstream tools make use of it, split them up into the appropriate matrices
head -n 3 logs/Velocyto/raw/matrix.mtx > logs/Velocyto/raw/spliced.mtx
head -n 3 logs/Velocyto/raw/matrix.mtx > logs/Velocyto/raw/unspliced.mtx
head -n 3 logs/Velocyto/raw/matrix.mtx > logs/Velocyto/raw/ambiguous.mtx
tail -n +4 logs/Velocyto/raw/matrix.mtx | cut -f 1,2,3 -d ' ' >> logs/Velocyto/raw/spliced.mtx
tail -n +4 logs/Velocyto/raw/matrix.mtx | cut -f 1,2,4 -d ' ' >> logs/Velocyto/raw/unspliced.mtx
tail -n +4 logs/Velocyto/raw/matrix.mtx | cut -f 1,2,5 -d ' ' >> logs/Velocyto/raw/ambiguous.mtx

#do various python'y things - cellranger 3 cell filter, scrublet, velocyto object
python3 /mnt/mapcloud/scripts/starsolo/postprocess.py

#reshuffle some files and run soupX
mv logs counts && mkdir logs && mv counts/Barcodes.stats logs && mv Log.out logs && mv Log.progress.out logs && rm -r fastq
Rscript /mnt/mapcloud/scripts/starsolo/soupx.R
#for whatever reason this creates this empty, unnecessary file. kick it
rm -f Rplots.pdf