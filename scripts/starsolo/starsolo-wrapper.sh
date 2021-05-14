#!/bin/bash
set -e

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

#acquire fastqs via standard 10x proceedings
mkdir fastq && cd fastq
printf '#!/bin/bash\nset -e\n\n----\n' > imeta.sh
imeta qu -z seq -d sample = $SAMPLE and type = cram and target = 1 >> imeta.sh
#this used to be -iget -K, but there have been some bizarre ghost errors early Q2 2021
#yielding irreproducible borked reads at a rare rate without scripts breaking
#so as a precaution, do a manual md5sum check instead, just in case it's a ghost irods error
#(can skip the -K for run time as will check md5sums later)
sed ':a;N;$!ba;s/----\ncollection:/iget/g' -i imeta.sh
sed ':a;N;$!ba;s/\ndataObj: /\//g' -i imeta.sh
#clean up the imeta by kicking out *#888.crams and yhumans
grep -v "#888.cram\|yhuman" imeta.sh > temp.sh && mv temp.sh imeta.sh
bash imeta.sh
for IRPATH in `grep "iget" imeta.sh | sed "s/iget //"`
do
	FID=`basename $IRPATH`
	MD5LOCAL=`md5sum $FID| cut -f -1 -d " "`
	MD5IRODS=`imeta ls -d $IRPATH md5 | grep 'value:' | sed 's/value: //'`
	if [ $MD5LOCAL != $MD5IRODS ]
	then
		echo "md5sum conflict encountered for $FID"
		exit 1
	fi
done
rm imeta.sh
#convert to FASTQ and collapse
parallel bash /mnt/mapcloud/scripts/10x/utils/cramfastq.sh ::: *.cram
cat *R1_001.fastq.gz > R1.fastq.gz
cat *R2_001.fastq.gz > R2.fastq.gz
rm *cram* && cd ..

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