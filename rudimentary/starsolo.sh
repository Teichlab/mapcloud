#!/bin/bash
set -e

REFERENCE=GRCh38
#v2, v3 for 3'; or v2-5 if 5'
CHEMISTRY=v2


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

for SAMPLE in 
do
	#set up nested folder architecture to match /archive/HCA irods design
	mkdir $SAMPLE && cd $SAMPLE && mkdir starsolo && cd starsolo
	
	#acquire fastqs via standard 10x proceedings
	mkdir fastq && cd fastq
	printf '#!/bin/bash\nset -e\n\n----\n' > imeta.sh
	imeta qu -z seq -d sample = $SAMPLE and type = cram and target = 1 >> imeta.sh
	sed ':a;N;$!ba;s/----\ncollection:/iget -K/g' -i imeta.sh
	sed ':a;N;$!ba;s/\ndataObj: /\//g' -i imeta.sh
	bash imeta.sh && rm imeta.sh && rm -f *#888.cram
	parallel bash /mnt/mapcloud/scripts/10x/cramfastq.sh ::: *.cram
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
	python3 /mnt/mapcloud/rudimentary/postprocess.py
	
	#reshuffle some files and run soupX
	mv logs counts && mkdir logs && mv counts/Barcodes.stats logs && mv Log.out logs && mv Log.progress.out logs && rm -r fastq
	Rscript /mnt/mapcloud/rudimentary/soupx.R
	#for whatever reason this creates this empty, unnecessary file. kick it
	rm -f Rplots.pdf
	
	#dump the thing on irods and reset
	cd ../.. && bash /mnt/mapcloud/scripts/irods-upload.sh 10X $SAMPLE
	rm -r $SAMPLE
done


#set up as follows
#(needed to regenerate star index as new star versions that support solo
#have some sort of slightly different index thing going on)
~/STAR-2.7.3a/bin/Linux_x86_64/STAR --runMode genomeGenerate --runThreadN 26 --genomeDir ~/STAR-2.7.3a/GRCh38 --genomeFastaFiles ~/cellranger/GRCh38/fasta/genome.fa --sjdbGTFfile ~/cellranger/GRCh38/genes/genes.gtf
mkdir ~/STAR-2.7.3a/whitelists && cd ~/STAR-2.7.3a/whitelists
cp ~/cellranger/cellranger-3.0.2/cellranger-cs/3.0.2/lib/python/cellranger/barcodes/3M-february-2018.txt.gz v3.txt.gz
gunzip v3.txt.gz
cp ~/cellranger/cellranger-3.0.2/cellranger-cs/3.0.2/lib/python/cellranger/barcodes/737K-august-2016.txt v2.txt
