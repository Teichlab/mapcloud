#!/bin/bash
set -e

#set your file type here: cram or fastq
FILETYPE=cram
#set your reference here: GRCh38 or mm10
REFERENCE=GRCh38
#set the cellranger version here: 2.0.2 for 3'v2, 2.1.1 for 5', 3.0.2 for 3'v3
VERSION=2.0.2
#set your VCF creation here (ignore if you're not doing the genotyping): 
#0 to use all the protein coding SNPs,
#a nonzero integer to use the SNPs for the top N expressed genes for your sample,
#a file path to use the SNPs from that file
VCF=1000

#loop over all the samples you want mapped
for SAMPLE in 
do
	#run the cellranger-wrapper.sh pipeline, providing cram/fastq input as argument 1
	#and GRCh38/mm10 reference choice as argument 2. argument 3 is the version to use.
	#arguments 4+5 form an imeta query pair, you can optionally provide a second one to 
	#refine the search if needed in 6+7. argument 4 will be used for naming the 
	#cellranger output files
	bash /mnt/mapcloud/scripts/10x/cellranger-wrapper.sh $FILETYPE $REFERENCE $VERSION sample $SAMPLE
	
	#genotyping script. comment out if unwanted
	#ONLY WORKS WITH A GRCh38 REFERENCE
	#bash /mnt/mapcloud/scripts/10x/genotyping.sh $SAMPLE $VCF
	
	#and now that we ran cellranger, we can toss it over to irods for storage
	bash /mnt/mapcloud/scripts/irods-upload.sh 10X $SAMPLE
	
	#and now that we copied over the results, time to burn the input/output to the ground and start anew
	rm -r fastq
	rm -r $SAMPLE
done

#iRODS diagnostic - how much space is in use?
iquest "%s" "select sum(DATA_SIZE) where COLL_NAME like '/archive/HCA/%'" | while read B dummy;do echo "CURRENT /archive/HCA USE: "$((B/(2**40)))" TiB";done
