#!/bin/bash
set -e

#set your file type here: cram or fastq
FILETYPE=cram
#set your reference here: GRCh38-VDJ
REFERENCE=GRCh38-VDJ
#set your VCF creation here (ignore if you're not doing the genotyping): 
#0 to use all the protein coding SNPs,
#a nonzero integer to use the SNPs for the top N expressed genes for your sample,
#a file path to use the SNPs from that file
VCF=1000

#loop over all the samples you want mapped
for SAMPLE in 
do
	#run the cellranger-vdj-wrapper.sh pipeline, providing cram/fastq input as argument 1
	#and GRCh38-VDJ reference choice as argument 2. arguments 3+4 form an imeta query pair
	#you can optionally provide a second one to refine the search if needed in 5+6
	#argument 4 will be used for naming the cellranger output files
	bash /mnt/mapcloud/scripts/10x/cellranger-vdj-wrapper.sh $FILETYPE $REFERENCE sample $SAMPLE
	
	#and now that we ran cellranger, we can toss it over to irods for storage
	bash /mnt/mapcloud/scripts/irods-upload.sh 10X-VDJ $SAMPLE
	
	#and now that we copied over the results, time to burn the input/output to the ground and start anew
	rm -r fastq
	rm -r $SAMPLE
done

#iRODS diagnostic - how much space is in use?
iquest "%s" "select sum(DATA_SIZE) where COLL_NAME like '/archive/HCA/%'" | while read B dummy;do echo "CURRENT /archive/HCA USE: "$((B/(2**40)))" TiB";done
