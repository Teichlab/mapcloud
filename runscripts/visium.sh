#!/bin/bash
set -e

#set your file type here: cram or fastq
FILETYPE=cram
#set your reference here: GRCh38 or mm10
REFERENCE=GRCh38-3.0.0
#set the metadata file here
META=meta.csv

#the metadata file should be a csv:
#	- first column sample ID
#	- second column image name (placed in the directory at the same level as the runscript)
#	- third column slide serial number
#	- fourth column spot location (A1/B1/C1/D1)
#	(- fifth column - optionally manual alignment if too few keypoints detected error happens otherwise)


#loop over all the samples you want mapped
for SAMPLE in 
do
	#run the cellranger-wrapper.sh pipeline, providing cram/fastq input as argument 1
	#and GRCh38/mm10 reference choice as argument 2. argument 3 is the version to use.
	#arguments 4+5 form an imeta query pair, you can optionally provide a second one to 
	#refine the search if needed in 6+7. argument 4 will be used for naming the 
	#cellranger output files
	bash /mnt/mapcloud/scripts/10x/spaceranger-wrapper.sh $FILETYPE $REFERENCE $META sample $SAMPLE
	
	#and now that we ran cellranger, we can toss it over to irods for storage
	bash /mnt/mapcloud/scripts/irods-upload.sh VISIUM $SAMPLE
	
	#and now that we copied over the results, time to burn the input/output to the ground and start anew
	rm -r fastq
	rm -r $SAMPLE
done

#iRODS diagnostic - how much space is in use?
iquest "%s" "select sum(DATA_SIZE) where COLL_NAME like '/archive/HCA/%'" | while read B dummy;do echo "CURRENT /archive/HCA USE: "$((B/(2**40)))" TiB";done

