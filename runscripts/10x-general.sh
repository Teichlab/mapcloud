#!/bin/bash
set -e

#set the cellranger version here, must be present as ~/cellranger/cellranger-$VERSION
VERSION=4.0.0
#set the cellranger command here - count for GEX, vdj for VDJ
COMMAND=count
#set your reference here, must be present as ~/cellranger/$REFERENCE
REFERENCE=GRCh38-3.0.0

#loop over all the samples you want mapped
for SAMPLE in 
do
	#run the general-cellranger-wrapper.sh pipeline
	#optionally provide TR/IG as the chain to force if running VDJ
	bash /mnt/mapcloud/scripts/10x/general-cellranger-wrapper.sh $VERSION $COMMAND $REFERENCE $SAMPLE
	
	#and now that we ran cellranger, we can toss it over to irods for storage
	#the destination folder depends on whether this is VDJ or not
	if [[ $COMMAND == "count" ]]
	then
		bash /mnt/mapcloud/scripts/irods-upload.sh 10X $SAMPLE
	else
		bash /mnt/mapcloud/scripts/irods-upload.sh 10X-VDJ $SAMPLE
	fi
	
	#and now that we copied over the results, time to burn the input/output to the ground and start anew
	rm -r fastq
	rm -r $SAMPLE
done

#iRODS diagnostic - how much space is in use?
iquest "%s" "select sum(DATA_SIZE) where COLL_NAME like '/archive/HCA/%'" | while read B dummy;do echo "CURRENT /archive/HCA USE: "$((B/(2**40)))" TiB";done
