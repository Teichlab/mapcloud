#!/bin/bash
set -e

#set your reference here: GRCh38 or mm10
REFERENCE=GRCh38
#set your barcodes file here
BARCODES=barcodes.txt

for SAMPLE in 
do
	#the st_pipeline wants STAR 2.5.4b on the path, the wrapper does this
	#but just in case it crashes out, add a clause to revert the path to normalcy
	PATHHOLD=$PATH
	bash /mnt/mapcloud/scripts/spatial/spatial-wrapper.sh $BARCODES $REFERENCE sample $SAMPLE || ( PATH=$PATHHOLD && exit 1 )
done
