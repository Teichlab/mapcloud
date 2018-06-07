#!/bin/bash
set -e

#set your reference here: GRCh38 or mm10
REFERENCE=GRCh38

#loop over all the samples you want mapped
for SAMPLE in 
do
	#run the star-htseq-wrapper.sh pipeline, providing GRCh38/mm10 reference choice 
	#as argument 1. arguments 2+3 form an imeta query pair
	#you can optionally provide a second one to refine the search if needed in 4+5
	#argument 3 will be used for naming the output folder
	bash /mnt/mapcloud/scripts/bulk/star-htseq-wrapper.sh $REFERENCE sample $SAMPLE
		
	#and now that we ran STAR+HTSeq, we can toss it over to irods for storage
	bash /mnt/mapcloud/scripts/irods-upload.sh BULK $SAMPLE
	
	#and now that we copied over the results, time to burn the output to the ground and start anew
	rm -r $SAMPLE
done
