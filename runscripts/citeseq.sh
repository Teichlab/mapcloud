#!/bin/bash
set -e

#set reference here. must be present as a cellranger index folder in ~/cellranger
REFERENCE=GRCh38

#SAMPLE here is a combination of a single CITE's GEX and CITE joined with a + sign
#(e.g. Imm_FLNG8965969+Imm_FLNG8965970)
for SAMPLE in 
do
	#run the starsolo wrapper
	bash /mnt/mapcloud/scripts/citeseq/citeseq-cellranger-wrapper.sh $SAMPLE $REFERENCE
	
	#dump the thing on irods and reset
	bash /mnt/mapcloud/scripts/irods-upload.sh 10X $SAMPLE
	rm -r $SAMPLE
	rm -r fastq_gex
	rm -r fastq_cite
done

#iRODS diagnostic - how much space is in use?
iquest "%s" "select sum(DATA_SIZE) where COLL_NAME like '/archive/HCA/%'" | while read B dummy;do echo "CURRENT /archive/HCA USE: "$((B/(2**40)))" TiB";done