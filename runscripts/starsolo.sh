#!/bin/bash
set -eo pipefail

#set reference here. must be present as a star index folder in ~/STAR-2.7.3a
REFERENCE=GRCh38-3.0.0
#choose 10x chemistry. v2, v3 for 3'; or v2-5 if 5'
CHEMISTRY=v2-5

for SAMPLE in 
do
	#run the starsolo wrapper
	bash /mnt/mapcloud/scripts/starsolo/starsolo-wrapper.sh $SAMPLE $REFERENCE $CHEMISTRY
	
	#dump the thing on irods and reset
	bash /mnt/mapcloud/scripts/irods-upload.sh 10X $SAMPLE
	rm -r $SAMPLE
done

#iRODS diagnostic - how much space is in use?
iquest "%s" "select sum(DATA_SIZE) where COLL_NAME like '/archive/HCA/%'" | while read B dummy;do echo "CURRENT /archive/HCA USE: "$((B/(2**40)))" TiB";done