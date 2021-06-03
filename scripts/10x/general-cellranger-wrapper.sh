#!/bin/bash
set -e

#a general wrapper capable of calling both count and vdj

#run with positional arguments like so
#cellranger version to use, located at ~/cellranger/cellranger-$VERSION
VERSION=$1
#the command to call - count or vdj
COMMAND=$2
#the reference to use, located at ~/cellranger/$REFERENCE
REFERENCE=$3
#the sample ID to process
SAMPLE=$4
#an optional fifth argument specifying the TR/IG chain for VDJ use only

#start by getting the FASTQs
bash /mnt/mapcloud/scripts/10x/utils/sample_fastq.sh $SAMPLE
mkdir fastq && mv *.fastq.gz fastq

#now run cellranger proper
#the syntax for passing the reference differs between count and vdj
if [[ $COMMAND == "count" ]]
then
	#count call
	~/cellranger/cellranger-$VERSION/cellranger $COMMAND --id=$SAMPLE --fastqs=fastq --transcriptome=/home/ubuntu/cellranger/$REFERENCE
else
	if [[ $# -eq 4 ]]
	then
		#vdj call
		~/cellranger/cellranger-$VERSION/cellranger $COMMAND --id=$SAMPLE --fastqs=fastq --reference=/home/ubuntu/cellranger/$REFERENCE
	fi
	if [[ $# -eq 5 ]]
	then
		#vdj call with the chain passed
		~/cellranger/cellranger-$VERSION/cellranger $COMMAND --id=$SAMPLE --fastqs=fastq --reference=/home/ubuntu/cellranger/$REFERENCE --chain $5
	fi
fi

#wipe out temporary files as those cause nothing but trouble
ls -d $SAMPLE/* | grep -v 'outs' | xargs rm -r

#leave the fastq input and complete output available to be dealt with in whatever manner