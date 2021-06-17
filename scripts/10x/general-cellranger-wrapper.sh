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
#an optional fifth argument specifying the TR/GD/IG chain for VDJ use only
#(specifying GD triggers dandelion postprocessing)

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
		#but is it GD?
		if [[ $5 == "GD" ]]
		then
			#cellranger won't recognise GD, need to set the chain to TR
			CHAIN=TR
		else
			CHAIN=$5
		fi
		#are we using primers?
		PRIMERS=""
		if [ -f "primers.txt" ]
		then
			PRIMERS="--inner-enrichment-primers=primers.txt"
		fi
		~/cellranger/cellranger-$VERSION/cellranger $COMMAND --id=$SAMPLE --fastqs=fastq --reference=/home/ubuntu/cellranger/$REFERENCE --chain $CHAIN $PRIMERS
	fi
fi

#wipe out temporary files as those cause nothing but trouble
ls -d $SAMPLE/* | grep -v 'outs' | xargs rm -r

#optional dandelion postprocessing of the GD
if [[ $# -eq 5 ]]
then
	if [[ $5 == "GD" ]]
	then
		#identify high confidence contigs, cell status be damned
		cut -f3-4 -d, $SAMPLE/outs/all_contig_annotations.csv | grep True | cut -f1 -d, > $SAMPLE/outs/hiconf.txt
		head -n 1 $SAMPLE/outs/all_contig_annotations.csv > $SAMPLE/outs/hiconf_contig_annotations.csv
		grep -f $SAMPLE/outs/hiconf.txt $SAMPLE/outs/all_contig_annotations.csv >> $SAMPLE/outs/hiconf_contig_annotations.csv 
		grep -A1 --no-group-separator -f $SAMPLE/outs/hiconf.txt $SAMPLE/outs/all_contig.fasta > $SAMPLE/outs/hiconf_contig.fasta
		#at this point we have hiconf files ready for dandelion ingestion and don't need hiconf.txt anymore
		rm $SAMPLE/outs/hiconf.txt && cd $SAMPLE
		singularity run -B $PWD ~/sc-dandelion_latest.sif dandelion-preprocess --chain TR --file_prefix hiconf --keep_trailing_hyphen_number
	fi
fi

#leave the fastq input and complete output available to be dealt with in whatever manner