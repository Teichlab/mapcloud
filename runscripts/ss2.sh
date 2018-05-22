#!/bin/bash
set -e

#set your reference here: GRCh38 or mm10
REFERENCE=GRCh38
#set your VCF creation here (ignore if you're not doing the genotyping): 
#0 to use all the protein coding SNPs,
#a nonzero integer to use the SNPs for the top N expressed genes for your sample,
#a file path to use the SNPs from that file
VCF=1000

#loop over run_lane combinations, e.g. 24013_1
for RUNLANE in 
do
	#call the pipeline proper
	#you can swap between GRCh38 and mm10, needs to be the first argument
	bash /mnt/mapcloud/scripts/ss2/star-htseq-wrapper.sh $REFERENCE $RUNLANE
	
	#compute TPMs. comment out if unwanted
	Rscript /mnt/mapcloud/scripts/ss2/tpm.R $RUNLANE $REFERENCE
	
	#genotyping script. comment out if unwanted
	#ONLY WORKS WITH A GRCh38 REFERENCE
	#bash /mnt/mapcloud/scripts/ss2/genotyping.sh $RUNLANE $VCF
	
	#this creates a $RUNLANE output folder. copy it over! specify where below
	bash /mnt/mapcloud/scripts/irods-upload.sh SS2 $SAMPLE
	
	#and now that we copied over the results, time to burn the input/output to the ground and start anew
	rm -r $RUNLANE
done
