#!/bin/bash
set -e

#set your remote path to transfer the results to here
REMOTEDIR=$SSHNAME@farm3-login.internal.sanger.ac.uk:/lustre/scratch117/cellgen/team205
#set your file type here: cram or fastq
FILETYPE=cram
#set your reference here: GRCh37 or GRCh38
REFERENCE=GRCh37
#set your VCF creation here (ignore if you're not doing the genotyping): 
#0 to use all the protein coding SNPs,
#a nonzero integer to use the SNPs for the top N expressed genes for your sample,
#a file path to use the SNPs from that file
VCF=1000

#assert that $SSHNAME needs to exist
if [ -z ${SSHNAME+x} ]
then
	echo "The \$SSHNAME variable is not set, cannot transfer files to server. Check https://github.com/Teichlab/mapcloud for setup details."
	exit 1
fi

#loop over all the samples you want mapped
for SAMPLE in 
do
	#run the cellranger-wrapper.sh pipeline, providing cram/fastq input as argument 1
	#and GRCh37/GRCh38 reference choice as argument 2. arguments 3+4 form an imeta query pair
	#you can optionally provide a second one to refine the search if needed in 5+6
	#argument 4 will be used for naming the cellranger output files
	bash /mnt/mapcloud/scripts/10x/cellranger-wrapper.sh $FILETYPE $REFERENCE sample $SAMPLE
	
	#EmptyDrops droplet calling, FCA style. comment out if unwanted
	#Rscript /mnt/mapcloud/scripts/10x/emptydrops.R $SAMPLE $REFERENCE
	
	#genotyping script. comment out if unwanted
	#ONLY WORKS WITH A GRCh37 REFERENCE
	#bash /mnt/mapcloud/scripts/10x/genotyping.sh $SAMPLE $VCF
	
	#and now that we ran cellranger, we can toss it over to the farm for safekeeping. specify where below
	sshpass -f ~/.sshpass rsync -Pr $SAMPLE $REMOTEDIR
	
	#and now that we copied over the results, time to burn the input/output to the ground and start anew
	rm -r fastq
	rm -r $SAMPLE
done
