#!/bin/bash
set -e

#set your remote path to transfer the results to here
REMOTEDIR=$SSHNAME@farm3-login.internal.sanger.ac.uk:/lustre/scratch117/cellgen/team205
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

#assert that the password file is there
if [ ! -f ~/.sshpass ]
then
	echo "~/.sshpass does not exist. Create it and put your SSH password there."
	exit 1
fi

#loop over run_lane combinations, e.g. 24013_1
for RUNLANE in 
do
	#call the pipeline proper
	#you can swap between GRCh37 and GRCh38, needs to be the first argument
	bash /mnt/mapcloud/scripts/ss2/star-htseq-wrapper.sh $REFERENCE $RUNLANE
	
	#genotyping script. comment out if unwanted
	#ONLY WORKS WITH A GRCh37 REFERENCE
	#bash /mnt/mapcloud/scripts/ss2/genotyping.sh $RUNLANE $VCF
	
	#this creates a $RUNLANE output folder. copy it over! specify where below
	sshpass -f ~/.sshpass rsync -Pr $RUNLANE $REMOTEDIR
	
	#and now that we copied over the results, time to burn the input/output to the ground and start anew
	rm -r $RUNLANE
done
