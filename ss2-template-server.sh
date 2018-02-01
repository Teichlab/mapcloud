#!/bin/bash
set -e

#set your remote path to transfer the results to here
REMOTEDIR=kp9@farm3-login.internal.sanger.ac.uk:/lustre/scratch117/cellgen/team205/10x-runs
#set your reference here: GRCh37 or GRCh38
REFERENCE=GRCh37

#loop over run_lane combinations, e.g. 24013_1
for RUNLANE in 
do
	#call the pipeline proper
	#you can swap between GRCh37 and GRCh38, needs to be the first argument
	bash /mnt/mapcloud/scripts/ss2/star-htseq-wrapper.sh $REFERENCE $RUNLANE
	
	#this creates a $RUNLANE output folder. copy it over! specify where below
	sshpass -f ~/.sshpass rsync -Pr $RUNLANE $REMOTEDIR
	
	#and now that we copied over the results, time to burn the input/output to the ground and start anew
	rm -r $RUNLANE
done
