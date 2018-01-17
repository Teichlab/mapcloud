#!/bin/bash
set -e

#loop over run_lane combinations, e.g. 24013_1
for RUNLANE in 
do
	#call the pipeline proper
	#you can swap between GRCh37 and GRCh38, needs to be the first argument
	bash /mnt/mapcloud/scripts/ss2/star-htseq.sh GRCh37 $RUNLANE
	
	#this creates a $RUNLANE output folder. copy it over! specify where below
	sshpass -f ~/.sshpass rsync -Pr $RUNLANE kp9@farm3-login.internal.sanger.ac.uk:/lustre/scratch117/cellgen/team205/
	
	#and now that we copied over the results, time to burn the input/output to the ground and start anew
	rm -r $RUNLANE
done
