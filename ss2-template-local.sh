#!/bin/bash
set -e

#loop over run_lane combinations, e.g. 24013_1
for RUNLANE in 
do
	#call the pipeline proper
	#you can swap between GRCh37 and GRCh38, needs to be the first argument
	bash /mnt/mapcloud/scripts/ss2/star-htseq-wrapper.sh GRCh37 $RUNLANE
	
	#the actual output lives in $RUNLANE/outs
	mv $RUNLANE/outs holder-output
		
	#and now that we copied over the results, time to burn the input/output to the ground and start anew
	rm -r $RUNLANE
	
	#...and safekeep the results
	mv holder-output $RUNLANE
done

