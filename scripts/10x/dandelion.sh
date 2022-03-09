#!/bin/bash
set -eo pipefail

#run with the sample name and location of dandelion container as positional arguments
SAMPLE=$1
DANDELION=$2

#go into the sample's output folder, so dandelion can see the outs, and run the pipeline
cd $SAMPLE
singularity run -B $PWD $DANDELION dandelion-preprocess --chain TR --file_prefix all --keep_trailing_hyphen_number --filter_to_high_confidence
#rename outs to dandelion to signify this is different
mv outs dandelion