#!/bin/bash
set -eo pipefail

#run with the sample name and location of dandelion container as positional arguments
SAMPLE=$1
DANDELION=$2

#identify high confidence contigs, cell status be damned
cut -f3-4 -d, $SAMPLE/outs/all_contig_annotations.csv | grep -i True | cut -f1 -d, > $SAMPLE/outs/hiconf.txt
head -n 1 $SAMPLE/outs/all_contig_annotations.csv > $SAMPLE/outs/hiconf_contig_annotations.csv
grep -f $SAMPLE/outs/hiconf.txt $SAMPLE/outs/all_contig_annotations.csv >> $SAMPLE/outs/hiconf_contig_annotations.csv 
grep -A1 --no-group-separator -f $SAMPLE/outs/hiconf.txt $SAMPLE/outs/all_contig.fasta > $SAMPLE/outs/hiconf_contig.fasta
#at this point we have hiconf files ready for dandelion ingestion and don't need hiconf.txt anymore
rm $SAMPLE/outs/hiconf.txt && cd $SAMPLE
singularity run -B $PWD $DANDELION dandelion-preprocess --chain TR --file_prefix hiconf --keep_trailing_hyphen_number
#rename outs to dandelion to signify this is different
mv outs dandelion