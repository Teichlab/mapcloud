#!/bin/bash
set -eo pipefail

#run with one positional argument - the sample ID
SAMPLE=$1

#requires TRUST4 cloned into /mnt/TRUST4

#create output stuff
mkdir $SAMPLE && cd $SAMPLE
mkdir trust && cd trust

#the location of this script is going to be useful for two script calls shortly
DIR=`dirname "$0"`

#reads please
bash $DIR/utils/sample_fastq.sh $SAMPLE

#while unnecessary, I feel safer collating the reads to a single file
cat *R1*.fastq.gz > R1.fastq.gz
cat *R2*.fastq.gz > R2.fastq.gz

#run trust. visium has 16bp CB + 12bp UMI, at least for these samples
#the barcodes are the same for visium v1 and visium v2 btw
/mnt/TRUST4/run-trust4 -t 25 -f /mnt/TRUST4/hg38_bcrtcr.fa --ref /mnt/TRUST4/human_IMGT+C.fa -u R2.fastq.gz --barcode R1.fastq.gz --barcodeRange 0 15 + --UMI R1.fastq.gz --umiRange 16 27 + --barcodeWhitelist ~/cellranger/spaceranger-1.1.0/lib/python/cellranger/barcodes/visium-v1.txt

#delete fastqs as they've served their purposes
rm *.fastq.gz

#post-processing, aiming for a 10x like formatting at the end
#normally you'd use trust's barcoderep-expand.py into trust-barcoderep-to-10x.pl
#however, this has two issues:
#1. it only expands the extra contigs for one of the chains for whatever reason
#2. it keeps the primary two chains paired
#meanwhile we want a situation where we have a fake barcode with one chain in it
#and have all of the identified chains present in the output
#as such, a slight modification was made to the python script to accomplish this
python3 $DIR/barcoderep-allchains.py -b TRUST_R2_barcode_report.tsv > TRUST_R2_barcode_report_allchains.tsv
/mnt/TRUST4/scripts/trust-barcoderep-to-10X.pl TRUST_R2_barcode_report_allchains.tsv 10x-allchains