#!/bin/bash
set -e

#turn the CRAM into FASTQ, using a collated BAM intermediate
samtools collate $1 collbams/$1
samtools fastq -1 fastq/$1-1.fastq.gz -2 fastq/$1-2.fastq.gz collbams/$1.bam
