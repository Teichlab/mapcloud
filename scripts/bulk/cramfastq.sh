#!/bin/bash
set -e

#turn the CRAM into FASTQ, using a collated BAM intermediate
samtools collate -Ou $1 $1.holder | samtools fastq -1 $1-1.fastq.gz -2 $1-2.fastq.gz -
