#!/bin/bash
set -e

#turn the CRAM into FASTQ, using a collated BAM intermediate
samtools collate -Ou $1 $1.holder | samtools fastq -1 fastq/$1-1.fastq.gz -2 fastq/$1-2.fastq.gz -
