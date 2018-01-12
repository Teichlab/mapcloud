#!/bin/bash
set -e

samtools view -b $1 | /home/ubuntu/biobambam2/bin/bamcollate2 collate=1 reset=1 resetaux=0 auxfilter=RG,BC,QT | samtools fastq -1 $1\_R1_001.fastq.gz -2 $1\_R2_001.fastq.gz --i1 $1\_I1_001.fastq.gz -n -i --index-format i8 -
