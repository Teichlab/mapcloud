#!/bin/bash
set -e

#discriminate between old and new CRAM files based on the presence of @SQ lines in the CRAM header
if [[ `samtools view -H $1 | grep '@SQ' | wc -l` == 0 ]]
then
	#new file, no alignment, can just go straight for FASTQ
	samtools fastq -1 $1\_R1_001.fastq.gz -2 $1\_R2_001.fastq.gz --i1 $1\_I1_001.fastq.gz -n -i --index-format i8 $1
else
	#old file with alignment, needs bamcollate2 incantation
	samtools view -b $1 | bamcollate2 collate=1 reset=1 resetaux=0 auxfilter=RG,BC,QT | samtools fastq -1 $1\_R1_001.fastq.gz -2 $1\_R2_001.fastq.gz --i1 $1\_I1_001.fastq.gz -n -i --index-format i8 -
fi
