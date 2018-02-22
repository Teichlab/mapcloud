#!/bin/bash
set -e

#sometimes BAMs are empty. this will make htseq error out
#only launch htseq if the BAM is nonempty in that case. leave a blank file otherwise
touch STAR/$1/HTSeq_count.txt
if [[ `samtools view STAR/$1/Aligned.sortedByCoord.out.bam | head -n 1 | wc -l` > 0 ]]
then
	htseq-count -f bam -r pos -s no -t exon -i gene_id STAR/$1/Aligned.sortedByCoord.out.bam \
	~/cellranger/$2/genes/genes.gtf > STAR/$1/HTSeq_count.txt
fi
