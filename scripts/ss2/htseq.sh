#!/bin/bash
set -e

htseq-count -f bam -r pos -s no -t exon -i gene_id STAR/$1/Aligned.sortedByCoord.out.bam \
/home/ubuntu/cellranger/GRCh37/genes/genes.gtf > STAR/$1/HTSeq_count.txt
