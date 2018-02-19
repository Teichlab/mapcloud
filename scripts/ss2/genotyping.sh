#!/bin/bash
set -e

#positional arguments: $1 RUN_LANE, $2 VCF info (0 for whole protein-coding file, int>0 for top N genes, path for custom file)

#enter the results subdirectory
cd $1

#process VCF accordingly - grab appropriate file, or create a new one
#ensure that it lives in snps.vcf at the end of it
if [ -f $2 ]
then
	#custom file. copy it over
	cp $2 snps.vcf
elif [ $2 == 0 ]
then
	#whole protein-coding SNP file. copy it over
	cp ~/vcf-protein-coding/dbsnp_138.hg19.protein_coding.recode.dedupped.vcf snps.vcf
else
	#create our own one. start by finding all protein_coding gene lines in the GTF
	grep -P 'protein_coding\tgene\t' ~/cellranger/GRCh37/genes/genes.gtf > genelines.gtf
	#find the top N expressed genes in our mapped sample
	Rscript /mnt/mapcloud/scripts/ss2/find-topgenes.R $2
	#and now find the gene lines for those genes
	grep -f topgenes.txt genelines.gtf > topgenelines.gtf
	#and now intersect that with the protein coding SNP file to get these genes' SNPs
	bedtools intersect -a ~/vcf-protein-coding/dbsnp_138.hg19.protein_coding.recode.dedupped.vcf -b topgenelines.gtf > snps.vcf
	#cleanup temp files - they all handily have gene in the name
	rm *gene*
fi

#genotyping proper!
mkdir genotyping-singlecell-vcfs && cd genotyping-singlecell-vcfs
#call bcftools mpileup | call on all the BAMs in the level above
parallel bash /mnt/mapcloud/scripts/ss2/mpileup.sh ::: `find .. -type f -name '*.bam'`
#do some editing to make the "sample name" more manageable
for FID in *.vcf
do
	sed "s/\.\.\/STAR\///" -i $FID
	sed "s/\/Aligned.sortedByCoord.out.bam//" -i $FID
	bgzip $FID && bcftools index $FID.gz
done
#merge the individual VCF files
cd ..
bcftools merge -Ov --threads `grep -c ^processor /proc/cpuinfo` -o outs/$1.vcf genotyping-singlecell-vcfs/*.vcf.gz

#snps.vcf has done its job
rm snps.vcf
