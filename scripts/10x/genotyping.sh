#!/bin/bash
set -e

#positional arguments: $1 sample name, $2 VCF info (0 for whole protein-coding file, int>0 for top N genes, path for custom file)

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
	Rscript /mnt/mapcloud/scripts/10x/find-topgenes.R $1 $2
	#and now find the gene lines for those genes
	grep -f topgenes.txt genelines.gtf > topgenelines.gtf
	#and now intersect that with the protein coding SNP file to get these genes' SNPs
	bedtools intersect -a ~/vcf-protein-coding/dbsnp_138.hg19.protein_coding.recode.dedupped.vcf -b topgenelines.gtf > snps.vcf
	#cleanup temp files - they all handily have gene in the name
	rm *gene*
fi

#genotyping proper!
#start by demultiplexing
mkdir out
samtools view -@ 4 -h $1/outs/possorted_genome_bam.bam | perl -nle 'use strict; use autodie; our %h; BEGIN{open(my$fh,q(<),shift@ARGV); my$od=shift@ARGV; $od//=q(); while(<$fh>){chomp; open(my$f2,"| samtools view -u - |bamstreamingmarkduplicates level=0 tag=UB | samtools view -b - > $od/$_.bam");$h{$_}=$f2; }close $fh}  if(/^@/){foreach my$fh (values %h){print {$fh} $_ }}elsif(m{\tCB:Z:(\S+)\b}){ my$fh=$h{$1}||next; print {$fh} $_;} END{close $_ foreach values %h; warn "closed BAMs\n"}' $1/outs/filtered_gene_bc_matrices/GRCh37/barcodes.tsv out
cd out
#now call bcftools mpileup | call
parallel bash /mnt/mapcloud/scripts/10x/mpileup.sh ::: *.bam
#do some editing to include the actual cell barcode in the VCF, otherwise it just stores the sample ID somehow
for FID in *.vcf
do
	sed "s/$1$/$1\_`basename $FID .vcf`/" -i $FID
	bgzip $FID && bcftools index $FID.gz
done
#merge the individual VCF files and copy them over to the output dump just in case
cd ..
bcftools merge -Ov --threads `grep -c ^processor /proc/cpuinfo` -o $1/outs/$1.vcf out/*.vcf.gz
rm out/*.bam*
mv out $1/genotyping-singlecell-vcfs

#snps.vcf has done its job
rm snps.vcf
