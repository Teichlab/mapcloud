#!/bin/bash
set -e

#set your file type here: cram or fastq
FILETYPE=cram
#set your reference here: GRCh37 or GRCh38
REFERENCE=GRCh37
#set your VCF creation here (ignore if you're not doing the genotyping): 
#0 to use all the protein coding SNPs,
#a nonzero integer to use the SNPs for the top N expressed genes for your sample,
#a file path to use the SNPs from that file
VCF=1000

#loop over all the samples you want mapped
for SAMPLE in 
do
	#run the cellranger-wrapper.sh pipeline, providing cram/fastq input as argument 1
	#and GRCh37/GRCh38 reference choice as argument 2. arguments 3+4 form an imeta query pair
	#you can optionally provide a second one to refine the search if needed in 5+6
	#argument 4 will be used for naming the cellranger output files
	bash /mnt/mapcloud/scripts/10x/cellranger-wrapper.sh $FILETYPE $REFERENCE sample $SAMPLE
	
	#EmptyDrops droplet calling, FCA style. comment out if unwanted
	#Rscript /mnt/mapcloud/scripts/10x/emptydrops.R $SAMPLE $REFERENCE
	
	#genotyping script. comment out if unwanted
	#ONLY WORKS WITH A GRCh37 REFERENCE
	#bash /mnt/mapcloud/scripts/10x/genotyping.sh $SAMPLE $VCF
	
	#save small selection of output locally instead of dumping the whole thing on the farm
	#save EmptyDrops/genotyping results too if those got made
	mkdir output-holder
	mv $SAMPLE/outs/filtered_gene_bc_matrices/$REFERENCE output-holder/cellranger
	mv $SAMPLE/outs/web_summary.html output-holder/cellranger
	if [ -d $SAMPLE/outs/final-count-matrix ]
	then
		mv $SAMPLE/outs/final-count-matrix output-holder
	fi
	if [ -f $SAMPLE/outs/$SAMPLE.vcf ]
	then
		mv $SAMPLE/outs/$SAMPLE.vcf output-holder
	fi
	
	#and now that we copied over the results, time to burn the input/output to the ground and start anew
	rm -r fastq
	rm -r $SAMPLE
	
	#...and reposition the output holder as the sample's final output
	mv output-holder $SAMPLE
done
