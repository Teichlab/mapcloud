#!/bin/bash
set -e

#helper function for downloading CRAMs and FASTQ'ing them
getfastq() {
	#positional arguments like so
	SAMPLE=$1
	OUTDIR=$2
	#commence the CRAM download
	printf '#!/bin/bash\nset -e\n\n----\n' > imeta.sh
	imeta qu -z seq -d sample = $SAMPLE and type = cram and target = 1 >> imeta.sh
	sed ':a;N;$!ba;s/----\ncollection:/iget -K/g' -i imeta.sh
	sed ':a;N;$!ba;s/\ndataObj: /\//g' -i imeta.sh
	bash imeta.sh
	rm imeta.sh
	rm -f *#888.cram
	#convert to FASTQ
	parallel bash /mnt/mapcloud/scripts/10x/cramfastq.sh ::: *.cram
	#rename to be cellranger compatible
	count=1
	for file in *I1_001.fastq.gz
	do
		mv "${file}" "${file/*.cram/$SAMPLE\_S$count\_L001}"
		file=$(sed 's/I1/R1/g' <<< $file)
		mv "${file}" "${file/*.cram/$SAMPLE\_S$count\_L001}"
		file=$(sed 's/R1/R2/g' <<< $file)
		mv "${file}" "${file/*.cram/$SAMPLE\_S$count\_L001}"
		(( count++ ))
	done
	mkdir $OUTDIR && mv *.fastq.gz $OUTDIR && rm *cram*
}

#positional arguments to whole script like so
#(SAMPLES is GEX-CITE together, e.g. Imm_FLNG8965969-Imm_FLNG8965970)
SAMPLES=$1
REFERENCE=$2
#also requires a 10x-compatible feature reference file called citeseq.csv in the run folder

#split up GEX and CITE on the + in $SAMPLES (GEX first, CITE second), then download stuff
GEX=`echo $SAMPLES | sed "s/-.*//"`
CITE=`echo $SAMPLES | sed "s/.*-//"`
getfastq $GEX fastq_gex
getfastq $CITE fastq_cite

#prepare cellranger's library file
echo "fastqs,sample,library_type" > libraries.csv
echo $(pwd)"/fastq_gex,"$GEX",Gene Expression" >> libraries.csv
echo $(pwd)"/fastq_cite,"$CITE,"Antibody Capture" >> libraries.csv

#run cellranger itself
~/cellranger/cellranger-3.0.2/cellranger count --id=$SAMPLES --libraries=libraries.csv --feature-ref=citeseq.csv --transcriptome=~/cellranger/$REFERENCE

#wipe out temporary files as those cause nothing but trouble
ls -d $SAMPLES/* | grep -v 'outs' | xargs rm -r

#run SoupX on the folder
Rscript /mnt/mapcloud/scripts/citeseq/soupx.R $SAMPLES

#leave the fastq input and complete output available to be dealt with in whatever manner