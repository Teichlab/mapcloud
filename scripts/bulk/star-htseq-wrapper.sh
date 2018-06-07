#!/bin/bash
set -e

#three (or five) positional arguments on launch, 2 onwards for use on the imeta call
#	* $1 - GRCh38 or mm10, which reference to use
#	* $2 - sample or sample_id (or id_run or whatever), which field are we checking
#	* $3 - the value that we want to fish out
#	* optionally, $4 and $5 for another pair

#create the output folder and go live there
mkdir $3 && cd $3

#ok, time to get the imeta ball rolling
#prepare the foundation for the imeta dump becoming a shell script with igets
#and then obtain the imeta using the value pair(s) provided on input
printf '#!/bin/bash\nset -e\n\n----\n' > imeta.sh
if [[ $# -eq 3 ]]
then
	imeta qu -z seq -d $2 = $3 and type = cram and target = 1 >> imeta.sh
fi
if [[ $# -eq 5 ]]
then
	imeta qu -z seq -d $2 = $3 and $4 = $5 and type = cram and target = 1 >> imeta.sh
fi

#a couple quick substitutions to actually turn the imeta dump into a bunch of igets
#turn all '----\ncollection:' into 'iget'
#turn all '\ndataObj: ' into '/'
sed ':a;N;$!ba;s/----\ncollection:/iget -K/g' -i imeta.sh
sed ':a;N;$!ba;s/\ndataObj: /\//g' -i imeta.sh

#time to get downloading
bash imeta.sh

#and now that it's downloaded, we can dispose of it
rm imeta.sh

#kick out the unmapped as we don't care; -f so that this doesn't break if this file doesn't exist
rm -f *#888.cram

#process the CRAMs into FASTQs. this is smartseq2-style, not 10x style
parallel bash /mnt/mapcloud/scripts/bulk/cramfastq.sh ::: *.cram

#collate all the read 1's and read 2's into single files
#no fear as to .gz, as the format explicitly supports collating compressed files
cat *-1.fastq.gz > R1.fastq.gz
cat *-2.fastq.gz > R2.fastq.gz

#run STAR, SS2-style, but just once.
mkdir STAR && cd STAR
~/cellranger/cellranger-2.0.2/STAR/2.5.1b/STAR --runMode alignReads \
	--runThreadN `grep -c ^processor /proc/cpuinfo` \
	--genomeDir ~/cellranger/$1/star \
	--readFilesIn ../R1.fastq.gz ../R2.fastq.gz \
	--outSAMtype BAM SortedByCoordinate \
	--outSAMattributes NH HI AS NM MD \
	--outReadsUnmapped Fastx \
	--outFilterMultimapNmax 20 \
	--outFilterType BySJout \
	--outFilterIntronMotifs RemoveNoncanonicalUnannotated \
	--readFilesCommand zcat

#and also just one very basic call of HTSeq
htseq-count -f bam -r pos -s no -t exon -i gene_id Aligned.sortedByCoord.out.bam \
	~/cellranger/$1/genes/genes.gtf > HTSeq.txt

#ok, time to prep the output (just HTSeq counts) and clean stuff up
cd .. && rm *.*
mkdir outs && grep -v '_' STAR/HTSeq.txt > outs/counts.txt
cp STAR/Log.final.out outs/Star.log.final.out
