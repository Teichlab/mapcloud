#!/bin/bash
set -e

#run by passing run_lane into the thing, like what dump_irods.sh expects
#but also pass GRCh37 or GRCh38 as the first argument to identify the reference

#so since we have what dump_irods.sh expects, let's use it
bash /mnt/mapcloud/scripts/ss2/dump_irods.sh $2

#at this stage, we have all the crams, corresponding .imeta files
#and a sampleInfo.txt file with conventionally named samples

#the sampleInfo file is an actual piece of output
#but it might be occasionally littered with yhuman garbanzo. lose that
cd $2
mkdir outs
grep -v 'yhuman' sampleInfo.txt > outs/sampleInfo.txt
rm sampleInfo.txt

#so, as mentioned, yhuman is not the best. just torch it if spotted
rm -f *yhuman*

#now we can loop over the CRAMs and do the mapping and stuff
#start off by making them be FASTQ
mkdir fastq
parallel bash /mnt/mapcloud/scripts/ss2/cramfastq.sh ::: *.cram

mkdir STAR
for FID in *.cram
do
	#set up STAR output directory and run STAR itself
	mkdir STAR/$FID && cd STAR/$FID
	~/cellranger/cellranger-2.0.2/STAR/2.5.1b/STAR --runMode alignReads \
	--runThreadN 28 \
	--genomeDir ~/cellranger/$1/star \
	--readFilesIn ../../fastq/$FID-1.fastq.gz ../../fastq/$FID-2.fastq.gz \
	--outSAMtype BAM SortedByCoordinate \
	--outSAMattributes NH HI AS NM MD \
	--outReadsUnmapped Fastx \
	--outFilterMultimapNmax 20 \
	--outFilterType BySJout \
	--outFilterIntronMotifs RemoveNoncanonicalUnannotated \
	--readFilesCommand zcat
	#and now back up to the original level!
	cd ../..
done

#and now that we've mapped what we needed to map, we can kick the FASTQs out
rm -r fastq
#and now parallelise HTSeq as it sees the need to import the index from scratch each time
#so this gives quite a performance boost
parallel bash /mnt/mapcloud/scripts/ss2/htseq.sh ::: *.cram

#the HTSeq reports feature some QC-like lines at the end, let's not use those when making the count matrix
#but first, let's obtain the gene IDs (from the first HTSeq results file) and prepare the file
grep -v '_' STAR/`ls *.cram | head -n 1`/HTSeq_count.txt | cut -f 1 > outs/counthold.txt
#and also prepare the thing that will hold cell names, I guess
printf 'ENSG' > outs/cellhold.txt

#now that this is prepared, we can commence the looping festivities!
for FID in *.cram
do
	#slap on the counts
	paste outs/counthold.txt <(grep -v '_' STAR/$FID/HTSeq_count.txt | cut -f 2) > outs/counthold2.txt
	mv outs/counthold2.txt outs/counthold.txt
	#slap on the cell ID
	cellname="$(sed 's/\..*//g' <(echo $FID))"
	printf "\t%s" $cellname >> outs/cellhold.txt
	#fish out the STAR unique mapping percentage and save it too, I guess
	mapprct="$(grep 'Uniquely mapped reads %' STAR/$FID/Log.final.out | sed 's/.*|//g' | tr -d '[:space:]')"
	printf "%s\t%s\n" $cellname $mapprct >> outs/uniqueMappedPercent.txt
done

#convert the cell names to the actual names, courtesy of a helper hard-wired python script
python3 /mnt/mapcloud/scripts/ss2/celltranslate.py

#and now make the actual output!
cat outs/cellhold.txt <(echo) outs/counthold.txt > outs/countMatrix.txt
cat outs/cellnamehold.txt <(echo) outs/counthold.txt > outs/countMatrixNames.txt

#and now to clean up the assorted holder debris that got made along the way
#along with the crams, and move imeta to its own folder
rm outs/*hold.txt
rm *.cram
mkdir imeta_dump
mv *.imeta imeta_dump
