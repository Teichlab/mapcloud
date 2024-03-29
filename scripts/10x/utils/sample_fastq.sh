#!/bin/bash
set -eo pipefail

#a script for downloading the CRAMs of a given sample, converting them to FASTQ, 
#and making the nomenclature be compatible with cellranger's expectations

#run with one positional argument - the sample ID
SAMPLE="$1"
#if a second argument is encountered, it's the library type
#(just the part after Chromium single cell)
if [[ $# -eq 2 ]]
then
	LIBRARY="Chromium single cell $2"
fi

#acquire FASTQs via standard 10x proceedings
#start by querying iRODS for the sample's CRAMs
echo "----" > holder.sh
#we may or may not need to include a library type
if [[ $# -eq 2 ]]
then
	imeta qu -z seq -d sample = "$SAMPLE" and library_type = "$LIBRARY" and target = 1 >> holder.sh
else
	imeta qu -z seq -d sample = "$SAMPLE" and type = cram and target = 1 >> holder.sh
fi
#this used to be -iget -K, but there have been some bizarre ghost errors early Q2 2021
#yielding irreproducible borked reads at a rare rate without scripts breaking
#so as a precaution, do a manual md5sum check instead, just in case it's a ghost irods error
#(can skip the -K for run time as will check md5sums later)
sed ':a;N;$!ba;s/----\ncollection:/iget/g' -i holder.sh
sed ':a;N;$!ba;s/\ndataObj: /\//g' -i holder.sh
#clean up the imeta by kicking out *#888.crams and yhumans
grep -v "#888.cram\|yhuman" holder.sh > temp.sh && mv temp.sh holder.sh
#make sure that it's all crams though
grep ".cram" holder.sh > temp.sh && mv temp.sh holder.sh

#actually download the files - turn the igets into a bash script
printf '#!/bin/bash\nset -e\n\n' > imeta.sh
cat holder.sh >> imeta.sh
rm holder.sh
bash imeta.sh

#do the md5sum check
for IRPATH in `grep "iget" imeta.sh | sed "s/iget //"`
do
	FID=`basename $IRPATH`
	MD5LOCAL=`md5sum $FID| cut -f -1 -d " "`
	MD5IRODS=`imeta ls -d $IRPATH md5 | grep 'value:' | sed 's/value: //'`
	if [ $MD5LOCAL != $MD5IRODS ]
	then
		echo "md5sum conflict encountered for $FID" 1>&2
		exit 1
	fi
done
rm imeta.sh

#CRAM to FASTQ conversion script is in the same folder as this one
DIR=`dirname "$0"`
#convert to FASTQ
parallel bash $DIR/cramfastq.sh ::: *.cram

#rename the resulting FASTQ files to be cellranger input friendly
#we're starting off with a file named like this: 22288_1#1.cram_I1_001.fastq.gz
#and we want to turn it into SAMPLE_S##_L00#_I1_001.fastq.gz
#we want the L to increment for each separate lot of CRAM files being processed
#and the S to increment for each separate file within that given lane
#(so from the example above, 22288_1 is the "L" and #1 is the "S")

#one global lane counter
lcount=1
for LANE in `ls *.cram | cut -f 1 -d "#" | sort | uniq`
do
	#reset sample counter per lane
	scount=1
	for FID in `ls $LANE*.cram`
	do
		#at this point, we've got FID to feature the full name of a CRAM
		#its corresponding reads will be FID_[I1/R1/R2]_001.fastq.gz
		#rename using our stored lcount and scount
		for READ in I1 R1 R2
		do
			mv $FID\_$READ\_001.fastq.gz $SAMPLE\_S$scount\_L00$lcount\_$READ\_001.fastq.gz
		done
		#bump sample counter
		(( scount++ ))
	done
	#bump lane counter
	(( lcount++ ))
done

#now that we have neatly named FASTQs, our job is done. remove CRAMs and call it a day
rm *.cram