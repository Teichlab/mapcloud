#!/bin/bash
set -e

#a script for downloading the CRAMs of a given sample, converting them to FASTQ, 
#and making the nomenclature be compatible with cellranger's expectations

#run with one positional argument - the sample ID
SAMPLE=$1

#acquire FASTQs via standard 10x proceedings
#start by querying iRODS for the sample's CRAMs
printf '#!/bin/bash\nset -e\n\n----\n' > imeta.sh
imeta qu -z seq -d sample = $SAMPLE and type = cram and target = 1 >> imeta.sh
#this used to be -iget -K, but there have been some bizarre ghost errors early Q2 2021
#yielding irreproducible borked reads at a rare rate without scripts breaking
#so as a precaution, do a manual md5sum check instead, just in case it's a ghost irods error
#(can skip the -K for run time as will check md5sums later)
sed ':a;N;$!ba;s/----\ncollection:/iget/g' -i imeta.sh
sed ':a;N;$!ba;s/\ndataObj: /\//g' -i imeta.sh
#clean up the imeta by kicking out *#888.crams and yhumans
grep -v "#888.cram\|yhuman" imeta.sh > temp.sh && mv temp.sh imeta.sh

#actually download the files
bash imeta.sh

#do the md5sum check
for IRPATH in `grep "iget" imeta.sh | sed "s/iget //"`
do
	FID=`basename $IRPATH`
	MD5LOCAL=`md5sum $FID| cut -f -1 -d " "`
	MD5IRODS=`imeta ls -d $IRPATH md5 | grep 'value:' | sed 's/value: //'`
	if [ $MD5LOCAL != $MD5IRODS ]
	then
		echo "md5sum conflict encountered for $FID"
		exit 1
	fi
done
rm imeta.sh

#convert to FASTQ
parallel bash /mnt/mapcloud/scripts/10x/utils/cramfastq.sh ::: *.cram

#rename the resulting FASTQ files to be cellranger input friendly
#we're starting off with a file named like this: 22288_1#1.cram_I1_001.fastq.gz
#and we want to turn it into SAMPLE_S##_L00#_I1_001.fastq.gz
#we want the L to increment for each separate lot of CRAM files being processed
#and the S to increment for each separate file within that given lane

#one global lane counter
lcount=1
for LANE in `ls *.cram | cut -f 1 -d "#"`
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