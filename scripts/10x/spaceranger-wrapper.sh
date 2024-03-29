#!/bin/bash
set -eo pipefail

#four (or six) positional arguments on launch, 3 onwards for use on the imeta call
#	* $1 - fastq or cram, depending on what file to get from irods
#	* $2 - GRCh38 or mm10, which reference to use
#	* $3 - metadata file, formatted as detailed in the visium runscript
#	* $4 - sample or sample_id (or id_run or whatever), which field are we checking
#	* $5 - the value that we want to fish out
#	* optionally, $6 and $7 for another pair

#ok, time to get the imeta ball rolling
#prepare the foundation for the imeta dump becoming a shell script with igets
#and then obtain the imeta using the value pair(s) provided on input
printf '#!/bin/bash\nset -e\n\n----\n' > imeta.sh
if [[ $# -eq 5 ]]
then
	imeta qu -z seq -d $4 = $5 and type = $1 and target = 1 >> imeta.sh
fi
if [[ $# -eq 7 ]]
then
	imeta qu -z seq -d $4 = $5 and $6 = $7 and type = $1 and target = 1 >> imeta.sh
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

if [[ $1 == cram ]]
then
	#process the CRAMs into FASTQs
	parallel bash /mnt/mapcloud/scripts/10x/utils/cramfastq.sh ::: *.cram

	#rename the resulting FASTQ files to be cellranger input friendly
	#we're starting off with a file named like this: 22288_1#1.cram_I1_001.fastq.gz
	#and we want to turn it into SAMPLENAME_S##_L001_I1_001.fastq.gz
	#the ## in the S doesn't matter, we just want however many unique values for however many unique files
	#while staying consistent across the same block of I1/R1/R2

	#so let's just pick up the I1's as the block reps, and swap everything up to the .cram
	#to samplename_S*COUNT*_L001. that'll do. but we have to keep a count going manually

	count=1
	for file in *I1_001.fastq.gz
	do
		mv "${file}" "${file/*.cram/$5\_S$count\_L001}"
		file=$(sed 's/I1/R1/g' <<< $file)
		mv "${file}" "${file/*.cram/$5\_S$count\_L001}"
		file=$(sed 's/R1/R2/g' <<< $file)
		mv "${file}" "${file/*.cram/$5\_S$count\_L001}"
		(( count++ ))
	done
else
	#rename the fastqs to a format that cellranger will get
	rename s/[^_]*_[^_]*_[^_]*_/$5\_/ *.fastq.gz
fi

#move the files into a subdirectory and clean the CRAMs up
mkdir fastq
mv *.fastq.gz fastq
rm -f *.cram

#spaceranger proper!
#start by stripping out appropriate stuff from the metadata
IMAGE=`grep $5 $3 | cut -f 2 -d ,`
if [[ `head -n 1 $3 | tr -cd , | wc -c` == 3 ]]
then
	#standard metadata file
	SLIDE=`grep $5 $3 | cut -f 3 -d ,`
	AREA=`grep $5 $3 | cut -f 4 -d ,`
	~/cellranger/spaceranger-1.1.0/spaceranger count --id=$5 --fastqs=fastq --transcriptome=/home/ubuntu/cellranger/$2 --reorient-images --image=$IMAGE --slide=$SLIDE --area=$AREA
elif [[ `head -n 1 $3 | tr -cd , | wc -c` == 4 ]]
then
	#extra column for manual alignment
	SLIDE=`grep $5 $3 | cut -f 3 -d ,`
	AREA=`grep $5 $3 | cut -f 4 -d ,`
	ALIGNMENT=`grep $5 $3 | cut -f 5 -d ,`
	~/cellranger/spaceranger-1.1.0/spaceranger count --id=$5 --fastqs=fastq --transcriptome=/home/ubuntu/cellranger/$2 --image=$IMAGE --slide=$SLIDE --area=$AREA --loupe-alignment=$ALIGNMENT
else
	#sometimes people drop the ball and don't save the slide serial number
	~/cellranger/spaceranger-1.1.0/spaceranger count --id=$5 --fastqs=fastq --transcriptome=/home/ubuntu/cellranger/$2 --reorient-images --image=$IMAGE --unknown-slide
fi

#wipe out temporary files as those cause nothing but trouble
ls -d $5/* | grep -v 'outs' | xargs rm -r

#leave the fastq input and complete output available to be dealt with in whatever manner
