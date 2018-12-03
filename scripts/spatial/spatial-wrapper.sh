#!/bin/bash
set -e

#accepts four positional arguments:
#	* $1 - barcodes file from ~/st_pipeline/ids
#	* $2 - reference: GRCh38 or mm10
#	* $3 - sample or sample_id (or id_run or whatever), which field are we checking
#	* $4 - the value that we want to fish out

#set up STAR 2.5.4b to be in the path to make st_pipeline happy
PATHHOLD=$PATH
PATH=$PATH:/home/ubuntu/STAR-2.5.4b/bin/Linux_x86_64

#ok, time to get the imeta ball rolling
#prepare the foundation for the imeta dump becoming a shell script with igets
#and then obtain the imeta using the value pair(s) provided on input
printf '#!/bin/bash\nset -e\n\n----\n' > imeta.sh
imeta qu -z seq -d $3 = $4 and type = cram and target = 1 >> imeta.sh

#a couple quick substitutions to actually turn the imeta dump into a bunch of igets
#turn all '----\ncollection:' into 'iget'
#turn all '\ndataObj: ' into '/'
sed ':a;N;$!ba;s/----\ncollection:/iget -K/g' -i imeta.sh
sed ':a;N;$!ba;s/\ndataObj: /\//g' -i imeta.sh

#time to get downloading
bash imeta.sh

#and now that it's downloaded, we can dispose of it
rm imeta.sh

#the CRAMs are like SS2 CRAMs
mkdir fastq
parallel bash /mnt/mapcloud/scripts/ss2/cramfastq.sh ::: *.cram
rm *.cram

#dump all the r1's into R1.fastq.gz and all the r2's into R2.fastq.gz
#gzip is cool so we can just cat them together
cat fastq/*-1.fastq.gz > fastq/R1.fastq.gz
cat fastq/*-2.fastq.gz > fastq/R2.fastq.gz
rm fastq/*-1.fastq.gz && rm fastq/*-2.fastq.gz

#and now the pipeline proper can commence!
mkdir $4
st_pipeline_run.py --expName $4 --ids $1 --ref-map ~/cellranger/$2/star --log-file log_file.txt --output-folder $4 --ref-annotation ~/cellranger/$2/genes/genes.gtf fastq/R1.fastq.gz fastq/R2.fastq.gz
mv log_file.txt $4 && cd $4
convertEnsemblToNames.py --annotation ~/cellranger/$2/genes/genes.gtf --output $4\_stdata_genenames.tsv $4\_stdata.tsv
mkdir qc && cd qc
st_qa.py --input-data ../$4\_stdata_genenames.tsv

#clean up, including reverting the path
cd ../.. && rm -r fastq
PATH=$PATHHOLD
