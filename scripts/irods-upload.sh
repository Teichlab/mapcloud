#!/bin/bash
set -e

#run with 10X/10X-VDJ/SS2 as the first argument and the sample/runlane as the second argument

iput -Kr $2 /archive/HCA/$1
if [ $1 == 'SS2' ]
then
	imeta addw -d /archive/HCA/$1/$2/%/% target mapcloud
	run=`echo $2 | sed -e 's/_.*//'`
	lane=`echo $2 | sed -e 's/.*_//'`
	imeta addw -d /archive/HCA/$1/$2/%/% id_run $run
	imeta addw -d /archive/HCA/$1/$2/%/% lane $lane
	#bit of an ugly hack as regular imeta queries take forever for SS2 for some reason
	studyid=`imeta ls -d /seq/$run/$2#1.cram | grep -A1 'attribute: study_id' | tail -n 1 | sed 's/value: //'`
else
	#needs a subfolder specified, as trying to add metadata when a different tagged subfolder exists fails
	SUB=`ls $2 | sed "s|^|/|"`
	imeta addw -d /archive/HCA/$1/$2$SUB/%/% target mapcloud
	imeta addw -d /archive/HCA/$1/$2$SUB/%/% sample $2
	studyid=`imeta qu -z seq -d sample = $2 and type = cram and target = 1 | head -n 2 | sed 's/collection:/ls -d/' | sed ':a;N;$!ba;s/\ndataObj: /\//g' | xargs imeta | grep -A1 'attribute: study_id' | tail -n 1 | sed 's/value: //'`
fi
ichmod -r own kp9#Sanger1 /archive/HCA/$1/$2
ichmod -r read ss_$studyid#archive /archive/HCA/$1/$2
