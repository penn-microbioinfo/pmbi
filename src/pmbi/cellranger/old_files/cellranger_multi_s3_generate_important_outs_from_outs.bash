#!/bin/bash


if [ -z $1 ]; then
	echo "ERROR: Must supply S3 URI as first positional argument."
	exit 1
fi

uri=$1
for i in $(aws s3 ls $uri | grep -oE "[0-9]+[_][0-9]+[a-z]*[.]multi[.]outs[.]tar[.]gz$"); do
	if [ ! -f ${i/.outs./.importantOuts.} ]; then
		aws s3 cp $(echo $uri | sed "s#/\$##")/${i} .
		bash /shared-ebs/microbioinfo-aws/scripts/cellranger_multi_pull_important_outs.bash
		rm $i
	fi
done
