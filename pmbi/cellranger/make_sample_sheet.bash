#!/bin/bash

sheet=$1

if [ -z $sheet ]; then
	echo "Must pass sheet as 1st positional argument."
	exit 1
fi
IFS=$'\n'; read -d '' -a sample_names < <(cat $sheet | cut -d , -f1,7,8 --output-delimiter='_')
IFS=$'\n'; read -d '' -a sample_indices < <(cat $sheet | cut -d , -f9,10 --output-delimiter "-" | sed 's/^/SI-/')
echo "Lane,Sample,Index"
for i in $(seq 1 $((${#sample_names[@]}-1))); do
	echo "*,${sample_names[i]},${sample_indices[i]}"
done
	
