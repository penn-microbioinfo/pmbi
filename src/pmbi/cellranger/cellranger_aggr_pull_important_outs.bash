#!/bin/bash

for i in $(ls *.tar); do 
    sample=$(echo $i | grep -Eo "[0-9]+[_][^.]+")
    rootdir=$(tar -tf $i | cut -d '/' -f1 | sort | uniq)
    if [[ "$rootdir" == "outs" ]]; then
        pathto="outs"
    elif [[ "$rootdir" == "$sample" ]]; then
        pathto="${sample}/outs"
    else
        echo "bad"
        exit 1
    fi
    
    IMPORTANT_OUTS=(\
       $pathto/per_sample_outs/${sample}/web_summary.html \
        $pathto/per_sample_outs/${sample}/count/sample_filtered_feature_bc_matrix.h5 \
        $pathto/multi/count/raw_feature_bc_matrix.h5 \
        $pathto/per_sample_outs/${sample}/vdj_t\
    )
    
    mkdir important_outs

    for imp_out in ${IMPORTANT_OUTS[@]}; do
        tar xvf $i $imp_out
        cp -r $imp_out important_outs/.
    done
    
    tar cvfz ${sample}.multi.importantOuts.tar.gz important_outs
    rm -r important_outs
    rm -r $rootdir
done
