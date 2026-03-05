#!/bin/bash

fragments=atac_fragments.tsv.gz
singlecell=atac_singlecell.csv
chromosomes=~/super1/ref/genomes/Mmul_10__SHIV.C.CH505.v2/Mmul_10_autosomes.txt
rm_bed=~/super1/ref/genomes/Mmul_10__SHIV.C.CH505.v2/GCF_003339765.1_Mmul_10_rm.bed
outdir=amulet
amulet_dir=~/dev/AMULET_wuv21
iscellidx=1

wd=/home/amsesk/super2/jayme_shiv/cr_arc_runs/
rundir_pat="^[0-9]$"

for rundir in $(readlink -f $(ls $wd | grep "$rundir_pat")); do
    outs=$rundir/outs
    if [[ ! -d $outs/$outdir ]]; then
        mkdir $outs/$outdir
    fi
    echo "Working on ${rundir}"
    $amulet_dir/AMULET.sh \
        $outs/$fragments \
        $outs/$singlecell \
        $chromosomes \
        $rm_bed \
        $outs/$outdir \
        $amulet_dir \
        --iscellidx $iscellidx
done
# amulet $fragments $singlecell $chromosomes $rm_bed $outdir $amulet_dir  --iscellidx $iscellidx
