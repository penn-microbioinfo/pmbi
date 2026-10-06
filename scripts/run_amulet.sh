#!/bin/bash

fragments=fragments.tsv.gz
singlecell=singlecell.csv
chromosomes=~/super1/ref/genomes/Mmul_10__SHIV.C.CH505.v2/Mmul_10_autosomes.txt
rm_bed=~/super1/ref/genomes/Mmul_10__SHIV.C.CH505.v2/GCF_003339765.1_Mmul_10_rm.bed
outdir=amulet
amulet_dir=~/dev/AMULET_wuv21

# For cellranger-arc
# iscellidx=1

# For cellranger-atac
iscellidx=9

wd=/home/amsesk/super2/jayme_shiv/cr_atac_runs_symlinks/
rundir_pat="^[0-9]+$"

for rundir in $(ls $wd | grep -E "$rundir_pat"); do
    rundir=${wd}/${rundir}
    outs=${rundir}/outs
    if [[ ! -d $outs/$outdir ]]; then
        mkdir $outs/$outdir
    fi
    echo "Working on ${rundir}"
    # echo " \
    $amulet_dir/AMULET.sh \
        $outs/$fragments \
        $outs/$singlecell \
        $chromosomes \
        $rm_bed \
        $outs/$outdir \
        $amulet_dir \
        --iscellidx $iscellidx
    # "
done
# amulet $fragments $singlecell $chromosomes $rm_bed $outdir $amulet_dir  --iscellidx $iscellidx
