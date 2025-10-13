#!/bin/bash

model_ls_pat=$1

for md in $(ls -d $model_ls_pat); do
    h5ad=$(echo $md | sed -r "s#[.]_n_latent_[0-9]+_model[/]##").h5ad
    margll_out=$(echo $md | sed "s/\\/$/_margll/")
    echo "$h5ad -- $margll_out"
    python ~/dev/pmbi/scripts/get_marginal_ll.py -a $h5ad -m $md -c 15000 > $margll_out
done
