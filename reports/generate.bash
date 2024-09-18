#!/bin/bash
quarto render reports/scvdj.qmd --to html && cp /home/ubuntu/dev/pmbi/reports/scvdj.html /srv/http/betts/coculture/reports/.
