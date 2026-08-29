#!/bin/bash
cd "/d/衰老研究/v3_pipeline"
echo "[$(date +%H:%M:%S)] frailty full-file stream start"
curl -s --retry 3 "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90020001-GCST90021000/GCST90020053/GCST90020053_buildGRCh37.tsv" | awk -F'\t' 'NR==1 || ($9+0) < 1e-5' > data/r5/Frailty_hits_p1e5.tsv
echo "[$(date +%H:%M:%S)] frailty done -> $(wc -l < data/r5/Frailty_hits_p1e5.tsv) lines"
