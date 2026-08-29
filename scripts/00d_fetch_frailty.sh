#!/bin/bash
cd "/d/衰老研究/v3_pipeline"
U="http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90020001-GCST90021000/GCST90020053/harmonised/34431594-GCST90020053-EFO_0009885.h.tsv.gz"
echo "[$(date +%H:%M:%S)] download harmonised frailty start"
curl -s --retry 3 -C - -o data/r5/raw/frailty_h.tsv.gz "$U"
echo "[$(date +%H:%M:%S)] downloaded $(ls -la data/r5/raw/frailty_h.tsv.gz | awk '{print $5}') bytes; filtering"
zcat data/r5/raw/frailty_h.tsv.gz | awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)h[$i]=i; print "variant_id\tchromosome\tbase_pair_location\teffect_allele\tother_allele\teffect_allele_frequency\tbeta\tstandard_error\tp_value"; next} ($h["p_value"]+0) < 1e-5 {print $h["variant_id"]"\t"$h["chromosome"]"\t"$h["base_pair_location"]"\t"$h["effect_allele"]"\t"$h["other_allele"]"\t"$h["effect_allele_frequency"]"\t"$h["beta"]"\t"$h["standard_error"]"\t"$h["p_value"]}' > data/r5/Frailty_hits_p1e5.tsv
rm -f data/r5/raw/frailty_h.tsv.gz
echo "[$(date +%H:%M:%S)] frailty done -> $(wc -l < data/r5/Frailty_hits_p1e5.tsv) lines"
