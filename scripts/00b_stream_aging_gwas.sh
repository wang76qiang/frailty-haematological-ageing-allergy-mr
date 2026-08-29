#!/bin/bash
# R5-01: 流式读取 5 个正交衰老 GWAS（文件按 p 升序排列），p<1e-5 落表，连续 5000 行越阈即停
set -u
cd "/d/衰老研究/v3_pipeline"
mkdir -p data/r5/logs
BASE="http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics"

stream_subset() {
  local name="$1" path="$2" kind="$3"
  local out="data/r5/${name}_hits_p1e5.tsv"
  if [ -s "$out" ]; then echo "[$(date +%H:%M:%S)] $name exists ($(wc -l <"$out") lines), skip"; return; fi
  echo "[$(date +%H:%M:%S)] streaming $name ..."
  if [ "$kind" = "gz" ]; then
    curl -s --retry 3 "$BASE/$path" | zcat 2>/dev/null | awk -F'\t' '
      NR==1 {print; next}
      ($9+0) < 1e-5 {print; consec=0; next}
      {consec++; if (consec>5000) exit}' > "$out"
  else
    curl -s --retry 3 "$BASE/$path" | awk -F'\t' '
      NR==1 {print; next}
      ($9+0) < 1e-5 {print; consec=0; next}
      {consec++; if (consec>5000) exit}' > "$out"
  fi
  echo "[$(date +%H:%M:%S)] $name -> $out ($(wc -l <"$out") lines)"
}

stream_subset GrimAge  "GCST90014001-GCST90015000/GCST90014288/GCST90014288_buildGRCh37.tsv.gz" gz
stream_subset Hannum   "GCST90014001-GCST90015000/GCST90014289/GCST90014289_buildGRCh37.tsv.gz" gz
stream_subset IEAA     "GCST90014001-GCST90015000/GCST90014290/GCST90014290_buildGRCh37.tsv.gz" gz
stream_subset PhenoAge "GCST90014001-GCST90015000/GCST90014292/GCST90014292_buildGRCh37.tsv.gz" gz
stream_subset Frailty  "GCST90020001-GCST90021000/GCST90020053/GCST90020053_buildGRCh37.tsv" plain
echo "[$(date +%H:%M:%S)] ALL DONE"
