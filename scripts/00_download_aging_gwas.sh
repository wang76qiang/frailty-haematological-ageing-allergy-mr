#!/bin/bash
# R5-01: 下载 5 个正交衰老 GWAS 汇总统计并流式过滤 p<1e-5 子集（原始大文件立即删除）
set -u
cd "/d/衰老研究/v3_pipeline"
mkdir -p data/r5/raw data/r5/logs
BASE="http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics"

process() {
  local name="$1" path="$2" kind="$3"
  local out="data/r5/${name}_hits_p1e5.tsv"
  if [ -s "$out" ]; then echo "[$(date +%H:%M:%S)] $name subset exists ($(wc -l < "$out") lines), skip"; return; fi
  local raw="data/r5/raw/${name}.tmp"
  echo "[$(date +%H:%M:%S)] downloading $name ..."
  curl -s --retry 3 --retry-delay 5 -C - -o "$raw" "$BASE/$path"
  echo "[$(date +%H:%M:%S)] filtering $name (p<1e-5) ..."
  if [ "$kind" = "gz" ]; then
    zcat "$raw" | awk 'NR==1 || ($9+0) < 1e-5' > "$out"
  else
    awk 'NR==1 || ($9+0) < 1e-5' "$raw" > "$out"
  fi
  rm -f "$raw"
  echo "[$(date +%H:%M:%S)] $name done -> $out ($(wc -l < "$out") lines)"
}

process GrimAge  "GCST90014001-GCST90015000/GCST90014288/GCST90014288_buildGRCh37.tsv.gz" gz
process Hannum   "GCST90014001-GCST90015000/GCST90014289/GCST90014289_buildGRCh37.tsv.gz" gz
process IEAA     "GCST90014001-GCST90015000/GCST90014290/GCST90014290_buildGRCh37.tsv.gz" gz
process PhenoAge "GCST90014001-GCST90015000/GCST90014292/GCST90014292_buildGRCh37.tsv.gz" gz
process Frailty  "GCST90020001-GCST90021000/GCST90020053/GCST90020053_buildGRCh37.tsv" plain
echo "[$(date +%H:%M:%S)] ALL DONE"
