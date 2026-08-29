#!/usr/bin/env python3
"""Stdlib-only Pan-UKBB streaming scanner (child process entry point).

Usage: _panukbb_scanner.py <component> <fname> <keys_file> <out_jsonl>
Scans data/real/panukbb/<fname> for chr:pos:ref:alt keys (EUR effect columns,
fallback meta_hq) and writes matched rows as JSON lines.  Intentionally free of
pandas/numpy so it stays tiny under concurrent multi-agent memory pressure.
"""

import gzip
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
PANU_DIR = os.path.join(ROOT, "data", "real", "panukbb")


def main():
    comp, fname, keys_file, out_jsonl = sys.argv[1:5]
    with open(keys_file) as fh:
        want = set(json.load(fh))
    full = os.path.join(PANU_DIR, fname)
    n_found = 0
    with gzip.open(full, "rt") as fh, open(out_jsonl, "w") as out:
        header = fh.readline().strip().split("\t")
        if "beta_EUR" in header:
            bcol, scol, pcol, acol = "beta_EUR", "se_EUR", "neglog10_pval_EUR", "af_EUR"
        else:
            bcol, scol, pcol, acol = ("beta_meta_hq", "se_meta_hq",
                                      "neglog10_pval_meta_hq", "af_meta_hq")
        idx = {"chrom": header.index("chr"), "pos": header.index("pos"),
               "ref": header.index("ref"), "alt": header.index("alt"),
               "beta": header.index(bcol), "se": header.index(scol),
               "p": header.index(pcol), "af": header.index(acol)}
        ncol = max(idx.values())
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= ncol:
                continue
            chrom = parts[idx["chrom"]].replace("chr", "")
            key = f"{chrom}:{parts[idx['pos']]}:{parts[idx['ref']]}:{parts[idx['alt']]}"
            if key not in want:
                continue
            try:
                beta = float(parts[idx["beta"]])
                se = float(parts[idx["se"]])
                logp = float(parts[idx["p"]])
                af = float(parts[idx["af"]])
            except ValueError:
                continue
            if not (se > 0):
                continue
            rec = {"snp": key, "chrom": chrom, "pos": int(parts[idx["pos"]]),
                   "ref": parts[idx["ref"]], "alt": parts[idx["alt"]],
                   "beta": beta, "se": se, "p": 10.0 ** (-logp),
                   "p_raw": logp, "af_alt": af, "component": comp}
            out.write(json.dumps(rec) + "\n")
            n_found += 1
    print(f"{comp}: {n_found} matches -> {out_jsonl}", flush=True)


if __name__ == "__main__":
    main()
