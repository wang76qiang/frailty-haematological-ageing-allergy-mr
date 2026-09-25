# Canonical LDSC re-computation (bulik/ldsc v1.0.1, Python-3 compatibility patches only)
# Working directory: D:\衰老研究\v3_pipeline\data\r6\canonical_ldsc
import os, subprocess, sys

BASE = r"D:\衰老研究\v3_pipeline\data\r6\canonical_ldsc"
RAW = os.path.join(BASE, "raw")
MUNG = os.path.join(BASE, "munged")
LD = os.path.join(BASE, "ldscores")
LOGS = os.path.join(BASE, "logs")
HM3 = r"D:\衰老研究\v3_pipeline\data\r6\ldsc_ref\w_hm3.snplist.gz"
PY = r"E:\ldsc_venv\Scripts\python.exe"
LDSC = r"E:\ldsc_work\ldsc"

env = dict(os.environ)
env["PYTHONPATH"] = LDSC + ";" + os.path.join(LDSC, "ldscore")

TRAITS = {
    # name: (rawfile, extra munge args)
    "GrimAge": ("GCST90014288_buildGRCh37.tsv.gz",
                ["--snp", "variant_id", "--a1", "A1", "--a2", "A2", "--p", "p_value",
                 "--signed-sumstats", "Effect,0", "--frq", "Freq1"]),
    "Hannum": ("GCST90014289_buildGRCh37.tsv.gz",
               ["--snp", "variant_id", "--a1", "A1", "--a2", "A2", "--p", "p_value",
                "--signed-sumstats", "Effect,0", "--frq", "Freq1"]),
    "IEAA": ("GCST90014290_buildGRCh37.tsv.gz",
             ["--snp", "variant_id", "--a1", "A1", "--a2", "A2", "--p", "p_value",
              "--signed-sumstats", "Effect,0", "--frq", "Freq1"]),
    "PhenoAge": ("GCST90014292_buildGRCh37.tsv.gz",
                 ["--snp", "variant_id", "--a1", "A1", "--a2", "A2", "--p", "p_value",
                  "--signed-sumstats", "Effect,0", "--frq", "Freq1"]),
    "Frailty": ("GCST90020053_buildGRCh37.tsv",
                ["--snp", "variant_id", "--a1", "effect_allele", "--a2", "other_allele",
                 "--p", "p_value", "--signed-sumstats", "beta,0",
                 "--frq", "effect_allele_frequency", "--N", "175226"]),
    "ALLERG_ASTHMA": (r"D:\衰老研究\v3_pipeline\data\real\finngen_full\finngen_R12_ALLERG_ASTHMA.gz",
                      ["--snp", "rsids", "--a1", "alt", "--a2", "ref", "--p", "pval",
                       "--signed-sumstats", "beta,0", "--frq", "af_alt", "--N", "283740"]),
    "ALLERG_RHINITIS": (r"D:\衰老研究\v3_pipeline\data\real\finngen_full\finngen_R12_ALLERG_RHINITIS.gz",
                        ["--snp", "rsids", "--a1", "alt", "--a2", "ref", "--p", "pval",
                         "--signed-sumstats", "beta,0", "--frq", "af_alt", "--N", "490219"]),
    "L12_ATOPIC": (r"D:\衰老研究\v3_pipeline\data\real\finngen_full\finngen_R12_L12_ATOPIC.gz",
                   ["--snp", "rsids", "--a1", "alt", "--a2", "ref", "--p", "pval",
                    "--signed-sumstats", "beta,0", "--frq", "af_alt", "--N", "464119"]),
}

def run(cmd, logfile):
    with open(logfile, "w") as lf:
        p = subprocess.run(cmd, stdout=lf, stderr=subprocess.STDOUT, env=env)
    return p.returncode

if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "munge"
    if which == "munge":
        for name, (raw, args) in TRAITS.items():
            rawpath = raw if os.path.isabs(raw) else os.path.join(RAW, raw)
            out = os.path.join(MUNG, name)
            if os.path.exists(out + ".sumstats.gz"):
                print("skip", name); continue
            cmd = [PY, "-W", "ignore", os.path.join(LDSC, "munge_sumstats.py"),
                   "--sumstats", rawpath, "--merge-alleles", HM3,
                   "--out", out] + args
            log = os.path.join(LOGS, "munge_%s.log" % name)
            rc = run(cmd, log)
            print(name, "rc=", rc)
    elif which == "rg":
        w = os.path.join(LD, "1000G_Phase3_weights_hm3_no_MHC", "weights.hm3_noMHC.")
        traits = ",".join(os.path.join(MUNG, t + ".sumstats.gz") for t in TRAITS)
        cmd = [PY, "-W", "ignore", os.path.join(LDSC, "ldsc.py"),
               "--rg", traits,
               "--ref-ld-chr", w, "--w-ld-chr", w,
               "--not-M-5-50",
               "--out", os.path.join(BASE, "out", "canonical_rg")]
        log = os.path.join(LOGS, "ldsc_rg.log")
        rc = run(cmd, log)
        print("rg rc=", rc)
    elif which == "rg_all":
        w = os.path.join(LD, "1000G_Phase3_weights_hm3_no_MHC", "weights.hm3_noMHC.")
        outcomes = ["ALLERG_ASTHMA", "ALLERG_RHINITIS", "L12_ATOPIC"]
        for exp in ["GrimAge", "Hannum", "IEAA", "PhenoAge", "Frailty"]:
            traits = ",".join([os.path.join(MUNG, exp + ".sumstats.gz")] +
                              [os.path.join(MUNG, o + ".sumstats.gz") for o in outcomes])
            cmd = [PY, "-W", "ignore", os.path.join(LDSC, "ldsc.py"),
                   "--rg", traits,
                   "--ref-ld-chr", w, "--w-ld-chr", w,
                   "--not-M-5-50",
                   "--out", os.path.join(BASE, "out", "canonical_rg_%s" % exp)]
            log = os.path.join(LOGS, "ldsc_rg_%s.log" % exp)
            rc = run(cmd, log)
            print(exp, "rc=", rc)
    elif which == "rg_mhc_all":
        # sensitivity: LD scores that include the MHC region (LDscore.1kgPhase3.hm3)
        w = os.path.join(LD, "LDscore", "LDscore.")
        outcomes = ["ALLERG_ASTHMA", "ALLERG_RHINITIS", "L12_ATOPIC"]
        for exp in ["Frailty"]:
            traits = ",".join([os.path.join(MUNG, exp + ".sumstats.gz")] +
                              [os.path.join(MUNG, o + ".sumstats.gz") for o in outcomes])
            cmd = [PY, "-W", "ignore", os.path.join(LDSC, "ldsc.py"),
                   "--rg", traits,
                   "--ref-ld-chr", w, "--w-ld-chr", w,
                   "--out", os.path.join(BASE, "out", "canonical_rg_mhc_%s" % exp)]
            log = os.path.join(LOGS, "ldsc_rg_mhc_%s.log" % exp)
            rc = run(cmd, log)
            print(exp, "rc=", rc)
    elif which == "rg_mhc":
        # sensitivity: LD scores that include the MHC region (LDscore.1kgPhase3.hm3)
        w = os.path.join(LD, "LDscore", "LDscore.")
        traits = ",".join(os.path.join(MUNG, t + ".sumstats.gz") for t in TRAITS)
        cmd = [PY, "-W", "ignore", os.path.join(LDSC, "ldsc.py"),
               "--rg", traits,
               "--ref-ld-chr", w, "--w-ld-chr", w,
               "--out", os.path.join(BASE, "out", "canonical_rg_mhc")]
        log = os.path.join(LOGS, "ldsc_rg_mhc.log")
        rc = run(cmd, log)
        print("rg_mhc rc=", rc)
