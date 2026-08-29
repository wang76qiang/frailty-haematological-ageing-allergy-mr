#!/usr/bin/env python3
"""
R5-13b: Three-source pharmacovigilance for type-2 allergy biologics.

Question: do drugs that block type-2 / alarmin axes (dupilumab, mepolizumab,
benralizumab, reslizumab, omalizumab, tezepelumab, lebrikizumab, tralokinumab)
show AGEING-related adverse-event signals (ageing / immunosenescence /
infection / malignancy) in spontaneous-report databases?  This is the safety
counterpart of R5-09 cis-pQTL drug-target MR: R5-09 asks whether lifelong
higher IL4/IL5/IL13/TSLP protein levels causally raise allergy risk (positive
controls); here we ask whether pharmacological blockade of the same axes is
accompanied by ageing-relevant AE reporting.

Sources
-------
(a) FAERS 2025Q4 (local ASCII): DRUG25Q4.txt / REAC25Q4.txt
(b) Canada Vigilance CVP extract (streamed from zip, one file at a time)
(c) ClinicalTrials.gov v2 JSON dumps (registration features; AE module if present)

Output
------
results/r5/tables/r5_pharmacovigilance_summary.csv   (long format, all sources)
results/r5/logs/r5_13b_summary.json
"""
import os
import re
import sys
import json
import glob
import zipfile

import numpy as np
import pandas as pd

BASE = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DATA_DIR = os.path.join(BASE, "data", "real")
OUT_TABLES = os.path.join(BASE, "results", "r5", "tables")
OUT_LOGS = os.path.join(BASE, "results", "r5", "logs")
CV_TMP = os.path.join(OUT_LOGS, "cv_tmp")
for d in (OUT_TABLES, OUT_LOGS, CV_TMP):
    os.makedirs(d, exist_ok=True)

FAERS_DRUG = os.path.join(DATA_DIR, "faers", "ASCII", "DRUG25Q4.txt")
FAERS_REAC = os.path.join(DATA_DIR, "faers", "ASCII", "REAC25Q4.txt")
CV_ZIP = os.path.join(DATA_DIR, "canada_vigilance", "extract_extrait.zip")
CV_PREFIX = "cvponline_extract_20241130/"
CT_DIR = os.path.join(DATA_DIR, "clinicaltrials")

CHUNK = 500_000
LOG = {"steps": []}

# drug -> (molecular target, R5-09 pQTL axis)
DRUG_TARGET = {
    "dupilumab":    ("IL4RA",  "IL4/IL13 axis (R5-09 positive control: IL4, IL13)"),
    "mepolizumab":  ("IL5",    "IL5 axis (R5-09 positive control: IL5)"),
    "reslizumab":   ("IL5",    "IL5 axis (R5-09 positive control: IL5)"),
    "benralizumab": ("IL5RA",  "IL5 axis (R5-09 positive control: IL5)"),
    "omalizumab":   ("IGHE",   "IgE (background; not in INTERVAL pQTL panel)"),
    "tezepelumab":  ("TSLP",   "TSLP axis (R5-09 positive control: TSLP)"),
    "lebrikizumab": ("IL13",   "IL13 axis (R5-09 positive control: IL13)"),
    "tralokinumab": ("IL13",   "IL13 axis (R5-09 positive control: IL13)"),
}
DRUGS = list(DRUG_TARGET)

# ageing-related AE keyword groups (MedDRA PT, lowercase word-boundary regex)
AE_GROUPS = {
    "ageing_general": [r"\bageing\b", r"\baging\b", r"\baged\b",
                       r"\bsenescen\w*\b", r"\belderly\b"],
    "immunosenescence": [r"immunosenescen\w*", r"immunodeficien\w*",
                         r"immunosuppress\w*", r"immune system disorder"],
    "infection": [r"\binfection\b", r"\binfections\b", r"\bsepsis\b",
                  r"\bseptic\w*\b", r"\bpneumonia\b", r"\btuberculosis\b",
                  r"\bherpes\w*\b", r"opportunistic infection"],
    "malignancy": [r"malignan\w*", r"\bcancer\b", r"carcinoma", r"lymphoma",
                   r"leukaemia", r"leukemia", r"melanoma", r"neoplasm",
                   r"sarcoma", r"myeloma", r"\btumou?r\b"],
}
AE_REGEX = {k: re.compile("|".join(v), re.I) for k, v in AE_GROUPS.items()}


def log(msg):
    print(msg, flush=True)
    LOG["steps"].append(msg)


def ror_2x2(a, b, c, d):
    """ROR with Haldane 0.5 correction and 95% CI."""
    ror = ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))
    se = np.sqrt(1 / (a + 0.5) + 1 / (b + 0.5) + 1 / (c + 0.5) + 1 / (d + 0.5))
    return ror, np.exp(np.log(ror) - 1.96 * se), np.exp(np.log(ror) + 1.96 * se)


def drug_report_sets_from_frame(df, name_col, id_col, drugs):
    sets = {d: set() for d in drugs}
    names = df[name_col].fillna("").str.lower()
    for d in drugs:
        m = names.str.contains(d, regex=False, na=False)
        if m.any():
            sets[d].update(df.loc[m, id_col].dropna().tolist())
    return sets


def ae_report_sets_from_frame(df, name_col, id_col):
    sets = {g: set() for g in AE_REGEX}
    names = df[name_col].fillna("").str.lower()
    for g, rgx in AE_REGEX.items():
        m = names.str.contains(rgx, na=False)
        if m.any():
            sets[g].update(df.loc[m, id_col].dropna().tolist())
    return sets


# --------------------------------------------------------------------------
# (a) FAERS
# --------------------------------------------------------------------------
def run_faers():
    log("FAERS: streaming DRUG25Q4/REAC25Q4 ...")
    drug_sets = {d: set() for d in DRUGS}
    n_ids = set()
    for chunk in pd.read_csv(FAERS_DRUG, sep="$", chunksize=CHUNK,
                             usecols=["primaryid", "drugname"], dtype=str,
                             low_memory=False, encoding="latin-1",
                             on_bad_lines="skip"):
        n_ids.update(chunk["primaryid"].dropna().tolist())
        part = drug_report_sets_from_frame(chunk, "drugname", "primaryid", DRUGS)
        for d in DRUGS:
            drug_sets[d] |= part[d]
    ae_sets = {g: set() for g in AE_REGEX}
    for chunk in pd.read_csv(FAERS_REAC, sep="$", chunksize=CHUNK,
                             usecols=["primaryid", "pt"], dtype=str,
                             low_memory=False, encoding="latin-1",
                             on_bad_lines="skip"):
        part = ae_report_sets_from_frame(chunk, "pt", "primaryid")
        for g in AE_REGEX:
            ae_sets[g] |= part[g]
    n_total = len(n_ids)
    log(f"FAERS N={n_total}; drug hits="
        f"{ {d: len(v) for d, v in drug_sets.items() if v} }")
    return build_ror_rows("FAERS", drug_sets, ae_sets, n_total), n_total


# --------------------------------------------------------------------------
# (b) Canada Vigilance
# --------------------------------------------------------------------------
def cv_extract(member):
    log(f"CV: extracting {member} ...")
    with zipfile.ZipFile(CV_ZIP) as z:
        z.extract(CV_PREFIX + member, path=CV_TMP)
    src = os.path.join(CV_TMP, CV_PREFIX + member)
    dst = os.path.join(CV_TMP, member)
    os.replace(src, dst)
    try:
        os.rmdir(os.path.join(CV_TMP, CV_PREFIX.rstrip("/")))
    except OSError:
        pass
    return dst


def cv_remove(path):
    try:
        os.remove(path)
        log(f"CV: deleted {path}")
    except OSError as e:
        log(f"WARN could not delete {path}: {e}")


def run_canada_vigilance():
    try:
        # N total from reports.txt (col0 = REPORT_ID)
        rp = os.path.join(CV_TMP, "reports.txt")
        if not os.path.exists(rp):
            rp = cv_extract("reports.txt")
        ids = set()
        for chunk in pd.read_csv(rp, sep="$", quotechar='"', header=None,
                                 usecols=[0], names=["rid"], chunksize=CHUNK,
                                 dtype=str, encoding="latin-1",
                                 on_bad_lines="skip", low_memory=False):
            ids.update(chunk["rid"].dropna().tolist())
        n_total = len(ids)
        cv_remove(rp)
        log(f"CV N={n_total}")

        # report_drug.txt: col1=REPORT_ID, col3=DRUGNAME
        rd = cv_extract("report_drug.txt")
        drug_sets = {d: set() for d in DRUGS}
        n = 0
        for chunk in pd.read_csv(rd, sep="$", quotechar='"', header=None,
                                 usecols=[1, 3], names=["rid", "drugname"],
                                 chunksize=CHUNK, dtype=str, encoding="latin-1",
                                 on_bad_lines="skip", low_memory=False):
            n += len(chunk)
            part = drug_report_sets_from_frame(chunk, "drugname", "rid", DRUGS)
            for d in DRUGS:
                drug_sets[d] |= part[d]
        cv_remove(rd)
        log(f"CV report_drug rows={n}; drug hits="
            f"{ {d: len(v) for d, v in drug_sets.items() if v} }")

        # reactions.txt: col1=REPORT_ID, col5=PT_NAME_ENG
        rx = cv_extract("reactions.txt")
        ae_sets = {g: set() for g in AE_REGEX}
        n = 0
        for chunk in pd.read_csv(rx, sep="$", quotechar='"', header=None,
                                 usecols=[1, 5], names=["rid", "pt"],
                                 chunksize=CHUNK, dtype=str, encoding="latin-1",
                                 on_bad_lines="skip", low_memory=False):
            n += len(chunk)
            part = ae_report_sets_from_frame(chunk, "pt", "rid")
            for g in AE_REGEX:
                ae_sets[g] |= part[g]
        cv_remove(rx)
        log(f"CV reactions rows={n}; AE hits="
            f"{ {g: len(v) for g, v in ae_sets.items() if v} }")
        return build_ror_rows("CanadaVigilance", drug_sets, ae_sets, n_total), n_total
    except Exception as e:
        import traceback
        log(f"FAILED CV: {e}")
        LOG["cv_error"] = traceback.format_exc()
        return [], 0


def build_ror_rows(source, drug_sets, ae_sets, n_total):
    rows = []
    for d, dset in drug_sets.items():
        target, axis = DRUG_TARGET[d]
        for g, gset in ae_sets.items():
            a = len(dset & gset)
            b = len(dset) - a
            c = len(gset) - a
            dd = max(n_total - a - b - c, 0)
            ror, lo, hi = ror_2x2(a, b, c, dd)
            rows.append({
                "source": source, "drug": d, "target": target,
                "r5_09_axis": axis, "ae_domain": g,
                "n_drug_reports": len(dset), "n_domain_reports": len(gset),
                "n_overlap_reports": a, "db_total_reports": n_total,
                "ror": round(ror, 4), "ror_lower": round(lo, 4),
                "ror_upper": round(hi, 4),
                # EMA-style rule: >=3 reports AND ROR>2 AND lower CI>1;
                # a=0 with a tiny domain inflates Haldane-corrected ROR spuriously
                "signal_ror_gt2": bool(a >= 3 and ror > 2 and lo > 1),
            })
    return rows


# --------------------------------------------------------------------------
# (c) ClinicalTrials.gov
# --------------------------------------------------------------------------
TRIAL_AE_RE = re.compile(
    r"hypersensitiv|anaphyla|urticaria|angioedema|pruritus|rash|eczema|"
    r"dermatitis|rhinitis|asthma|bronchospasm|wheez|allerg", re.I)


def run_trials():
    """Registration-feature description (+ AE module if ever present)."""
    rows = []
    for path in sorted(glob.glob(os.path.join(CT_DIR, "*_trials.json"))):
        drug = os.path.basename(path).replace("_trials.json", "")
        rec = {"source": "ClinicalTrials.gov", "drug": drug,
               "target": DRUG_TARGET.get(drug, ("", ""))[0],
               "r5_09_axis": DRUG_TARGET.get(drug, ("", ""))[1],
               "ae_domain": "n/a", "n_drug_reports": 0, "n_domain_reports": 0,
               "n_overlap_reports": 0, "db_total_reports": 0,
               "ror": np.nan, "ror_lower": np.nan, "ror_upper": np.nan,
               "signal_ror_gt2": False,
               "note": ""}
        try:
            with open(path) as f:
                data = json.load(f)
            studies = data.get("studies", [])
            rec["db_total_reports"] = len(studies)          # = n registered trials
            n_ae_mod, n_allergy_ev = 0, 0
            for s in studies:
                rs = s.get("resultsSection")
                if not rs:
                    continue
                ae = rs.get("adverseEventsModule")
                if not ae:
                    continue
                n_ae_mod += 1
                for bucket in ("seriousEvents", "otherEvents"):
                    for ev in ae.get(bucket, []) or []:
                        if TRIAL_AE_RE.search(str(ev.get("term", ""))):
                            n_allergy_ev += 1
            rec["n_drug_reports"] = n_ae_mod
            rec["n_overlap_reports"] = n_allergy_ev
            rec["note"] = ("no resultsSection/adverseEventsModule in dump; "
                           "registration counts only" if n_ae_mod == 0
                           else f"{n_ae_mod} trials with AE module")
        except Exception as e:
            rec["note"] = f"parse error: {e}"
        rows.append(rec)
    return rows


def main():
    try:
        rows = []
        faers_rows, faers_n = run_faers()
        rows += faers_rows
        cv_rows, cv_n = run_canada_vigilance()
        rows += cv_rows
        rows += run_trials()

        df = pd.DataFrame(rows)
        df.to_csv(os.path.join(OUT_TABLES, "r5_pharmacovigilance_summary.csv"),
                  index=False)

        # concise signals for the report
        sig = df[(df["source"].isin(["FAERS", "CanadaVigilance"])) &
                 (df["signal_ror_gt2"])]
        LOG["summary"] = {
            "faers_N": faers_n, "cv_N": cv_n,
            "n_rows": len(df),
            "signals": sig[["source", "drug", "ae_domain", "ror",
                            "ror_lower", "ror_upper"]].to_dict("records"),
        }
        with open(os.path.join(OUT_LOGS, "r5_13b_summary.json"), "w") as f:
            json.dump(LOG, f, indent=2, default=str)
        log(f"DONE R5-13b: {len(df)} rows, {len(sig)} ROR>2 signals")
    except Exception as e:
        import traceback
        log(f"FAILED R5-13b: {e}")
        with open(os.path.join(OUT_LOGS, "r5_13b_summary.json"), "w") as f:
            json.dump({"failed": str(e), "trace": traceback.format_exc(),
                       "steps": LOG["steps"]}, f, indent=2)
        raise


if __name__ == "__main__":
    main()
