#!/usr/bin/env python3
"""
R5-13: Three-source pharmacovigilance triangulation.

Sources
-------
(a) Canada Vigilance (CVP online extract) -> drug x PT ROR
(b) Existing FAERS results (results/tables/faers_allergy_signals.csv)
(c) ClinicalTrials.gov v2 JSON dumps (adverseEventsModule)

Outputs
-------
results/r5/tables/r5_canada_vigilance_ror.csv
results/r5/tables/r5_two_database_concordance.csv
results/r5/tables/r5_trial_ae_summary.csv
results/r5/tables/r5_target_vs_ae_direction.csv
results/r5/logs/r5_13_summary.json

Disk-safe: extracts one big CV file at a time into results/r5/logs/cv_tmp,
streams it, then deletes it before extracting the next.
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

CV_ZIP = os.path.join(DATA_DIR, "canada_vigilance", "extract_extrait.zip")
CV_PREFIX = "cvponline_extract_20241130/"
FAERS_CSV = os.path.join(BASE, "results", "tables", "faers_allergy_signals.csv")
RECON_CSV = os.path.join(BASE, "results", "r1", "tables",
                         "r1_drug_target_faers_reconciliation.csv")
CT_DIR = os.path.join(DATA_DIR, "clinicaltrials")
CONFIG_YAML = os.path.join(BASE, "config.yaml")

CHUNK = 500_000
LOG = {"steps": []}


def log(msg):
    print(msg, flush=True)
    LOG["steps"].append(msg)


def load_config():
    import yaml
    with open(CONFIG_YAML, encoding="utf-8") as f:
        cfg = yaml.safe_load(f)
    return cfg["faers"]["target_drugs"], cfg["faers"]["allergy_pts"]


# canonical PT groups: the 15 config PTs collapse into 11 MedDRA-level groups;
# regex synonyms extend matching (lowercase, word-boundary).
PT_PATTERNS = {
    "hypersensitivity": [r"\bhypersensitivity\b", r"\bdrug hypersensitivity\b",
                         r"\ballergic reaction\b", r"\bhypersensitivity reaction\b",
                         r"\ballergy\b"],
    "anaphylaxis": [r"\banaphylaxis\b", r"\banaphylactic reaction\b",
                    r"\banaphylactic shock\b", r"\banaphylactoid reaction\b",
                    r"\banaphylactic\b"],
    "urticaria": [r"\burticaria\b"],
    "angioedema": [r"\bangioedema\b", r"\bangioneurotic oedema\b"],
    "pruritus": [r"\bpruritus\b", r"\bitching\b"],
    "rash": [r"\brash\b", r"\bdrug eruption\b"],
    "eczema": [r"\beczema\b", r"\batopic dermatitis\b", r"\bdermatitis atopic\b"],
    "allergic_rhinitis": [r"\ballergic rhinitis\b", r"\brhinitis allergic\b",
                          r"\brhinitis\b", r"\bhay ?fever\b",
                          r"\bseasonal allergy\b"],
    "asthma": [r"\basthma\b", r"\basthmatic\b"],
    "bronchospasm": [r"\bbronchospasm\b"],
    "wheezing": [r"\bwheez\w*\b"],
}
PT_REGEX = {k: re.compile("|".join(v), re.I) for k, v in PT_PATTERNS.items()}


def extract_member(member):
    """Extract a single zip member to CV_TMP and return its path."""
    log(f"Extracting {member} ...")
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


def remove_file(path):
    try:
        os.remove(path)
        log(f"Deleted {path} to save disk")
    except OSError as e:
        log(f"WARN could not delete {path}: {e}")


def count_reports(path):
    """reports.txt col0 = REPORT_ID; count unique."""
    ids = set()
    for chunk in pd.read_csv(path, sep="$", quotechar='"', header=None,
                             usecols=[0], names=["rid"], chunksize=CHUNK,
                             encoding="latin-1", on_bad_lines="skip",
                             dtype=str, low_memory=False):
        ids.update(chunk["rid"].dropna().tolist())
    log(f"reports.txt unique REPORT_ID = {len(ids)}")
    return ids


def scan_report_drug(path, drugs):
    """report_drug.txt: col1=REPORT_ID, col3=DRUGNAME -> {drug: set(report_id)}."""
    drug_reports = {d: set() for d in drugs}
    pat = {d: re.compile(re.escape(d), re.I) for d in drugs}
    n = 0
    for chunk in pd.read_csv(path, sep="$", quotechar='"', header=None,
                             usecols=[1, 3], names=["rid", "drugname"],
                             chunksize=CHUNK, encoding="latin-1",
                             on_bad_lines="skip", dtype=str, low_memory=False):
        n += len(chunk)
        chunk["drugname"] = chunk["drugname"].fillna("").str.lower()
        for d in drugs:
            m = chunk["drugname"].str.contains(d, regex=False, na=False)
            if m.any():
                drug_reports[d].update(chunk.loc[m, "rid"].dropna().tolist())
        if n % 2_000_000 < CHUNK:
            hit = sum(len(v) for v in drug_reports.values())
            log(f"  report_drug {n/1e6:.1f}M rows, matched rows so far={hit}")
    total = {d: len(v) for d, v in drug_reports.items()}
    log(f"report_drug scan done ({n} rows). Reports per drug: "
        f"{ {d: total[d] for d in drugs if total[d]} }")
    return drug_reports


def scan_reactions(path):
    """reactions.txt: col1=REPORT_ID, col5=PT_NAME_ENG -> {pt: set(report_id)}."""
    pt_reports = {pt: set() for pt in PT_REGEX}
    n = 0
    for chunk in pd.read_csv(path, sep="$", quotechar='"', header=None,
                             usecols=[1, 5], names=["rid", "pt"],
                             chunksize=CHUNK, encoding="latin-1",
                             on_bad_lines="skip", dtype=str, low_memory=False):
        n += len(chunk)
        names = chunk["pt"].fillna("").str.lower()
        for pt, rgx in PT_REGEX.items():
            m = names.str.contains(rgx, na=False)
            if m.any():
                pt_reports[pt].update(chunk.loc[m, "rid"].dropna().tolist())
        if n % 2_000_000 < CHUNK:
            log(f"  reactions {n/1e6:.1f}M rows scanned")
    log(f"reactions scan done ({n} rows). PT hit reports: "
        f"{ {pt: len(v) for pt, v in pt_reports.items() if v} }")
    return pt_reports


def ror_table(drug_reports, pt_reports, n_total):
    """Drug x PT 2x2 tables -> ROR with Haldane 0.5 correction."""
    rows = []
    for d, drep in drug_reports.items():
        nd = len(drep)
        for pt, prep in pt_reports.items():
            a = len(drep & prep)
            b = nd - a
            c = len(prep) - a
            dd = n_total - a - b - c
            if dd < 0:
                dd = 0
            ror = ((a + 0.5) * (dd + 0.5)) / ((b + 0.5) * (c + 0.5))
            se = np.sqrt(1 / (a + 0.5) + 1 / (b + 0.5) + 1 / (c + 0.5) + 1 / (dd + 0.5))
            lo, hi = np.exp(np.log(ror) - 1.96 * se), np.exp(np.log(ror) + 1.96 * se)
            rows.append({"drug": d, "pt": pt, "a": a, "b": b, "c": c, "d": dd,
                         "n_drug_reports": nd, "n_pt_reports": len(prep),
                         "ror": ror, "ror_lower": lo, "ror_upper": hi,
                         "signal_ror_gt2": bool(ror > 2 and lo > 1)})
    return pd.DataFrame(rows)


def run_canada_vigilance(drugs):
    try:
        # 1) reports.txt -> N total (already extracted; else extract)
        reports_path = os.path.join(CV_TMP, "reports.txt")
        if not os.path.exists(reports_path):
            reports_path = extract_member("reports.txt")
        report_ids = count_reports(reports_path)
        n_total = len(report_ids)
        remove_file(reports_path)

        # 2) report_drug.txt -> drug report sets
        rd_path = extract_member("report_drug.txt")
        drug_reports = scan_report_drug(rd_path, drugs)
        remove_file(rd_path)

        # 3) reactions.txt -> PT report sets
        rx_path = extract_member("reactions.txt")
        pt_reports = scan_reactions(rx_path)
        remove_file(rx_path)

        ror = ror_table(drug_reports, pt_reports, n_total)
        ror.to_csv(os.path.join(OUT_TABLES, "r5_canada_vigilance_ror.csv"),
                   index=False)
        n_hit = int((ror["n_drug_reports"] > 0).groupby(ror["drug"]).max().sum())
        LOG["cv"] = {
            "n_total_reports": n_total,
            "drugs_with_reports": {d: len(v) for d, v in drug_reports.items() if v},
            "pt_hit_reports": {pt: len(v) for pt, v in pt_reports.items() if v},
        }
        log(f"CV ROR table written: {len(ror)} drug x PT rows; "
            f"{n_hit}/{len(drugs)} drugs with >=1 report")
        return ror
    except Exception as e:
        import traceback
        log(f"FAILED CV step: {e}")
        LOG["cv_error"] = traceback.format_exc()
        return pd.DataFrame()


def run_concordance(cv_ror):
    """Merge CV max-ROR per drug with FAERS signals; flag ROR>2 in both."""
    try:
        faers = pd.read_csv(FAERS_CSV)
        faers = faers.rename(columns={"ror": "faers_ror",
                                      "ror_lower": "faers_ror_lower",
                                      "ror_upper": "faers_ror_upper",
                                      "signal": "faers_signal",
                                      "n_reports": "faers_n_reports",
                                      "n_allergy_events": "faers_n_events"})
        if not cv_ror.empty:
            idx = cv_ror.groupby("drug")["ror"].idxmax()
            cv_best = cv_ror.loc[idx, ["drug", "pt", "ror", "ror_lower", "ror_upper", "a"]]
            cv_best = cv_best.rename(columns={"pt": "cv_top_pt", "ror": "cv_max_ror",
                                              "ror_lower": "cv_ror_lower",
                                              "ror_upper": "cv_ror_upper",
                                              "a": "cv_top_pt_n"})
        else:
            cv_best = pd.DataFrame(columns=["drug", "cv_top_pt", "cv_max_ror",
                                            "cv_ror_lower", "cv_ror_upper",
                                            "cv_top_pt_n"])
        m = faers.merge(cv_best, on="drug", how="left")
        # robust CV signal: ROR>2 AND CI lower>1 AND >=1 observed event
        # (avoids Haldane 0.5-correction artefacts on zero-count cells)
        m["cv_signal"] = ((m["cv_max_ror"] > 2) & (m["cv_ror_lower"] > 1)
                          & (m["cv_top_pt_n"].fillna(0) > 0))
        m["faers_ror_gt2"] = m["faers_ror"] > 2

        def verdict(r):
            if r["faers_ror_gt2"] and r["cv_signal"]:
                return "concordant_signal"
            if r["faers_ror_gt2"] and pd.isna(r["cv_max_ror"]):
                return "faers_only_cv_no_data"
            if r["faers_ror_gt2"]:
                return "faers_only"
            if r["cv_signal"]:
                return "cv_only"
            if pd.isna(r["cv_max_ror"]):
                return "neither_cv_no_data"
            return "neither"

        m["concordance"] = m.apply(verdict, axis=1)
        m["concordant"] = m["concordance"] == "concordant_signal"
        m = m.sort_values(["concordant", "faers_ror"], ascending=[False, False])
        m.to_csv(os.path.join(OUT_TABLES, "r5_two_database_concordance.csv"),
                 index=False)
        n_conc = int(m["concordant"].sum())
        LOG["concordance"] = {
            "n_concordant": n_conc,
            "concordant_drugs": m.loc[m["concordant"], "drug"].tolist(),
            "verdict_counts": m["concordance"].value_counts().to_dict(),
        }
        log(f"Concordance: {n_conc} drugs with ROR>2 in both databases")
        return m
    except Exception as e:
        import traceback
        log(f"FAILED concordance step: {e}")
        LOG["concordance_error"] = traceback.format_exc()
        return pd.DataFrame()


TRIAL_AE_RE = re.compile(
    r"hypersensitiv|anaphyla|urticaria|angioedema|pruritus|rash|eczema|"
    r"dermatitis|rhinitis|asthma|bronchospasm|wheez|allerg", re.I)


def run_trials():
    """Summarise allergy-related AEs in ClinicalTrials.gov JSON dumps."""
    rows = []
    for path in sorted(glob.glob(os.path.join(CT_DIR, "*_trials.json"))):
        drug = os.path.basename(path).replace("_trials.json", "")
        rec = {"drug": drug, "n_trials": 0, "n_with_results": 0,
               "n_with_ae_module": 0, "n_allergy_events": 0,
               "n_allergy_affected": 0, "allergy_terms": "", "note": ""}
        try:
            with open(path) as f:
                data = json.load(f)
            studies = data.get("studies", [])
            rec["n_trials"] = len(studies)
            terms = {}
            for s in studies:
                rs = s.get("resultsSection")
                if not rs:
                    continue
                rec["n_with_results"] += 1
                ae = rs.get("adverseEventsModule")
                if not ae:
                    continue
                rec["n_with_ae_module"] += 1
                for bucket in ("seriousEvents", "otherEvents"):
                    for ev in ae.get(bucket, []) or []:
                        term = str(ev.get("term", ""))
                        if not TRIAL_AE_RE.search(term):
                            continue
                        aff = sum(int(st_.get("numAffected", 0) or 0)
                                  for st_ in ev.get("stats", []) or [])
                        terms[term] = terms.get(term, 0) + max(aff, 1)
            if terms:
                rec["n_allergy_events"] = len(terms)
                rec["n_allergy_affected"] = int(sum(terms.values()))
                rec["allergy_terms"] = "; ".join(
                    f"{t}({n})" for t, n in
                    sorted(terms.items(), key=lambda x: -x[1])[:15])
            if rec["n_with_ae_module"] == 0:
                rec["note"] = "no resultsSection/adverseEventsModule in dump"
        except Exception as e:
            rec["note"] = f"parse error: {e}"
        rows.append(rec)
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(OUT_TABLES, "r5_trial_ae_summary.csv"), index=False)
    LOG["trials"] = {
        "n_drugs": len(df),
        "total_trials": int(df["n_trials"].sum()),
        "with_ae_module": int(df["n_with_ae_module"].sum()),
    }
    log(f"CT.gov: {len(df)} drugs, {df['n_trials'].sum()} trials, "
        f"{df['n_with_ae_module'].sum()} with AE module")
    return df


def run_target_vs_ae(concordance):
    """Text-level target-vs-AE direction table built on the R1 reconciliation."""
    try:
        recon = pd.read_csv(RECON_CSV)
    except Exception as e:
        log(f"WARN reconciliation table missing: {e}; using FAERS-only skeleton")
        recon = pd.DataFrame()
    if concordance is None or concordance.empty:
        concordance = pd.DataFrame(columns=["drug", "cv_max_ror", "concordance"])
    cv_map = concordance.set_index("drug")["cv_max_ror"].to_dict() \
        if "cv_max_ror" in concordance else {}
    conc_map = concordance.set_index("drug")["concordance"].to_dict() \
        if "concordance" in concordance else {}

    rows = []
    for _, r in recon.iterrows():
        drugs = [d.strip() for d in str(r.get("FAERS_drugs", "")).split(",")
                 if d.strip() and d.strip() != "nan"]
        cv_vals = [cv_map[d] for d in drugs if d in cv_map and pd.notna(cv_map[d])]
        cv_max = max(cv_vals) if cv_vals else np.nan
        conc_drugs = [d for d in drugs
                      if conc_map.get(d) == "concordant_signal"]
        mr_dir = r.get("MR_direction", "")
        faers_sig = bool(r.get("FAERS_signal", False))
        faers_ror = r.get("FAERS_max_ROR", np.nan)
        # direction verdict (text level)
        if mr_dir == "risk-enhancing" and r.get("MR_significant") is True:
            verdict = ("MR predicts target activation increases allergy risk; "
                       "expect inhibitor to LOWER allergy AE reporting")
        elif mr_dir == "protective" and r.get("MR_significant") is True:
            verdict = ("MR predicts target activation lowers allergy risk; "
                       "inhibitor may INCREASE allergy AE reporting")
        else:
            verdict = "MR null: no directional AE expectation"
        if faers_sig:
            verdict += (" | FAERS elevated (ROR>2) - likely indication/"
                        "channeling confounding")
        if conc_drugs:
            verdict += f" | CV concordant: {', '.join(conc_drugs)}"
        rows.append({
            "gene": r.get("gene"), "MR_OR": r.get("MR_OR"),
            "MR_p": r.get("MR_p"), "MR_significant": r.get("MR_significant"),
            "MR_direction": mr_dir,
            "FAERS_signal": faers_sig, "FAERS_max_ROR": faers_ror,
            "FAERS_drugs": r.get("FAERS_drugs"),
            "CV_max_ROR": cv_max,
            "CV_concordant_drugs": "; ".join(conc_drugs),
            "direction_verdict": verdict,
        })
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(OUT_TABLES, "r5_target_vs_ae_direction.csv"),
              index=False)
    log(f"Target-vs-AE direction table written: {len(df)} targets")
    return df


def main():
    drugs, cfg_pts = load_config()
    log(f"Loaded config: {len(drugs)} drugs, {len(cfg_pts)} config PTs "
        f"(-> {len(PT_PATTERNS)} canonical groups with synonyms)")
    LOG["config"] = {"drugs": drugs, "config_pts": cfg_pts,
                     "canonical_pt_groups": list(PT_PATTERNS)}

    cv_ror = run_canada_vigilance(drugs)
    conc = run_concordance(cv_ror)
    run_trials()
    run_target_vs_ae(conc)

    with open(os.path.join(OUT_LOGS, "r5_13_summary.json"), "w") as f:
        json.dump(LOG, f, indent=2, default=str)
    log("DONE R5-13")


if __name__ == "__main__":
    main()
