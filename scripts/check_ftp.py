import json, urllib.request, re

# PhenoAge 三个版本样本量
for a in ["GCST90014292","GCST90014298","GCST90014304"]:
    with urllib.request.urlopen(f"https://www.ebi.ac.uk/gwas/rest/api/studies/{a}", timeout=30) as r:
        d = json.load(r)
    print(a, "|", str(d.get("initialSampleSize","?"))[:60], "| cohorts:", str(d.get("cohort","?"))[:40])

def ftp_files(rng, acc):
    url = f"http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/{rng}/{acc}/"
    try:
        with urllib.request.urlopen(url, timeout=60) as r:
            html = r.read().decode()
        files = re.findall(r'href="([^"?/][^"]*)"', html)
        files = [f for f in files if not f.startswith("..")]
        print(acc, "->", files)
    except Exception as e:
        print(acc, "FTP ERROR", e)

print()
ftp_files("GCST90014001-GCST90015000", "GCST90014288")
ftp_files("GCST90014001-GCST90015000", "GCST90014289")
ftp_files("GCST90014001-GCST90015000", "GCST90014290")
ftp_files("GCST90014001-GCST90015000", "GCST90014292")
ftp_files("GCST90020001-GCST90021000", "GCST90020053")
