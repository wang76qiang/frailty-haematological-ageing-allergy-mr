import json, urllib.request

accs = ["GCST90014284","GCST90014285","GCST90014286","GCST90014287",
        "GCST90014288","GCST90014289","GCST90014290",
        "GCST90020053","GCST009856","GCST90278198"]
for a in accs:
    url = f"https://www.ebi.ac.uk/gwas/rest/api/studies/{a}"
    try:
        with urllib.request.urlopen(url, timeout=30) as r:
            d = json.load(r)
        trait = d.get("diseaseTrait", {}).get("trait", "?")
        n = str(d.get("initialSampleSize", "?"))[:45]
        rep = str(d.get("replicateSampleSize", "?"))[:30]
        print(f"{a} | {trait[:50]:50s} | {n:45s} | rep:{rep}")
    except Exception as e:
        print(f"{a} | ERROR {e}")
