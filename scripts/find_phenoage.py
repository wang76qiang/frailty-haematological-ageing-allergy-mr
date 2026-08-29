import json, urllib.request

# 1) 取 GCST90014288 的 pubmedId
with urllib.request.urlopen("https://www.ebi.ac.uk/gwas/rest/api/studies/GCST90014288", timeout=30) as r:
    d = json.load(r)
pmid = d.get("publicationInfo", {}).get("pubmedId")
print("pubmedId:", pmid, "| title:", d.get("publicationInfo", {}).get("title", "")[:80])

# 2) 按 pubmedId 列出该论文全部 study
url = f"https://www.ebi.ac.uk/gwas/rest/api/studies/search/findByPublicationIdPubmedId?pubmedId={pmid}&size=100"
with urllib.request.urlopen(url, timeout=60) as r:
    dd = json.load(r)
ss = dd.get("_embedded", {}).get("studies", [])
print("n_studies in paper:", len(ss))
for s in sorted(ss, key=lambda x: x["accessionId"]):
    print(s["accessionId"], "|", s.get("diseaseTrait", {}).get("trait", "?")[:60], "|", "fullPval:", s.get("fullPvalueSet"))
