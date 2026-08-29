import json
d = json.load(open(r'D:\衰老研究\v3_pipeline\data\r5\telomere_studies.json'))
ss = d.get('_embedded', {}).get('studies', [])
print('n_studies', len(ss))
for s in ss:
    pub = s.get('publicationInfo', {}) or {}
    print(s['accessionId'], '|',
          str(s.get('initialSampleSize', '?'))[:55], '|',
          'snp:', s.get('snpCount'), '|',
          'fullPvalSet:', s.get('fullPvalueSet'), '|',
          str(pub.get('title', '?'))[:50], '|',
          pub.get('publicationDate'))
