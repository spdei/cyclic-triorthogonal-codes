#!/usr/bin/env python3
"""Cross-check the binary classification against the independent Magma search (magma/*.csv).

For every odd length n <= 763 present in either, compares the sets of maximal admissible spectral
supports up to multiplier equivalence, their dimensions, triple-evenness, and exact distances
(Magma computes these for n <= 227)."""
import csv, glob, json, os, sys
from math import gcd
ROOT = os.path.join(os.path.dirname(__file__), '..')
def canon(K, n): return min(tuple(sorted((u*x) % n for x in K)) for u in range(1, n) if gcd(u, n) == 1)
def parseK(s): return [int(t) for t in s.strip('{} ').split(',') if t.strip()]
magma = {}
for f in sorted([os.path.join(ROOT, 'magma', 'binary_search.csv')]):
    for r in csv.DictReader(open(f)):
        n = int(r['n'])
        if n <= 763:
            key = canon(parseK(r['Roots_K']), n)
            old = magma.setdefault(n, {}).get(key)
            if old is None or int(r['Exact_dZ']) > int(old['Exact_dZ']):
                magma[n][key] = r
mine = {int(k): v for k, v in json.load(open(os.path.join(ROOT, 'results', 'binary_results.json'))).items()}
mine = {n: v for n, v in mine.items() if v['codes']}
bad = []
if set(magma) != set(mine): bad.append(('lengths differ', sorted(set(magma) ^ set(mine))))
agree_d = 0
for n in sorted(set(magma) & set(mine)):
    Y = {tuple(c['S']): c for c in mine[n]['codes']}
    if set(magma[n]) != set(Y): bad.append(('supports', n)); continue
    for key, r in magma[n].items():
        if int(r['Dimension']) != Y[key]['dim']: bad.append(('dim', n))
        if r['Is_Triply_Even'].strip() != 'true': bad.append(('not TE', n))
        dz = int(r['Exact_dZ'])
        if dz > 0:
            if Y[key]['A']['lo'] == Y[key]['A']['hi'] == dz: agree_d += 1
            else: bad.append(('distance', n, dz, Y[key]['A']))
total = sum(len(v) for v in magma.values())
print(f"{len(mine)} lengths, {total} supports compared; {agree_d} exact distances agree")
print("ALL AGREE" if not bad else f"DISAGREEMENTS: {bad}")
sys.exit(1 if bad else 0)
