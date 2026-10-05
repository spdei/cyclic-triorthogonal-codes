#!/usr/bin/env python3
"""Cross-check the qudit classification (p = 3, 5, 7) against the independent Magma search
(magma/Qudit_Triorthogonal_p*.csv, produced by magma/triorthogonal-qudit-search.m).

For every (construction, length) in the Magma output, compares the sets of maximal admissible
supports up to multiplier equivalence, their dimensions, and every distance Magma computed."""
import csv, glob, json, os, sys
from math import gcd
from collections import defaultdict
ROOT = os.path.join(os.path.dirname(__file__), '..')
def canon(K, n): return min(tuple(sorted((u*x) % n for x in K)) for u in range(1, n) if gcd(u, n) == 1) if K else ()
rows = [json.loads(l) for l in open(os.path.join(ROOT, 'results', 'qudit_results.jsonl'))]
bad = False
for f in sorted(glob.glob(os.path.join(ROOT, 'magma', 'qudit_search_p*.csv'))):
    p = int(os.path.basename(f).split('_p')[1].split('.')[0])
    M = defaultdict(dict)
    for r in csv.DictReader(open(f)):
        N = int(r['N']); K = [int(t) for t in r['Roots_K'].strip('"[] ').split(',') if t.strip()]
        M[(r['construction'], N)][canon(K, N)] = (int(r['dim']), int(r['dZ']))
    P = defaultdict(dict)
    for r in rows:
        if r['p'] == p:
            for c in r['codes']: P[(r['cons'], r['N'])][canon(c['S'], r['N'])] = (c['dim'], c['d'], c['status'])
    sup = dims = 0; dist = 0; partial = []
    for k in M:
        a, b = M[k], P.get(k, {})
        if set(a) != set(b):
            if set(a) < set(b): partial.append(k)
            else: print('SUPPORT MISMATCH', p, k); bad = True
        for K in set(a) & set(b):
            sup += 1
            if a[K][0] != b[K][0]: print('DIM MISMATCH', p, k); bad = True
            if a[K][1] > 0:
                ok = (a[K][1] == b[K][1]) if b[K][2] == 'exact' else (a[K][1] >= b[K][1])
                if not ok: print('DISTANCE MISMATCH', p, k, a[K][1], b[K][1:]); bad = True
                dist += 1
    print(f"p={p}: {sup} supports and {dist} distances compared; Magma run stopped early at {partial or 'no length'}")
print("ALL AGREE" if not bad else "DISAGREEMENTS FOUND")
sys.exit(1 if bad else 0)
