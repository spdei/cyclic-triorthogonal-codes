#!/usr/bin/env python3
"""Classify cyclic [[n,1,d]]_p triorthogonal codes (p odd), constructions A and B.

usage: python scripts/search_odd.py P NMAX [--out results_pP.jsonl]

For each length n in data/cosets_pP.txt with n <= NMAX: enumerate all maximal admissible
spectral supports up to multipliers, re-check triorthogonality directly, and compute a
certified distance (exact with verified witness, or a rigorous lower bound).
"""
import argparse, json, os, sys, time
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from triortho.cyclic import load_cosets, code_basis
from triortho.search import maximal_supports
from triortho.checks import tri_check
from triortho.distance import distance_cyclic


def run_length(n, cosets, p):
    out = {}
    for cons in 'AB':
        t0 = time.time()
        orbs, complete = maximal_supports(n, cosets, cons)
        recs = []
        for key, S in orbs:
            B = code_basis(n, cosets, S, p)
            if cons == 'A':
                assert tri_check(np.vstack([np.ones(n, dtype=np.int64), B]), p, 1)
                length = n
            else:
                assert tri_check(B, p, 0)
                length = n - 1
            d, st = distance_cyclic(cons, B, p)
            recs.append(dict(S=list(key), dim=int(len(B)), n=length, d=d, status=st))
        best = max(recs, key=lambda r: (r['d'], r['status'] == 'exact')) if recs else None
        out[cons] = dict(p=p, cons=cons, N=n, complete=complete, n_orbits=len(recs),
                         codes=recs, best=best, secs=round(time.time() - t0, 1))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('p', type=int)
    ap.add_argument('nmax', type=int)
    ap.add_argument('--out', default=None)
    a = ap.parse_args()
    root = os.path.join(os.path.dirname(__file__), '..')
    data = load_cosets(os.path.join(root, 'data', f'cosets_p{a.p}.txt'))
    outpath = a.out or f'results_p{a.p}.jsonl'
    with open(outpath, 'w') as f:
        for (pp, n), cosets in sorted(data.items()):
            if pp != a.p or n > a.nmax:
                continue
            res = run_length(n, cosets, a.p)
            for cons in 'AB':
                r = res[cons]
                f.write(json.dumps(r) + '\n')
                b = r['best']
                if b:
                    print(a.p, cons, n, 'orbits', r['n_orbits'],
                          f"[[{b['n']},1,{b['d']}{'' if b['status'] == 'exact' else '+'}]]", flush=True)


if __name__ == '__main__':
    main()
