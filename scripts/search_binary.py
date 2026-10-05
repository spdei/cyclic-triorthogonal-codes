#!/usr/bin/env python3
"""Classify binary cyclic triorthogonal codes of odd length (constructions A and B share the
supports with 0 not in 3S). Two-stage distances: a cheap certified lower bound for every
support, then a deep certified search on the best few.

usage: python scripts/search_binary.py NMIN NMAX [--time-limit 60] [--out results_p2.jsonl]
"""
import argparse, json, os, sys, time
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from triortho.cyclic import load_cosets, code_basis
from triortho.search import maximal_supports
from triortho.checks import tri_check
from triortho.distance import bdist


def triply_even(B):
    Bi = np.array(B, dtype=np.int64)
    if (Bi.sum(1) % 8).any() or ((Bi @ Bi.T) % 4).any():
        return False
    return all(not ((((Bi * Bi[i])) @ Bi.T) % 2).any() for i in range(len(Bi)))


def run_length(n, cosets, tlim):
    t0 = time.time()
    orbs, complete = maximal_supports(n, cosets, 'B', time_limit=tlim)
    codes = []
    for key, S in orbs:
        B = code_basis(n, cosets, S, 2)
        assert tri_check(B, 2, 0)
        MA = np.vstack([B, np.ones(n, dtype=np.int64)])
        LA, UA, _ = bdist(MA, 0, captab=4e5, caplook=4e6, wmaxr=0, trials=0)
        LB, UB, _ = bdist(B, 1, captab=4e5, caplook=4e6, wmaxr=0, trials=0)
        codes.append(dict(S=list(key), dim=int(len(B)), triply_even=triply_even(B),
                          A=dict(n=n, lo=LA, hi=UA), B=dict(n=n - 1, lo=LB - 1, hi=(UB - 1 if UB > 0 else -1)),
                          _B=B))
    for k, mode in (('A', 0), ('B', 1)):
        if not codes:
            break
        best = max(c[k]['lo'] for c in codes)
        for c in [c for c in codes if c[k]['hi'] != c[k]['lo'] and c[k]['lo'] >= best][:4]:
            M = np.vstack([c['_B'], np.ones(n, dtype=np.int64)]) if k == 'A' else c['_B']
            L, U, _ = bdist(M, mode, captab=8e6, caplook=5e7, wmaxr=15, trials=8)
            c[k] = dict(n=n, lo=L, hi=U) if k == 'A' else dict(n=n - 1, lo=L - 1, hi=(U - 1 if U > 0 else -1))
    for c in codes:
        c.pop('_B')
    return dict(N=n, complete=complete, n_orbits=len(codes), codes=codes, secs=round(time.time() - t0, 1))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('nmin', type=int)
    ap.add_argument('nmax', type=int)
    ap.add_argument('--time-limit', type=float, default=60)
    ap.add_argument('--out', default='results_p2.jsonl')
    a = ap.parse_args()
    root = os.path.join(os.path.dirname(__file__), '..')
    data = load_cosets(os.path.join(root, 'data', 'cosets_p2.txt'))
    with open(a.out, 'a') as f:
        for (p, n), cosets in sorted(data.items()):
            if not (a.nmin <= n <= a.nmax):
                continue
            r = run_length(n, cosets, a.time_limit)
            f.write(json.dumps(r) + '\n')
            f.flush()
            if r['codes']:
                bA = max(r['codes'], key=lambda c: c['A']['lo'])
                print(n, 'orbits', r['n_orbits'], '' if r['complete'] else '(INCOMPLETE)',
                      f"A: d>={bA['A']['lo']} (witness {bA['A']['hi']})", flush=True)


if __name__ == '__main__':
    main()
