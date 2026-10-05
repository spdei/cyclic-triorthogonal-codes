#!/usr/bin/env python3
"""Reproduce Tables 1 and 2 of the paper and the other numbers quoted in its text.

usage: python scripts/paper_tables.py
"""
import itertools, os, sys
from math import gcd, log
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, '..')
sys.path.insert(0, ROOT)
sys.path.insert(0, HERE)
from triortho.cyclic import load_cosets, code_basis
from triortho.gf import cyclotomic_data
from triortho.linalg import rref
from threshold import performance
import make_tables as mt

COSETS = {}
def cosets(p, n):
    if p not in COSETS:
        COSETS[p] = load_cosets(os.path.join(ROOT, 'data', f'cosets_p{p}.txt'))
    return COSETS[p][(p, n)]


def code_A(p, n, S, cos=None):
    """Construction A: logical row 1, stabilizers the cyclic code with support S."""
    return code_basis(n, cos or cosets(p, n), S, p), np.ones((1, n), dtype=np.int64)


def code_B(p, N, S, cos=None):
    """Construction B: the cyclic space with support S, punctured at coordinate 0."""
    R, _ = rref(code_basis(N, cos or cosets(p, N), S, p), p)
    return R[1:, 1:], R[:1, 1:]


# ------------------------------------------------------------------------------------- Table 1
def rm_codes(p, nmax=3000):
    """Punctured GRM_p(r,m), largest admissible r (3r < m(p-1)): list of (length, distance)."""
    out, m = [], 1
    while p ** m - 1 <= nmax:
        r = max((t for t in range(1, m * (p - 1)) if 3 * t < m * (p - 1)), default=0)
        if r:
            a, b = divmod(m * (p - 1) - r - 1, p - 1)        # dual GRM_p(rho, m), rho = a(p-1)+b
            out.append((p ** m - 1, (p - b) * p ** (m - a - 1) - 1))
        m += 1
    return out


def cyclic_records(p, binary, qudit):
    """{n_quantum: (best certified lower bound, constructions)} over the search."""
    rec = {}
    for r in qudit:
        if r['p'] == p and r['codes']:
            n, d = r['codes'][0]['n'], mt.best(r['codes'])[0]
            old = rec.get(n, (0, ''))
            rec[n] = (d, r['cons']) if d > old[0] else (d, old[1] + ', ' + r['cons']) if d == old[0] else old
    return rec


def table1(binary, qudit):
    print("Table 1: shortest length of an [[n,1,d']]_p code with d' >= d")
    print("      " + "".join(f"{'p=%d' % p:>27}" for p in (3, 5, 7)))
    print("  d " + "      RM  conc.   cyclic    " * 3)
    cols = {}
    for p in (3, 5, 7):
        rm = rm_codes(p)
        cyc = cyclic_records(p, binary, qudit)
        cyc_col = {d: min(((n, c) for n, (dd, c) in cyc.items() if dd >= d), default=None) for d in range(2, 9)}
        base = set(rm) | {(n, cyc[n][0]) for n, _ in filter(None, cyc_col.values())}
        conc = {}
        for k in (2, 3):
            for combo in itertools.combinations_with_replacement(sorted(base), k):
                n = int(np.prod([c[0] for c in combo])); d = int(np.prod([c[1] for c in combo]))
                for dd in range(2, 9):
                    if d >= dd and n < conc.get(dd, 10 ** 9):
                        conc[dd] = n
        cols[p] = {d: (min((n for n, dist in rm if dist >= d), default=None), conc.get(d), cyc_col[d])
                   for d in range(2, 9)}
    for d in range(2, 9):
        row = f"  {d} "
        for p in (3, 5, 7):
            r, c, y = cols[p][d]
            row += f"{r if r else '-':>8}{c if c else '-':>7}" + (f"{y[0]:>6} ({y[1]})" if y else "     -").ljust(14)
        print(row)
    print()


# ------------------------------------------------------------------------------------- Table 2
def maximal_support_containing(p, n, S, binary, qudit):
    units = [u for u in range(1, n) if gcd(u, n) == 1]
    recs = binary[n]['codes'] if p == 2 else next(r for r in qudit if r['p'] == p and r['cons'] == 'A' and r['N'] == n)['codes']
    for c in recs:
        if any({(u * s) % n for s in S} <= set(c['S']) for u in units):
            return c['S']


def digit_sum(j, p):
    s = 0
    while j:
        s, j = s + j % p, j // p
    return s


def rm_support(p, m, r):
    n = p ** m - 1
    return {j for j in range(1, n) if digit_sum(j, p) <= r}


def show(label, H0, H1, p):
    r = performance(H0, H1, p)
    n, d = H1.shape[1], r['dZ']
    print(f"  {label:<26} gamma={log(n) / log(d):5.2f}  eps_out ~ {r['lead']:.4g} eps^{d}"
          f"  threshold={r['threshold']:.3f}  P_acc(0.01)={r['stats'](0.01)[0]:.2f}")


def table2(binary, qudit):
    print("Table 2: single-round performance under twirled depolarizing noise")
    rows = [
        ("[[19,1,3]]_3 (B)", 3, code_B(3, 20, {2, 5, 6, 14, 15, 18})),
        ("[[40,1,4]]_3 (A)", 3, code_A(3, 40, {1, 2, 3, 5, 6, 9, 14, 15, 18, 27})),
        ("[[8,1,2]]_3 (RM)", 3, code_A(3, 8, rm_support(3, 2, 1))),
        ("[[26,1,2]]_3 (RM)", 3, code_A(3, 26, rm_support(3, 3, 1))),
        ("[[26,1,2]]_3 (A)", 3, code_A(3, 26, maximal_support_containing(3, 26, rm_support(3, 3, 1), binary, qudit))),
        ("[[11,1,3]]_5 (B)", 5, code_B(5, 12, {2, 3, 9, 10})),
        ("[[24,1,3]]_5 (RM)", 5, code_A(5, 24, rm_support(5, 2, 2))),
        ("[[4,1,2]]_5 (RS)", 5, code_A(5, 4, {1})),
        ("[[16,1,4]]_5 (RS o RS)", 5, code_A(5, 16, {1, 4, 5, 9, 13})),
        ("[[7,1,3]]_7 (B)", 7, code_B(7, 8, {1, 4, 7})),
        ("[[12,1,3]]_7 (A)", 7, code_A(7, 12, {1, 2, 7})),
        ("[[6,1,2]]_7 (RS)", 7, code_A(7, 6, {1})),
    ]
    for label, p, (H0, H1) in rows:
        show(label, H0, H1, p)
    print()


# ------------------------------------------------------------------------------ numbers in text
def text_numbers(binary, qudit):
    print("Numbers quoted in the text")
    rm = rm_support(2, 6, 1)                                     # Section 7, n = 63
    show("RM(1,6)* [[63,1,3]]_2", *code_A(2, 63, rm), 2)
    show("cyclic [[63,1,3]]_2", *code_A(2, 63, maximal_support_containing(2, 63, rm, binary, qudit)), 2)
    p = 13                                                       # Section 8, p = 13
    cos = cyclotomic_data(14, p, seed=1)
    show("[[13,1,5]]_13 (Prop. 14)", *code_B(p, 14, set(range(5, 10)), cos), p)
    show("[[12,1,4]]_13 (RS)", *code_A(p, 12, {1, 2, 3}, cyclotomic_data(12, p, seed=1)), p)
    print()


def main():
    binary, qudit = mt.load()
    mt.apply_concatenation(binary, qudit)
    table1(binary, qudit)
    table2(binary, qudit)
    text_numbers(binary, qudit)


if __name__ == '__main__':
    main()
