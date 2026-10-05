import os, sys, random
import numpy as np
ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
sys.path.insert(0, ROOT)
from triortho.cyclic import load_cosets

_DATA = {}
def cosets(p, n):
    if p not in _DATA:
        _DATA[p] = load_cosets(os.path.join(ROOT, 'data', f'cosets_p{p}.txt'))
    return _DATA[p][(p, n)]

def lengths(p, nmax):
    if p not in _DATA:
        _DATA[p] = load_cosets(os.path.join(ROOT, 'data', f'cosets_p{p}.txt'))
    return sorted(n for (q, n) in _DATA[p] if q == p and n <= nmax)

def random_support(p, n, rng, include_zero=True):
    cs = [c for c, _ in cosets(p, n) if include_zero or c != (0,)]
    S = set()
    for c in cs:
        if rng.random() < 0.5:
            S |= set(c)
    return S
