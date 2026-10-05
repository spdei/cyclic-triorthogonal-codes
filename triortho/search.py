"""Enumeration of maximal admissible spectral supports.

Construction A (cyclic stabilizers, all-ones logical row), admissible iff
    0 not in S, 0 not in 2S, 0 not in 3S.
Construction B (punctured cyclic triorthogonal space), admissible iff 0 not in 3S.
For p = 2 the two conditions coincide (S is closed under doubling, so S c 2S c 3S).
Both are closed under taking subsets of S, so a depth-first search over unions of cosets,
pruned by the condition, enumerates every admissible support; we keep the maximal ones.
"""
import time
from math import gcd


class _Bits:
    """Subsets of Z_n as Python ints, with fast sumsets via cyclic rotation."""
    def __init__(self, n):
        self.n, self.full = n, (1 << n) - 1
    def rot(self, x, a):
        n = self.n
        return ((x << a) | (x >> (n - a))) & self.full if a else x
    def elems(self, x):
        out, i = [], 0
        while x:
            if x & 1:
                out.append(i)
            x >>= 1
            i += 1
        return out
    def neg(self, x):
        return sum(1 << ((-e) % self.n) for e in self.elems(x))
    def double(self, x):          # 2S as a bitmask
        two = 0
        for a in self.elems(x):
            two |= self.rot(x, a)
        return two


def admissible_mask(Sm, bits, cons):
    if Sm & 1 and cons == 'A':                       # 0 in S
        return False
    two = bits.double(Sm)
    negS = bits.neg(Sm)
    if two & negS:                                   # 0 in 3S  <=>  -c in 2S for some c in S
        return False
    if cons == 'A' and (two & 1):                    # 0 in 2S
        return False
    return True


def admissible(S, n, cons):
    bits = _Bits(n)
    return admissible_mask(sum(1 << s for s in S), bits, cons)


def maximal_supports(n, cosets, cons, time_limit=None):
    """All maximal admissible supports for construction cons in {'A','B'}, up to multipliers.
    Returns (list of (canonical tuple, set S)), complete_flag)."""
    bits = _Bits(n)
    cs = [c for c, _ in cosets if not (cons == 'A' and c == (0,))]
    cm = [sum(1 << e for e in c) for c in cs]
    t0 = time.time()
    res, complete = [], [True]
    def dfs(i, Sm):
        if time_limit is not None and time.time() - t0 > time_limit:
            complete[0] = False
            return
        if i == len(cs):
            if Sm and all((Sm | m) == Sm or not admissible_mask(Sm | m, bits, cons) for m in cm):
                res.append(Sm)
            return
        T = Sm | cm[i]
        if admissible_mask(T, bits, cons):
            dfs(i + 1, T)
        dfs(i + 1, Sm)
    dfs(0, 0)
    units = [u for u in range(1, n) if gcd(u, n) == 1]
    seen = {}
    for Sm in res:
        E = bits.elems(Sm)
        key = min(tuple(sorted((u * e) % n for e in E)) for u in units)
        seen.setdefault(key, set(E))
    return sorted(seen.items()), complete[0]
