#!/usr/bin/env python3
"""Exact distillation performance of a triorthogonal code under twirled depolarizing noise.

Inputs rho = (1-e)|M><M| + e/(p-1) sum_{j!=0} Z^j|M><M|Z^-j. For the standard protocol
(measure stabilizers, postselect), an error Z(j) is accepted iff j in A^perp, and is a logical
error iff j not in H^perp. Hence
    P_acc(e) = sum_{j in A^perp} P(j),   e_out(e) = sum_{j in A^perp \\ H^perp} P(j) / P_acc,
with P(j) = (1-e)^{n-w} (e/(p-1))^w. Weight distributions of A^perp and H^perp are obtained
exactly from those of A and H by the MacWilliams identity.
"""
import itertools, sys, os
from fractions import Fraction
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))


def weight_dist(G, p, n):
    G = np.array(G, dtype=np.int64).reshape(-1, n) % p
    k = len(G)
    A = [0] * (n + 1)
    if k == 0:
        A[0] = 1
        return A
    coeffs = np.array(list(itertools.product(range(p), repeat=k)), dtype=np.int64)
    for chunk in range(0, len(coeffs), 200000):
        V = (coeffs[chunk:chunk + 200000] @ G) % p
        w = np.count_nonzero(V, axis=1)
        for x, c in zip(*np.unique(w, return_counts=True)):
            A[int(x)] += int(c)
    return A


def macwilliams(A, p, n):
    """Weight distribution of the dual, from W_{C^perp}(x,y) = W_C(x+(p-1)y, x-y)/|C|."""
    size = sum(A)
    B = [Fraction(0)] * (n + 1)
    from math import comb
    for w, a in enumerate(A):
        if not a:
            continue
        # (x+(p-1)y)^{n-w} (x-y)^w : coefficient of y^t
        for i in range(n - w + 1):
            ci = comb(n - w, i) * (p - 1) ** i
            for j in range(w + 1):
                B[i + j] += a * ci * comb(w, j) * (-1) ** j
    return [b / size for b in B]


def performance(H0, H1, p):
    n = np.array(H1).shape[-1]
    H = np.vstack([np.array(H1).reshape(-1, n), np.array(H0).reshape(-1, n)])
    Aperp = macwilliams(weight_dist(H0, p, n), p, n)
    Hperp = macwilliams(weight_dist(H, p, n), p, n)
    bad = [a - h for a, h in zip(Aperp, Hperp)]
    def stats(e):
        q = e / (p - 1)
        Pacc = sum(float(a) * (1 - e) ** (n - w) * q ** w for w, a in enumerate(Aperp))
        Pbad = sum(float(b) * (1 - e) ** (n - w) * q ** w for w, b in enumerate(bad))
        return Pacc, Pbad / Pacc
    dZ = next(w for w, b in enumerate(bad) if b)
    lead = float(bad[dZ]) / (p - 1) ** dZ
    lo, hi = 1e-6, (p - 1) / p - 1e-9
    if stats(hi)[1] < hi:
        thr = hi
    else:
        f = lambda e: stats(e)[1] - e
        # threshold: first root of e_out(e) = e above 0
        es = np.linspace(1e-4, hi, 4000)
        root = None
        for a, b in zip(es[:-1], es[1:]):
            if f(a) < 0 <= f(b):
                for _ in range(80):
                    m = (a + b) / 2
                    if f(m) < 0: a = m
                    else: b = m
                root = (a + b) / 2
                break
        thr = root
    return dict(n=n, k=len(np.array(H1).reshape(-1, n)), dZ=dZ, lead=lead, threshold=thr,
                Pacc_lin=-float(sum(Aperp[1:2])), stats=stats)
