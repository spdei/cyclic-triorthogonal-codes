"""Distances of cyclic triorthogonal codes.

Certified computation (C tools in csrc/):
  dwit  (any p): min weight of x with x_0 = 1 and sum_i x_i b_i = target, by meet-in-the-middle
                 with exact reconstruction of every hash match, so every reported weight comes
                 with a witness; the enumeration is exhaustive, so absence below is rigorous.
  bdist (p = 2): the same, bit-packed, plus a randomized split phase for witnesses above the cap.

Brute-force references (for tests): bruteforce_x0 and dz_definition.
"""
import itertools, os, subprocess
from math import comb
import numpy as np
from .linalg import nullspace

HERE = os.path.dirname(os.path.abspath(__file__))
BUILD = os.path.join(os.path.dirname(HERE), 'build')


def _binary(name):
    path = os.path.join(BUILD, name)
    if not os.path.exists(path):
        os.makedirs(BUILD, exist_ok=True)
        src = os.path.join(os.path.dirname(HERE), 'csrc', name + '.c')
        subprocess.run(['gcc', '-O2', '-o', path, src], check=True, capture_output=True)
    return path


def dwit(M, p, targets, wmax):
    """Witness x (numpy vector) of minimal weight <= wmax, or None."""
    M = np.array(M, dtype=np.int64) % p
    K, N = M.shape
    inp = f"{p} {K} {N} {wmax} {len(targets)}\n" + "\n".join(" ".join(map(str, r)) for r in M) \
        + "\n" + "\n".join(" ".join(map(str, t)) for t in targets)
    out = subprocess.run([_binary('dwit')], input=inp, capture_output=True, text=True).stdout.split()
    if out[0] == '-1':
        return None
    x = np.zeros(N, dtype=np.int64)
    for tok in out[1:]:
        i, c = map(int, tok.split(':'))
        x[i] = c
    return x


def choose_wmax(N, p, nt, cap_tab=1.2e7, cap_look=6e7, hard=12):
    w = 1
    while w < hard:
        h1 = w // 2
        h2 = w - h1
        if comb(N - 1, h1) * (p - 1) ** h1 > cap_tab or nt * comb(N - 1, h2) * (p - 1) ** h2 > cap_look:
            break
        w += 1
    return w


def distance_cyclic(cons, B, p, wmax=None):
    """Distance of construction A (H = [1; B], length n) or B (B punctured at 0, length n-1)
    for a cyclic code with generator matrix B. Returns (d, status) with status 'exact' or
    'lower' (rigorous lower bound)."""
    B = np.array(B, dtype=np.int64) % p
    n = B.shape[1]
    if cons == 'B':
        targets = [[0] * len(B)]
        w = choose_wmax(n, p, 1) if wmax is None else wmax
        x = dwit(B, p, targets, w)
        if x is None:
            return w, 'lower'                    # d(W^perp) > w  =>  d >= w
        assert x[0] != 0 and not ((B @ x) % p).any()
        return int(np.count_nonzero(x)) - 1, 'exact'
    M = np.vstack([B, np.ones(n, dtype=np.int64)])
    targets = [[0] * len(B) + [t] for t in range(1, p)]
    w = choose_wmax(n, p, p - 1) if wmax is None else wmax
    x = dwit(M, p, targets, w)
    if x is None:
        return w + 1, 'lower'
    assert x[0] != 0 and not ((B @ x) % p).any() and x.sum() % p
    return int(np.count_nonzero(x)), 'exact'


def bdist(M, mode, captab=1.5e7, caplook=2e8, wmaxr=15, trials=40):
    """Binary search; mode 0 = construction A (last row of M is all-ones), 1 = construction B.
    Returns (L, U, witness): no witness of weight < L exists; U is the weight of a verified
    witness (or -1)."""
    M = np.array(M, dtype=np.int64) % 2
    K, n = M.shape
    inp = f"{n} {K} {mode} {int(captab)} {int(caplook)} {wmaxr} {trials}\n" + "\n".join(" ".join(map(str, r)) for r in M)
    out = subprocess.run([_binary('bdist')], input=inp, capture_output=True, text=True).stdout.split()
    L, U = int(out[0]), int(out[1])
    supp = list(map(int, out[2:]))
    x = None
    if U > 0:
        x = np.zeros(n, dtype=np.int64)
        x[supp] = 1
        tgt = np.zeros(K, dtype=np.int64)
        if mode == 0:
            tgt[-1] = 1
        assert ((M @ x) % 2 == tgt).all() and x[0] == 1 and len(supp) == U
    if mode == 0 and L % 2 == 0:
        L += 1                                    # only odd weights are possible in A
    return L, U, x


# ---------------------------------------------------------------- brute-force references
def bruteforce_x0(M, p, targets, wmax=None):
    """Minimal weight of x with x_0 = 1 and M x in targets, by exhaustive enumeration."""
    M = np.array(M, dtype=np.int64) % p
    K, N = M.shape
    T = {tuple(np.array(t) % p) for t in targets}
    best = None
    rest = np.array(list(itertools.product(range(p), repeat=N - 1)), dtype=np.int64).reshape(-1, N - 1)
    X = np.hstack([np.ones((len(rest), 1), dtype=np.int64), rest])
    S = (X @ M.T) % p
    for x, s in zip(X, S):
        if tuple(s) in T:
            w = int(np.count_nonzero(x))
            if best is None or w < best:
                best = w
    return best


def dz_definition(H0, H1, p):
    """d_Z = min wt(v), v in A^perp minus H^perp, straight from the definition."""
    H0 = np.array(H0, dtype=np.int64).reshape(-1, np.array(H1).shape[-1]) % p
    H1 = np.array(H1, dtype=np.int64).reshape(-1, H0.shape[1]) % p
    n = H1.shape[1]
    N = nullspace(H0, p) if len(H0) else np.eye(n, dtype=np.int64)
    best = None
    for c in itertools.product(range(p), repeat=len(N)):
        v = (np.array(c, dtype=np.int64) @ N) % p
        if ((H1 @ v) % p).any():
            w = int(np.count_nonzero(v))
            if best is None or w < best:
                best = w
    return best


def dz_general(H0, H1, p, wmax):
    """Certified d_Z = min wt(v): H0 v = 0, H1 v != 0, for a code without assumed symmetry:
    minimum over the first support coordinate j and the logical row i with h^i.v != 0."""
    n = H1.shape[1]
    best, wit = None, None
    for i in range(len(H1)):
        M = np.vstack([H0, H1[i:i + 1]]) if len(H0) else H1[i:i + 1]
        targets = [[0] * len(H0) + [t] for t in range(1, p)]
        for j in range(n):
            order = [j] + [l for l in range(n) if l != j]
            cap = wmax if best is None else min(wmax, best - 1)
            if cap < 1:
                break
            x = dwit(M[:, order], p, targets, cap)
            if x is not None:
                v = np.zeros(n, dtype=np.int64); v[order] = x
                assert not ((H0 @ v) % p).any() and (H1 @ v % p).any()
                w = int(np.count_nonzero(v))
                if best is None or w < best:
                    best, wit = w, v
    if best is None:
        return wmax + 1, 'lower', None
    return best, 'exact', wit
