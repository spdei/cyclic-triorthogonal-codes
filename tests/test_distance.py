"""C distance tools vs. brute force, and the distance formulas of Propositions A and B vs.
the definition d_Z = min wt(A^perp minus H^perp)."""
from _common import *
from triortho.distance import dwit, bdist, bruteforce_x0, dz_definition, distance_cyclic
from triortho.cyclic import code_basis
from triortho.linalg import rref

def test_dwit_matches_bruteforce():
    rng = np.random.default_rng(5)
    for p, N in ((2, 14), (3, 10), (5, 7), (7, 6)):
        for _ in range(25):
            K = int(rng.integers(1, N - 1))
            M = rng.integers(0, p, (K, N))
            targets = [[0] * K] if rng.random() < 0.5 else [list(rng.integers(0, p, K))]
            ref = bruteforce_x0(M, p, targets)
            x = dwit(M, p, targets, N)
            if ref is None:
                assert x is None
            else:
                assert x is not None and int(np.count_nonzero(x)) == ref
                assert x[0] != 0 and tuple((M @ x) % p) in {tuple(t) for t in targets}

def test_bdist_matches_bruteforce():
    rng = np.random.default_rng(6)
    for _ in range(40):
        N = int(rng.integers(6, 15)); K = int(rng.integers(1, N - 2))
        M = rng.integers(0, 2, (K, N))
        for mode in (0, 1):
            MM = np.vstack([M, np.ones(N, dtype=np.int64)]) if mode == 0 else M
            tgt = [[0] * K + [1]] if mode == 0 else [[0] * K]
            ref = bruteforce_x0(MM, 2, tgt)
            L, U, x = bdist(MM, mode, wmaxr=N, trials=50)
            if ref is not None:
                assert L <= ref and (U == -1 or U >= ref)
                if U > 0 and U == L:
                    assert U == ref
                if U == -1:
                    assert L <= ref
            else:
                assert U == -1

def _construction(cons, B, p):
    """Explicit (H0, H1) for construction A ([1;C]) or B (C punctured at coordinate 0)."""
    n = B.shape[1]
    if cons == 'A':
        return B, np.ones((1, n), dtype=np.int64)
    R, piv = rref(B, p)
    assert piv[0] == 0
    return R[1:, 1:], R[:1, 1:]

def test_distance_formulas_match_definition():
    """Every maximal admissible support, plus random sub-unions of cosets (admissibility is
    closed under subsets), for small lengths: certified distance == brute-force definition."""
    from triortho.search import maximal_supports
    rng = random.Random(9)
    count = 0
    for p, ns in ((2, [7, 15, 21, 31]), (3, [8, 10, 13, 16, 20]), (5, [8, 12, 13, 16]), (7, [8, 12, 16])):
        for n in ns:
            cs = cosets(p, n)
            for cons in 'AB':
                orbs, _ = maximal_supports(n, cs, cons)
                cases = []
                for _, S in orbs:
                    cases.append(S)
                    parts = [set(c) for c, _ in cs if set(c) <= S]
                    for _ in range(3):
                        T = set().union(*[c for c in parts if rng.random() < 0.6]) if parts else set()
                        if T:
                            cases.append(T)
                for S in cases:
                    B = code_basis(n, cs, S, p)
                    H0, H1 = _construction(cons, B, p)
                    if p ** (H1.shape[1] - len(H0)) > 3e5:
                        continue
                    ref = dz_definition(H0, H1, p)
                    d, st = distance_cyclic(cons, B, p, wmax=n)
                    assert st == 'exact' and d == ref, (p, n, cons, sorted(S), d, ref)
                    count += 1
    assert count >= 40, count
