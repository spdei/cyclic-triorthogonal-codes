"""Regression: reproduce headline results end to end (search + certified distance)."""
from _common import *
from triortho.cyclic import code_basis
from triortho.search import maximal_supports
from triortho.distance import distance_cyclic, bdist

def best(p, n, cons):
    orbs, complete = maximal_supports(n, cosets(p, n), cons)
    assert complete
    return max((distance_cyclic(cons, code_basis(n, cosets(p, n), S, p), p)[0] for _, S in orbs), default=0)

def test_odd_prime_records():
    assert best(3, 20, 'B') == 3                       # [[19,1,3]]_3
    assert all(best(3, n, 'B') < 3 for n in lengths(3, 19))   # nothing smaller with d >= 3
    assert best(3, 8, 'A') == 2                        # [[8,1,2]]_3  (Campbell-Anwar-Browne)
    assert best(3, 40, 'A') == 4                       # [[40,1,4]]_3
    assert best(5, 12, 'B') == 3                       # [[11,1,3]]_5
    assert best(5, 16, 'A') == 4                       # [[16,1,4]]_5
    assert best(7, 8, 'B') == 3                        # [[7,1,3]]_7
    assert best(7, 29, 'B') == 5                       # [[28,1,5]]_7
    assert best(7, 6, 'A') == 2                        # [[6,1,2]]_7 (Reed-Solomon)

def test_qubit_records():
    for n, d in ((15, 3), (85, 5), (127, 7), (223, 9)):
        orbs, complete = maximal_supports(n, cosets(2, n), 'B')
        assert complete
        dists = []
        for _, S in orbs:
            B = code_basis(n, cosets(2, n), S, 2)
            L, U, _ = bdist(np.vstack([B, np.ones(n, dtype=np.int64)]), 0)
            dists.append((L, U))
        assert max(L for L, U in dists) == d and any(L == U == d for L, U in dists), (n, dists)

def test_mds_families():
    """Proposition (two MDS families): construction A, n = p-1, S = {1..r}, 3r < p-1 gives
    [[p-1,1,r+1]]_p; construction B, n = p+1, S = interval of radius r around n/2, 6r < p+1,
    gives [[p,1,2r+1]]_p."""
    from triortho.gf import cyclotomic_data
    from triortho.checks import cube_orthogonal, tri_check
    from triortho.linalg import rref
    from triortho.distance import dz_general
    for p in (7, 11, 13):
        n = p + 1; h = n // 2; r = max(t for t in range(p) if 6 * t < n)
        W = code_basis(n, cyclotomic_data(n, p, seed=1), {(h + t) % n for t in range(-r, r + 1)}, p)
        assert cube_orthogonal(W, p)
        R, _ = rref(W, p); H0, H1 = R[1:, 1:], R[:1, 1:]
        assert tri_check(np.vstack([H1, H0]), p, 1)
        assert dz_general(H0, H1, p, wmax=2 * r + 1)[:2] == (2 * r + 1, 'exact')
        m = p - 1; ra = max(t for t in range(1, p) if 3 * t < m)
        C = code_basis(m, cyclotomic_data(m, p, seed=1), set(range(1, ra + 1)), p)
        ones = np.ones((1, m), dtype=np.int64)
        assert tri_check(np.vstack([ones, C]), p, 1)
        assert dz_general(C, ones, p, wmax=ra + 1)[:2] == (ra + 1, 'exact')
