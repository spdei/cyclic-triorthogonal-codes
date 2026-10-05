"""Exhaustive transversal-gate checks (exact phase arithmetic mod 8, 9, p) on the headline
codes, and random tests of the transversality theorems: tri_check <=> transversal gate acts
as a non-Clifford tensor-product logical gate up to a Pauli."""
from _common import *
from triortho.cyclic import code_basis
from triortho.checks import tri_check, transversal_check
from triortho.linalg import rref

def punctured(p, N, S):
    B = code_basis(N, cosets(p, N), S, p)
    R, piv = rref(B, p)
    return R[:, 1:]                         # first row logical (R[0,0]=1), others stabilizers

def test_qutrit_19():
    H = punctured(3, 20, {2, 5, 6, 14, 15, 18})
    assert tri_check(H, 3, 1)
    r = transversal_check(H, 1, 3, gate='U')
    assert r['status'] in ('exact', 'pauli') and [r['logical'][(x,)] for x in range(3)] == [0, 8, 1]
    assert transversal_check(H, 1, 3, gate='V')['status'] == 'fail'      # not self-orthogonal

def test_trivial_carry_example():
    H = np.array([[2]])
    rV = transversal_check(H, 1, 3, 'V'); rU = transversal_check(H, 1, 3, 'U')
    assert [rV['logical'][(x,)] for x in range(3)] == [0, 2, 1]      # not V^eps = (0,2,4)
    assert [rU['logical'][(x,)] for x in range(3)] == [0, 8, 1]

def test_qudit_headline_codes():
    for p, N, S in ((5, 12, {2, 3, 9, 10}), (7, 8, {1, 4, 7})):
        H = punctured(p, N, S)
        assert tri_check(H, p, 1)
        r = transversal_check(H, 1, p)
        assert r['status'] == 'exact'                                # no correction for p >= 5
        assert all(r['logical'][(x,)] == (-x ** 3) % p for x in range(p))   # logical U^{-1}

def test_qubit_15():
    n = 15
    B = code_basis(n, cosets(2, n), {1, 2, 4, 8}, 2)
    H = np.vstack([np.ones(n, dtype=np.int64), B])
    assert tri_check(H, 2, 1)
    r = transversal_check(H, 1, 2)
    assert r['status'] == 'exact' and r['logical'][(1,)] == n % 8     # triply even: exact T^n

def test_theorem_random_small_matrices():
    """For random small H over F_3 and F_5 (k = 1): tri_check(H) holds iff the transversal
    cubic gate preserves the code space up to a Pauli and acts as a non-Clifford gate."""
    rng = np.random.default_rng(13)
    hits = {True: 0, False: 0}
    for p, n in ((3, 6), (3, 7), (5, 5)):
        for _ in range(300):
            m = int(rng.integers(1, 4))
            H = rng.integers(0, p, (m, n))
            if len(rref(H, p)[1]) < m:
                continue
            ok = tri_check(H, p, 1)
            r = transversal_check(H, 1, p, 'U')
            good = r['status'] != 'fail'
            if good:
                L = [r['logical'][(x,)] for x in range(p)]
                q = 9 if p == 3 else p
                # non-Clifford: cubic part nonzero, i.e. lambda != 0
                if p == 3:
                    good = (L[1] % 3) != 0
                else:
                    good = any((L[x] - L[1] * x ** 3) % p for x in range(p)) is False and L[1] % p != 0
            assert ok == good, (p, H.tolist(), r)
            hits[ok] += 1
    assert hits[True] > 3 and hits[False] > 50

def test_binary_triorthogonal_needs_clifford_correction():
    """Bravyi-Haah triorthogonality (p = 2) only guarantees transversal T up to Clifford
    corrections: a weight-2 stabilizer gives phase dependence 2a (mod 8), an S-type phase,
    which is not a Pauli. The checker must report this as 'fail' (not Pauli-correctable),
    while a triply-even code is exact."""
    H = np.array([[1, 1, 1], [0, 1, 1]])                     # phase 2a, plus a cross term in x
    assert tri_check(H, 2, 1)
    assert transversal_check(H, 1, 2)['status'] == 'fail'
    Hd = np.array([[1, 0, 0], [0, 1, 1]])                    # disjoint: phase exactly 2a, x-independent
    assert tri_check(Hd, 2, 1)
    assert transversal_check(Hd, 1, 2)['status'] == 'fail'   # linear but an S-type phase, not Pauli
    H4 = np.array([[1, 1, 1, 1, 1], [0, 1, 1, 1, 1]])          # weight-4 stabilizer: Pauli (Z) phase
    r = transversal_check(H4, 1, 2)
    assert tri_check(H4, 2, 1) and r['status'] == 'pauli'

def test_V_gate_needs_self_orthogonality():
    """For p = 3, Z(w) V^{otimes n} acts diagonally on the code space iff every triple of rows
    involving a stabilizer row has vanishing trilinear form AND the stabilizer space is
    orthogonal to the whole row space (condition L). Checked on random small matrices."""
    from triortho.checks import generalised_check
    rng = np.random.default_rng(23)
    hits = {True: 0, False: 0}
    for n in (3, 4, 5, 6):
        for _ in range(1500):
            k = int(rng.integers(1, 3)); r = int(rng.integers(1, 3))
            H = rng.integers(0, 3, (k + r, n))
            if len(rref(H, 3)[1]) < k + r:
                continue
            A, Hall = H[k:], H
            L = not ((A @ Hall.T) % 3).any()
            expect = generalised_check(H, 3, k) and L
            got = transversal_check(H, k, 3, 'V')['status'] != 'fail'
            assert expect == got, H.tolist()
            hits[expect] += 1
    assert hits[True] >= 5 and hits[False] > 100, hits
