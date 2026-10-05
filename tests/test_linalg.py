from _common import *
from triortho.linalg import rref, rank, nullspace, in_rowspan

def test_rref_and_nullspace():
    rng = np.random.default_rng(1)
    for p in (2, 3, 5, 7):
        for _ in range(60):
            m, n = rng.integers(1, 9), rng.integers(1, 12)
            M = rng.integers(0, p, (m, n))
            if rng.random() < 0.3:                     # force dependencies
                M = np.vstack([M, (M[0] * 2 + M[-1]) % p])
            R, piv = rref(M, p)
            r = len(piv)
            assert R.shape[0] == r
            for i, c in enumerate(piv):                 # reduced echelon structure
                assert R[i, c] == 1 and np.count_nonzero(R[:, c]) == 1
            for row in M:                               # same row space
                assert in_rowspan(row, R, p)
            for row in R:
                assert in_rowspan(row, M, p)
            N = nullspace(M, p)
            assert N.shape[0] == n - r
            assert not ((M @ N.T) % p).any()
            assert rank(N, p) == n - r if len(N) else True
