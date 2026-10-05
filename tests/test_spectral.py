"""Lemma (spectral support): Schur powers have support tS, duals have support Z_n \\ (-S),
coordinate sums vanish iff 0 not in S. Checked numerically on random supports."""
from _common import *
from triortho.cyclic import code_basis, sumset
from triortho.linalg import rank, nullspace

def schur_span(B, t, p):
    rows = [B[i] for i in range(len(B))]
    if t == 2:
        P = [(B[i] * B[j]) % p for i in range(len(B)) for j in range(i, len(B))]
    else:
        P = [(B[i] * B[j] * B[k]) % p for i in range(len(B)) for j in range(i, len(B)) for k in range(j, len(B))]
    return np.array(P, dtype=np.int64).reshape(-1, B.shape[1])

def test_spectral_lemma():
    rng = random.Random(7)
    for p, ns in ((2, [7, 15, 21, 31, 45]), (3, [8, 10, 13, 20, 26]), (5, [8, 12, 13, 24]), (7, [8, 12, 16, 24])):
        for n in ns:
            for _ in range(6):
                S = random_support(p, n, rng)
                if not S:
                    continue
                B = code_basis(n, cosets(p, n), S, p)
                for t in (2, 3):
                    assert rank(schur_span(B, t, p), p) == len(sumset(S, n, t)), (p, n, t)
                D = nullspace(B, p)
                Sd = set(range(n)) - {(-s) % n for s in S}
                B2 = code_basis(n, cosets(p, n), Sd, p) if Sd else np.zeros((0, n), dtype=np.int64)
                assert len(D) == len(B2) and rank(np.vstack([D, B2]), p) == len(D)
                assert (not (B.sum(1) % p).any()) == (0 not in S)
