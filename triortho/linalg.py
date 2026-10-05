"""Linear algebra over the prime field F_p, using numpy int64 arrays."""
import numpy as np


def rref(M, p):
    """Reduced row echelon form of M over F_p.

    Returns (R, pivots): R has rank(M) rows, pivot entries equal to 1, and every pivot
    column is zero outside its pivot row.
    """
    M = np.array(M, dtype=np.int64) % p
    if M.ndim != 2 or M.shape[0] == 0:
        return M.reshape(0, M.shape[-1] if M.ndim == 2 else 0), []
    r, piv = 0, []
    rows, cols = M.shape
    for c in range(cols):
        if r == rows:
            break
        nz = np.nonzero(M[r:, c])[0]
        if len(nz) == 0:
            continue
        q = r + nz[0]
        if q != r:
            M[[r, q]] = M[[q, r]]
        M[r] = (M[r] * pow(int(M[r, c]), -1, p)) % p
        f = M[:, c].copy()
        f[r] = 0
        M = (M - np.outer(f, M[r])) % p
        piv.append(c)
        r += 1
    return M[:r], piv


def rank(M, p):
    return len(rref(M, p)[1])


def nullspace(M, p):
    """Basis (as rows) of {x : M x = 0} over F_p."""
    M = np.array(M, dtype=np.int64) % p
    n = M.shape[1]
    R, piv = rref(M, p)
    free = [j for j in range(n) if j not in set(piv)]
    B = np.zeros((len(free), n), dtype=np.int64)
    for t, f in enumerate(free):
        B[t, f] = 1
        for row, pc in zip(R, piv):
            B[t, pc] = (-row[f]) % p
    return B


def in_rowspan(v, M, p):
    """True if vector v lies in the row space of M."""
    if len(M) == 0:
        return not (np.array(v) % p).any()
    return rank(np.vstack([M, v]), p) == rank(M, p)
