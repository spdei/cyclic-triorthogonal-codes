"""Direct checks of the definitions.

* tri_check        : conditions (T1),(T2) of the triorthogonality definition (any prime p,
                     optional weights). For p = 2 this is the Bravyi-Haah condition with the
                     row partition by weight parity.
* cube_orthogonal  : tau vanishes identically on a space (a triorthogonal space).
* transversal_check: exhaustive verification of the transversal gate on all code-state
                     components, with the exact phase arithmetic (mod 8, 9 or p).
"""
import itertools
import numpy as np


def tau_tensor_vanishes(B, p, weights=None, skip_logical_diag_k=0):
    """True iff sum_l w_l B[a,l]B[b,l]B[c,l] = 0 (mod p) for all a<=b<=c, except
    (i,i,i) with i < skip_logical_diag_k."""
    B = np.array(B, dtype=np.int64) % p
    w = np.ones(B.shape[1], dtype=np.int64) if weights is None else np.array(weights, dtype=np.int64) % p
    m = len(B)
    for i in range(m):
        P = ((B * B[i] * w) % p) @ B.T % p          # P[j,l] = tau(i,j,l)
        for j in range(i, m):
            for l in range(j, m):
                if P[j, l] and not (i == j == l and i < skip_logical_diag_k):
                    return False
    return True


def tri_check(H, p, k, weights=None):
    """(T1) and (T2) for H whose first k rows are logical."""
    H = np.array(H, dtype=np.int64) % p
    if not tau_tensor_vanishes(H, p, weights, skip_logical_diag_k=k):
        return False
    w = np.ones(H.shape[1], dtype=np.int64) if weights is None else np.array(weights, dtype=np.int64) % p
    return all(int(((H[i] ** 3) * w).sum() % p) != 0 for i in range(k))


def cube_orthogonal(W, p, weights=None):
    return tau_tensor_vanishes(W, p, weights)


def _phase_table(p, gate):
    """Exponent table phi[j] and modulus q for the single-qudit gate exp(2 pi i phi(j)/q)."""
    if p == 2:
        return [0, 1], 8                                  # T = diag(1, e^{i pi/4})
    if p == 3:
        if gate == 'U':
            return [0, 1, 8], 9                           # zeta^{l(j)^3} = diag(1,z,z^-1)
        if gate == 'V':
            return [0, 1, 2], 9                           # zeta^{l(j)}  = diag(1,z,z^2)
        raise ValueError(gate)
    return [(j ** 3) % p for j in range(p)], p            # omega^{j^3}


def transversal_check(H, k, p, gate='U', weights=None, max_states=3 ** 11):
    """Apply the transversal gate (gate^{w_l} on qudit l) to every component x*H1 + a*H0 of
    every code state and analyse the phases.

    Returns dict(status, logical) where status is
        'exact'  : phases independent of the stabilizer coefficients a,
        'pauli'  : stabilizer dependence equals that of a fixed Pauli Z(w),
        'fail'   : the code space is not preserved (up to a Pauli),
    and logical[x] is the exponent (mod q) acquired by |x-bar> after the correction."""
    H = np.array(H, dtype=np.int64) % p
    m, n = H.shape
    phi, q = _phase_table(p, gate)
    phi = np.array(phi, dtype=np.int64)
    w = np.ones(n, dtype=np.int64) if weights is None else np.array(weights, dtype=np.int64)
    if p ** m > max_states:
        raise ValueError("code too large for exhaustive check")
    H1, H0 = H[:k], H[k:]
    r = len(H0)
    A = np.array(list(itertools.product(range(p), repeat=r)), dtype=np.int64).reshape(p ** r, r)
    stab = (A @ H0) % p if r else np.zeros((1, n), dtype=np.int64)
    unit = q // p
    status, logical, L_ref = 'exact', {}, None
    for x in itertools.product(range(p), repeat=k):
        y = (np.array(x, dtype=np.int64) @ H1) % p if k else np.zeros(n, dtype=np.int64)
        C = (stab + y) % p
        E = (phi[C] * w).sum(axis=1) % q
        diff = (E - E[0]) % q
        if diff.any():
            if (diff % unit).any():
                return dict(status='fail', logical=None)
            # must be linear in a with a coefficient vector independent of x
            L = np.array([diff[p ** (r - 1 - j)] for j in range(r)], dtype=np.int64)  # e_j rows
            if ((A @ L) % q != diff).any():
                return dict(status='fail', logical=None)
            if L_ref is None:
                L_ref = L
            elif (L_ref != L).any():
                return dict(status='fail', logical=None)
            status = 'pauli'
        elif L_ref is not None and L_ref.any():
            return dict(status='fail', logical=None)
        else:
            L_ref = np.zeros(r, dtype=np.int64) if L_ref is None else L_ref
        logical[x] = int(E[0])
    return dict(status=status, logical=logical)


def generalised_check(H, p, k, weights=None):
    """Generalised triorthogonality: tau_w(a,b,c) = 0 whenever at least one of a,b,c is a
    stabilizer row (rows k..m-1); the logical rows are unconstrained."""
    H = np.array(H, dtype=np.int64) % p
    w = np.ones(H.shape[1], dtype=np.int64) if weights is None else np.array(weights, dtype=np.int64) % p
    m = len(H)
    for i in range(m):
        P = ((H * H[i] * w) % p) @ H.T % p
        for j in range(m):
            for l in range(m):
                if P[j, l] and max(i, j, l) >= k:
                    return False
    return True
