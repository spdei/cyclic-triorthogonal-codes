"""Independent pure-Python cyclotomic data vs. the GAP/PARI data files."""
from _common import *
from triortho.gf import cyclotomic_data, is_irreducible, pmul
from triortho.cyclic import code_basis
from triortho.linalg import rank

CASES = [(2, 7), (2, 15), (2, 21), (2, 31), (2, 45), (3, 8), (3, 13), (3, 20), (3, 40),
         (5, 12), (5, 13), (5, 24), (7, 8), (7, 29), (7, 32)]

def test_minimal_polynomials_match_data_files():
    for p, n in CASES:
        mine = cyclotomic_data(n, p, seed=n)
        prod = [1]
        for c, mp in mine:
            assert len(mp) - 1 == len(c) and mp[-1] == 1
            assert is_irreducible(mp, p) or len(mp) == 2
            prod = pmul(prod, mp, p)
        assert prod == [(-1) % p] + [0] * (n - 1) + [1]          # product = x^n - 1
        ref = cosets(p, n)
        assert sorted(tuple(mp) for _, mp in mine) == sorted(tuple(mp) for _, mp in ref)
        assert sorted(len(c) for c, _ in mine) == sorted(len(c) for c, _ in ref)

def test_codes_from_independent_data_are_cyclic():
    rng = random.Random(3)
    for p, n in CASES:
        data = cyclotomic_data(n, p, seed=n)
        for _ in range(5):
            S = set()
            for c, _ in data:
                if rng.random() < 0.5:
                    S |= set(c)
            B = code_basis(n, data, S, p)
            assert len(B) == len(S) and rank(B, p) == len(S)
            Bs = np.roll(B, 1, axis=1)
            assert rank(np.vstack([B, Bs]), p) == len(S)        # shift invariant
