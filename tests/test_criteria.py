"""Propositions A and B (sumset criteria) vs. direct trilinear checks, and the search."""
from _common import *
from triortho.cyclic import code_basis, sumset
from triortho.checks import tri_check, cube_orthogonal
from triortho.search import admissible, maximal_supports

def test_sumset_criteria_both_directions():
    rng = random.Random(11)
    seen = {True: 0, False: 0}
    for p, ns in ((2, [7, 15, 21, 31, 45, 63]), (3, [8, 10, 13, 20, 26, 40]), (5, [8, 12, 13, 24]), (7, [8, 12, 16, 24, 29])):
        for n in ns:
            for _ in range(15):
                S = random_support(p, n, rng)
                if not S:
                    continue
                B = code_basis(n, cosets(p, n), S, p)
                # construction B: W triorthogonal space  <=>  0 not in 3S
                okB = cube_orthogonal(B, p)
                assert okB == (0 not in sumset(S, n, 3)) == admissible(S, n, 'B')
                seen[okB] += 1
                # construction A: [1; C] triorthogonal  <=>  0 not in S, 2S, 3S
                if p == 2 or n % p:
                    H = np.vstack([np.ones(n, dtype=np.int64), B])
                    expect = 0 not in S and 0 not in sumset(S, n, 2) and 0 not in sumset(S, n, 3)
                    if 0 not in S:
                        assert tri_check(H, p, 1) == expect == admissible(S, n, 'A')
    assert seen[True] > 5 and seen[False] > 5        # both outcomes exercised

def test_search_is_complete_and_maximal():
    """Compare maximal_supports with brute force over all unions of cosets (small n)."""
    from itertools import combinations
    from triortho.cyclic import canonical_support
    for p, n in ((2, 15), (2, 21), (3, 13), (3, 20), (5, 12), (7, 12)):
        cs = cosets(p, n)
        for cons in 'AB':
            pool = [c for c, _ in cs if not (cons == 'A' and c == (0,))]
            adm = []
            for r in range(1, len(pool) + 1):
                for sub in combinations(pool, r):
                    S = set().union(*map(set, sub))
                    if admissible(S, n, cons):
                        adm.append(S)
            maxi = [S for S in adm if not any(S < T for T in adm)]
            brute = sorted({canonical_support(S, n) for S in maxi})
            mine, complete = maximal_supports(n, cs, cons)
            assert complete and [k for k, _ in mine] == brute, (p, n, cons)
