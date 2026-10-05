"""The appendix tables (scripts/make_tables.py).

1. Proposition 18 against the data: every concatenation of two construction-A supports lies in a
   maximal support, and no exact distance contradicts d = d1*d2 (make_tables aborts otherwise).
2. The strict rule: an entry is printed exact only if every support at that length is exact.
   Regressions: n = 465 was once printed as [[465,1,7]], although two of its supports are
   concatenations of distance exactly 9.
3. The concatenation [[15,1,3]]_2 o [[31,1,3]]_2 is a cyclic code with support pi^{-1}(S1) u 15*S2.
"""
from _common import *
sys.path.insert(0, os.path.join(ROOT, 'scripts'))
import make_tables as mt
from triortho.cyclic import code_basis
from triortho.linalg import rank


def _processed():
    binary, qudit = mt.load()
    changes = mt.apply_concatenation(binary, qudit)
    return binary, qudit, changes


def test_concatenation_bound():
    binary, qudit, changes = _processed()
    exact = {(p, N) for p, N, _, _, kind, d in changes if kind == 'exact' and d == 9}
    assert {(2, 465), (2, 675)} <= exact, changes


def test_strict_rule():
    binary, qudit, _ = _processed()
    tex = mt.render(binary, qudit)
    for n, r in binary.items():                       # exact => every support exact
        if r['codes'] and mt.best([c['A'] for c in r['codes']])[1] == mt.best([c['A'] for c in r['codes']])[0]:
            assert all(c['A']['lo'] == c['A']['hi'] for c in r['codes']), n
    for r in qudit:
        if r['codes']:
            lo, hi = mt.best(r['codes'])
            assert (lo == hi) == all(c['status'] == 'exact' or c['hi'] == c['lo'] for c in r['codes'])
    assert '[[465,1,\\ge9]]' in tex and '[[675,1,\\ge9]]' in tex and '[[341,1,\\ge9]]' in tex
    assert '[[465,1,7]]' not in tex
    for rec in ('[[15,1,3]]', '[[127,1,7]]', '[[223,1,9]]', '[[225,1,9]]'):   # records stay exact
        assert rec in tex, rec


def test_concatenation_15_31_is_cyclic():
    import json
    d = json.load(open(os.path.join(ROOT, 'results', 'binary_results.json')))
    for kin, kout in ((15, 31), (31, 15)):
        Sin, Sout = d[str(kin)]['codes'][0]['S'], d[str(kout)]['codes'][0]['S']
        Cin, Cout = code_basis(kin, cosets(2, kin), Sin, 2), code_basis(kout, cosets(2, kout), Sout, 2)
        n = kin * kout; rows = []
        for blk in range(kout):                      # block = residue class mod kout
            for c in Cin:
                v = np.zeros(n, dtype=np.int64); v[blk + kout * np.arange(kin)] = c; rows.append(v)
        for s in Cout:
            v = np.zeros(n, dtype=np.int64)
            for blk in np.nonzero(s)[0]:
                v[blk + kout * np.arange(kin)] = 1
            rows.append(v)
        E = np.array(rows) % 2; r = rank(E, 2)
        T = sorted({j for j in range(n) if j % kin in set(Sin)} | {(kin * s) % n for s in Sout})
        B = code_basis(n, cosets(2, n), T, 2)
        assert r == len(B) == rank(np.vstack([E, B]), 2)
        assert rank(np.vstack([E, np.roll(E, 1, axis=1)]), 2) == r
