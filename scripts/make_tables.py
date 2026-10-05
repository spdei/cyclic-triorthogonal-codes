#!/usr/bin/env python3
"""Generate the appendix of best codes at each length (LaTeX) from results/.

usage: python scripts/make_tables.py [--out FILE]      (default: print to stdout)

Rules (tested in tests/test_tables.py):
  * An entry is exact only if the distance of every maximal support at that length is known
    exactly. Otherwise it is printed as >= (best lower bound), or as a--b when every support
    has a finite upper bound.
  * Concatenation (Proposition 18): if a maximal support contains pi^{-1}(S1) u k1*S2 for
    construction-A supports S1, S2 of distances d1, d2, its distance is >= d1*d2, with equality
    when it equals that support. Lower bounds are raised accordingly (to a fixpoint).
"""
import argparse, json, os, sys
from math import gcd

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
NAMES = {3: 'Qutrits', 5: 'Ququints', 7: '$p=7$'}


def load():
    """Return (binary, qudit) results. Each support gets 'lo' and 'hi' (hi = -1 if unknown)."""
    binary = {int(n): r for n, r in json.load(open(os.path.join(ROOT, 'results', 'binary_results.json'))).items()}
    qudit = [json.loads(l) for l in open(os.path.join(ROOT, 'results', 'qudit_results.jsonl'))]
    for r in qudit:
        for c in r['codes']:
            c['lo'], c['hi'] = c['d'], (c['d'] if c['status'] == 'exact' else -1)
    return binary, qudit


def a_tables(binary, qudit):
    """{p: {N: [dist record]}} for construction A; the records are shared with the results."""
    tabs = {2: {n: [dict(S=c['S'], rec=c['A']) for c in r['codes']] for n, r in binary.items() if r['codes']}}
    for r in qudit:
        if r['cons'] == 'A' and r['codes']:
            tabs.setdefault(r['p'], {})[r['N']] = [dict(S=c['S'], rec=c) for c in r['codes']]
    return tabs


def apply_concatenation(binary, qudit):
    """Raise construction-A bounds using Proposition 18. Returns a list of the changes made."""
    changes = []
    for p, tab in a_tables(binary, qudit).items():
        changed = True
        while changed:
            changed = False
            for k1 in tab:
                for k2 in tab:
                    N = k1 * k2
                    if N not in tab:
                        continue
                    units = [u for u in range(1, N) if gcd(u, N) == 1]
                    orbits = [[frozenset((u * s) % N for s in c['S']) for u in units] for c in tab[N]]
                    for c1 in tab[k1]:
                        for c2 in tab[k2]:
                            r1, r2 = c1['rec'], c2['rec']
                            T = frozenset({j for j in range(N) if j % k1 in set(c1['S'])} |
                                          {(k1 * s) % N for s in c2['S']})
                            hits = [i for i, orb in enumerate(orbits) if any(T <= M for M in orb)]
                            if not hits:
                                sys.exit(f"p={p}: concatenation {k1}x{k2} lies in no maximal support")
                            for i in hits:
                                rec = tab[N][i]['rec']
                                equal = any(T == M for M in orbits[i])
                                if equal and r1['lo'] == r1['hi'] and r2['lo'] == r2['hi']:
                                    d = r1['lo'] * r2['lo']                 # exact: d = d1*d2
                                    if (rec['hi'] >= 0 and rec['hi'] != d) or rec['lo'] > d:
                                        sys.exit(f"p={p} n={N}: recorded {rec} contradicts d={d}")
                                    if rec['hi'] != d:
                                        rec['lo'] = rec['hi'] = d
                                        changed = True
                                        changes.append((p, N, k1, k2, 'exact', d))
                                else:
                                    b = r1['lo'] * r2['lo']                 # lower bound d >= d1*d2
                                    if 0 <= rec['hi'] < b:
                                        sys.exit(f"p={p} n={N}: exact distance {rec['hi']} < bound {b}")
                                    if rec['lo'] < b:
                                        rec['lo'] = b
                                        changed = True
                                        changes.append((p, N, k1, k2, 'lower', b))
    return changes


def best(recs):
    """Best distance over a length's supports: (lo, hi), hi = -1 if it is only a lower bound."""
    lo = max(r['lo'] for r in recs)
    if all(r['hi'] >= 0 for r in recs):
        return lo, max(r['hi'] for r in recs)
    return lo, -1


def fmt(lo, hi):
    if hi == lo:
        return f"{lo}"
    if hi < 0:
        return f"\\ge{lo}"
    return f"{lo}\\text{{--}}{hi}"


def table(entries, cols):
    out = ["\\begin{center}\\small", "\\begin{tabular}{" + "l" * cols + "}"]
    for i in range(0, len(entries), cols):
        row = entries[i:i + cols]
        out.append(" & ".join(row + [""] * (cols - len(row))) + " \\\\")
    return out + ["\\end{tabular}\\end{center}", ""]


def render(binary, qudit):
    L = ["% Generated by scripts/make_tables.py from results/; do not edit by hand.",
         "\\section{Best code at each length}\\label{app:tables}",
         "For each construction and length, we list the code with the largest distance among all maximal "
         "supports, followed in parentheses by the number of maximal supports up to multiplier "
         "equivalence. Lengths at which every code has distance $1$ are omitted. An entry is exact only "
         "if the distance of every maximal support at that length is known exactly; otherwise $\\ge d$ "
         "gives the best rigorous lower bound, and $a\\text{--}b$ means that the best distance lies "
         "between $a$ and $b$. Lower bounds also use Proposition~\\ref{prop:concat}: a support "
         "containing the concatenation of two construction-A supports, of distances $d_1$ and $d_2$, "
         "has distance at least $d_1d_2$, with equality if it equals that support.",
         "",
         "\\paragraph{Qubits.} Both constructions use the same supports, all of them triply even. Each "
         "entry gives the construction-A code $[[n,1,d]]$ and, after ``B:'', the distance of the "
         "construction-B code of length $n-1$."]
    ents = []
    for n in sorted(binary):
        r = binary[n]
        if r['codes']:
            a, b = best([c['A'] for c in r['codes']]), best([c['B'] for c in r['codes']])
            ents.append(f"$[[{n},1,{fmt(*a)}]]$, B: ${fmt(*b)}$ ({r['n_orbits']})")
    L += table(ents, 3)
    for p in (3, 5, 7):
        for cons in 'AB':
            rs = sorted((r for r in qudit if r['p'] == p and r['cons'] == cons), key=lambda r: r['N'])
            ents = []
            for r in rs:
                if r['codes']:
                    lo, hi = best(r['codes'])
                    if lo >= 2:
                        ents.append(f"$[[{r['codes'][0]['n']},1,{fmt(lo, hi)}]]$ ({r['n_orbits']})")
            L.append(f"\\paragraph{{{NAMES[p]}, construction {cons}.}} "
                     f"{sum(r['n_orbits'] for r in rs)} maximal supports in total.")
            L += table(ents, 4)
    return "\n".join(L)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', default=None)
    out = ap.parse_args().out
    binary, qudit = load()
    for ch in apply_concatenation(binary, qudit):
        print("concatenation: p=%d n=%d (%dx%d) %s %d" % ch, file=sys.stderr)
    tex = render(binary, qudit)
    if out:
        open(out, 'w').write(tex + "\n")
    else:
        print(tex)


if __name__ == '__main__':
    main()
