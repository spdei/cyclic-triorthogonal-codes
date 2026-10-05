"""Pure-Python finite-field arithmetic, used to build cyclotomic data independently of GAP/PARI.

Polynomials over F_p are lists of ints, lowest degree first. Elements of GF(p^m) are such
polynomials reduced modulo a fixed irreducible f of degree m.
"""
import random
from math import gcd


def _trim(a):
    while len(a) > 1 and a[-1] == 0:
        a.pop()
    return a


def padd(a, b, p):
    n = max(len(a), len(b))
    return _trim([((a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0)) % p for i in range(n)])


def psub(a, b, p):
    return padd(a, [(-x) % p for x in b], p)


def pmul(a, b, p):
    r = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                r[i + j] = (r[i + j] + x * y) % p
    return _trim(r)


def pdivmod(a, b, p):
    a = list(a)
    b = _trim(list(b))
    inv = pow(b[-1], -1, p)
    q = [0] * max(1, len(a) - len(b) + 1)
    while len(a) >= len(b) and any(a):
        coef = (a[-1] * inv) % p
        sh = len(a) - len(b)
        q[sh] = coef
        for i, y in enumerate(b):
            a[sh + i] = (a[sh + i] - coef * y) % p
        _trim(a)
        if len(a) < len(b) or (len(a) == 1 and a[0] == 0):
            break
    return _trim(q), _trim(a)


def pmod(a, b, p):
    return pdivmod(a, b, p)[1]


def pgcd(a, b, p):
    a, b = _trim(list(a)), _trim(list(b))
    while any(b):
        a, b = b, pmod(a, b, p)
    inv = pow(a[-1], -1, p)
    return [(x * inv) % p for x in a]


def ppowmod(a, e, f, p):
    r, base = [1], pmod(a, f, p)
    while e:
        if e & 1:
            r = pmod(pmul(r, base, p), f, p)
        base = pmod(pmul(base, base, p), f, p)
        e >>= 1
    return r


def prime_factors(n):
    fs, d = set(), 2
    while d * d <= n:
        while n % d == 0:
            fs.add(d)
            n //= d
        d += 1
    if n > 1:
        fs.add(n)
    return sorted(fs)


def is_irreducible(f, p):
    """Rabin's test for monic f of degree m over F_p."""
    m = len(f) - 1
    x = [0, 1]
    xp = ppowmod(x, p ** m, f, p)
    if _trim(psub(xp, x, p)) != [0]:
        return False
    for q in prime_factors(m):
        h = psub(ppowmod(x, p ** (m // q), f, p), x, p)
        if len(pgcd(f, h, p)) != 1:
            return False
    return True


def random_irreducible(m, p, rng):
    while True:
        f = [rng.randrange(p) for _ in range(m)] + [1]
        if f[0] != 0 and is_irreducible(f, p):
            return f


def mult_order(a, n):
    k, x = 1, a % n
    while x != 1:
        x = (x * a) % n
        k += 1
    return k


def cyclotomic_cosets(n, p):
    seen, out = set(), []
    for r in range(n):
        if r in seen:
            continue
        c, s = [], r
        while s not in seen:
            seen.add(s)
            c.append(s)
            s = (s * p) % n
        out.append(tuple(sorted(c)))
    return out


def cyclotomic_data(n, p, seed=0):
    """Cosets of Z_n under multiplication by p, each with the minimal polynomial over F_p of
    alpha^s (s in the coset), for one fixed primitive n-th root of unity alpha.
    Returns a list of (coset, coefficients-low-to-high)."""
    assert gcd(n, p) == 1
    rng = random.Random(seed)
    m = mult_order(p, n)
    f = random_irreducible(m, p, rng) if m > 1 else [0, 1]
    N = p ** m - 1
    # field multiplication in GF(p^m) = F_p[x]/(f)
    mul = (lambda a, b: pmod(pmul(a, b, p), f, p)) if m > 1 else (lambda a, b: [(a[0] * b[0]) % p])
    def power(a, e):
        r, base = [1], a
        while e:
            if e & 1:
                r = mul(r, base)
            base = mul(base, base)
            e >>= 1
        return r
    # element of order exactly n
    while True:
        beta = _trim([rng.randrange(p) for _ in range(m)])
        if beta == [0]:
            continue
        alpha = power(beta, N // n)
        if all(power(alpha, n // q) != [1] for q in prime_factors(n)) and power(alpha, n) == [1]:
            break
    out = []
    for c in cyclotomic_cosets(n, p):
        # minimal polynomial: prod_{t in c} (X - alpha^t), coefficients in GF(p^m)
        poly = [[1]]                       # list of field elements, lowest degree first
        for t in c:
            root = power(alpha, t)
            neg = [(-x) % p for x in root]
            new = [[0] for _ in range(len(poly) + 1)]
            for i, a in enumerate(poly):   # multiply by (X + neg)
                new[i + 1] = padd(new[i + 1], a, p)
                new[i] = padd(new[i], mul(a, neg), p)
            poly = new
        coeffs = []
        for a in poly:
            a = _trim(list(a))
            assert len(a) == 1, "minimal polynomial must have coefficients in F_p"
            coeffs.append(a[0])
        out.append((c, coeffs))
    return out
