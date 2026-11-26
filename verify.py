#!/usr/bin/env python3
"""Re-check every counterexample in this repository.  Python 3 standard library only.

    python verify.py            all checks, a few minutes
    python verify.py --quick    skip the slow brute-force cross-checks (about a minute)

Each check enumerates ALL admissible choices of the "sparse" terms (powers of 2, 3, ...) and
shows that what is left over is never of the required shape.  For the forms below this is a
finite computation, so a passing check is a proof that the number has no representation.

Membership tests used
  x^2 + 5y^2, polygonal sums : plain search (loop over one variable, integer square root).
  x^2 + 2y^2                 : plain search (default) and, independently, the classical
                               criterion "every prime == 5 or 7 (mod 8) divides m to an even
                               power" (Z[sqrt(-2)] is a principal ideal domain).
  primes                     : explicit divisors (and, for A303932, square roots of 11) from the
                               certificate files, so no primality test is trusted; a
                               deterministic Miller-Rabin test is used for recounts only.
  x^2 + y^2 (A300510,        : the certificate files give, for every number, a divisor d == 3
  A300448) and x^2 + 2y^2      (mod 4), respectively d == 5 or 7 (mod 8), coprime to its
  (A301452)                    cofactor, which needs no primality test; the recounts factor
                               every number completely (Miller-Rabin and Pollard's rho) and use
                               the criterion "every prime == 3 (mod 4), respectively == 5 or 7
                               (mod 8), divides v to an even power".
"""
import sys
from math import gcd, isqrt

QUICK = "--quick" in sys.argv


def is_square(v):
    if v < 0:
        return False
    r = isqrt(v)
    return r * r == v


# ---------------------------------------------------------------- x^2 + b*y^2 by plain search
def has_x2_by2(m, b):
    y = 0
    while b * y * y <= m:
        if is_square(m - b * y * y):
            return True
        y += 1
    return False


# ---------------------------------------------------------------- x^2 + 2y^2 by the prime criterion
_SMALL = None


def small_primes(limit=22000):
    global _SMALL
    if _SMALL is None:
        s = bytearray([1]) * (limit + 1)
        s[0:2] = b"\x00\x00"
        for i in range(2, isqrt(limit) + 1):
            if s[i]:
                s[i * i::i] = bytearray(len(s[i * i::i]))
        _SMALL = [i for i in range(limit + 1) if s[i]]
    return _SMALL


def is_prime(n):
    """Deterministic Miller-Rabin with the first 13 primes as bases, valid for n < 3.3e24
    (Sorenson and Webster, 2017)."""
    if n < 2:
        return False
    bases = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41)
    for q in bases:
        if n % q == 0:
            return n == q
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for a in bases:
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def has_x2_2y2_criterion(m):
    """m = x^2 + 2y^2  iff  every prime q == 5, 7 (mod 8) divides m to an even power (m < 1e13)."""
    if m == 0:
        return True
    while m % 2 == 0:
        m //= 2
    for p in small_primes():
        if p * p > m:
            break
        if m % p == 0:
            e = 0
            while m % p == 0:
                m //= p
                e += 1
            if p % 8 in (5, 7) and e % 2:
                return False
    if m == 1:
        return True
    if m < 22000 ** 2 or is_prime(m):          # a single prime is left
        return m % 8 in (1, 3)
    r = isqrt(m)                               # m = p*p or p*q with p, q > 22000 (since m < 22000^3)
    if r * r == m:
        return True
    return m % 8 in (1, 3) and _two_good_primes(m)


def _two_good_primes(m):
    """m = p*q, p < q primes above 22000: are both == 1, 3 (mod 8)?  Fermat factorisation is enough here."""
    a = isqrt(m) + 1
    while not is_square(a * a - m):
        a += 1
    b = isqrt(a * a - m)
    return (a - b) % 8 in (1, 3) and (a + b) % 8 in (1, 3)


# ---------------------------------------------------------------- the individual claims
def count_a302984(n, member):
    """Number of pairs (z, w) with n - 2^z - 5*2^w of the form x^2 + 2y^2."""
    good = pairs = 0
    z = 0
    while 2 ** z <= n:
        w = 0
        while 2 ** z + 5 * 2 ** w <= n:
            good += member(n - 2 ** z - 5 * 2 ** w)
            pairs += 1
            w += 1
        z += 1
    return good, pairs


def count_a302982(n):
    good = pairs = 0
    z = 0
    while 2 ** z <= n:
        w = 0
        while 2 ** z + 3 * 2 ** w <= n:
            good += has_x2_by2(n - 2 ** z - 3 * 2 ** w, 5)
            pairs += 1
            w += 1
        z += 1
    return good, pairs


def gen_poly_tester(A, B):
    """is v = d*(A*d + B)/2 for some integer d (positive, zero or negative)?"""
    def t(v):
        D = B * B + 8 * A * v
        if not is_square(D):
            return False
        r = isqrt(D)
        return (-B + r) % (2 * A) == 0 or (-B - r) % (2 * A) == 0
    return t


def count_smooth_square(n, p, q, f, g_test):
    """Number of (a, b, c) with n - (p^a q^b)^2 - f(c) passing g_test (c over all integers)."""
    good = 0
    A = 1
    while A * A <= n:
        Bv = A
        while Bv * Bv <= n:
            m = n - Bv * Bv
            for sign in (1, -1):
                c = 0 if sign == 1 else -1
                while f(c) <= m:
                    good += g_test(m - f(c))
                    c += sign
            Bv *= q
        A *= p
    return good


def check_cert_a303932(N, path, mode):
    """mode 'nr': no prime p with 11 a non-residue; mode 'qr': no prime p with 11 a residue."""
    cert = {}
    for line in open(path):
        if line.startswith("#") or line.startswith("k,"):
            continue
        k, m, p, kind, w = line.strip().split(",")
        cert[(int(k), int(m))] = (int(p), kind, int(w))
    need = []
    k = 0
    while 2 ** k < N:
        m = 0
        while 2 ** k + 3 ** m < N:
            need.append((k, m))
            m += 1
        k += 1
    assert set(need) == set(cert) and len(need) == len(cert)
    for (k, m) in need:
        p, kind, w = cert[(k, m)]
        assert p == N - 2 ** k - 3 ** m and p >= 1
        if kind == "small":
            assert p == 1
        elif kind == "divisor":
            assert 1 < w < p and p % w == 0
        elif kind == "sqrt":                   # 11 is a square mod p
            assert mode == "nr" and w * w % p == 11 % p and p != 11
        elif kind == "nonres":                 # Euler: if p is prime, 11 is not a square mod p
            assert mode == "qr" and p > 2 and p % 2 == 1 and p != 11 and pow(11, (p - 1) // 2, p) == p - 1
        else:
            raise AssertionError(kind)
    return len(need), sum(1 for v in cert.values() if v[1] in ("sqrt", "nonres"))


def recount_a303932(N):
    """Independent recount with a primality test: (primes p, with 11 a residue, with 11 a non-residue).
    The primes 2 and 11 are counted as residues, which can only make Conjecture 1 easier."""
    primes = res = nonres = 0
    k = 0
    while 2 ** k < N:
        m = 0
        while 2 ** k + 3 ** m < N:
            p = N - 2 ** k - 3 ** m
            if p >= 2 and is_prime(p):
                primes += 1
                if p in (2, 11) or pow(11, (p - 1) // 2, p) == 1:
                    res += 1
                else:
                    nonres += 1
            m += 1
        k += 1
    return primes, res, nonres


# ---------------------------------------------------------------- A157242: p + 2^x + 11*2^y, p == 5 (mod 6)
def pairs_p2k(N, k):
    """All (x, y) with x, y >= 1 and 2^x + k*2^y < N."""
    out = []
    y = 1
    while k * 2 ** y < N:
        x = 1
        while 2 ** x + k * 2 ** y < N:
            out.append((x, y))
            x += 1
        y += 1
    return out


def count_a157242(N, k=11, c=5):
    """Number of (x, y), x, y >= 1, such that N - 2^x - k*2^y is a prime congruent to c mod 6."""
    cnt = 0
    for x, y in pairs_p2k(N, k):
        p = N - 2 ** x - k * 2 ** y
        if p % 6 == c and is_prime(p):
            cnt += 1
    return cnt


def check_cert_a157242(N, path, k=11, c=5):
    """Every pair (x, y) is listed; p is in the wrong class mod 6 or has an explicit proper divisor."""
    cert = {}
    for line in open(path):
        if line.startswith("#") or line.startswith("x,"):
            continue
        x, y, p, pm, d = map(int, line.split(","))
        cert[(x, y)] = (p, pm, d)
    need = pairs_p2k(N, k)
    assert set(need) == set(cert) and len(need) == len(cert)
    inclass = 0
    for (x, y) in need:
        p, pm, d = cert[(x, y)]
        assert p == N - 2 ** x - k * 2 ** y and p >= 1 and pm == p % 6
        if pm == c:
            inclass += 1
            assert p == 1 or (1 < d < p and p % d == 0)
    return len(need), inclass


# ---------------------------------------------------------------- squares minus m*4^k (A300510, A300448)
_FACTORS = {}


def _rho(n):
    """A proper factor of the odd composite n (Pollard's rho with fixed starting values)."""
    c = 1
    while True:
        x = y = 2
        g = 1
        while g == 1:
            x = (x * x + c) % n
            y = (y * y + c) % n
            y = (y * y + c) % n
            g = gcd(abs(x - y), n)
        if g != n:
            return g
        c += 1


def factorize(n):
    """dict prime -> exponent for 0 < n < 3.3e24 (trial division, Miller-Rabin, Pollard's rho)."""
    if n in _FACTORS:
        return _FACTORS[n]
    out, m = {}, n
    for q in small_primes():
        if q * q > m:
            break
        while m % q == 0:
            out[q] = out.get(q, 0) + 1
            m //= q
    stack = [m] if m > 1 else []
    while stack:
        t = stack.pop()
        if is_prime(t):
            out[t] = out.get(t, 0) + 1
        else:
            d = _rho(t)
            stack += [d, t // d]
    _FACTORS[n] = out
    return out


def bad_prime(q, form):
    """The primes that must divide a number x^2 + y^2 (form 1) or x^2 + 2y^2 (form 2) to an even power."""
    return q % 4 == 3 if form == 1 else q % 8 in (5, 7)


def represented(v, form=1):
    """Is v >= 0 of the form x^2 + y^2 (form 1) or x^2 + 2y^2 (form 2)?  Prime criterion."""
    if v == 0:
        return True
    return all(e % 2 == 0 for q, e in factorize(v).items() if q != 2 and bad_prime(q, form))


def sq_values(n, mults):
    """All v = n^2 - m*4^k with m in mults, k >= 0 and m*4^k <= n^2."""
    out = []
    for m in mults:
        t = m
        while t <= n * n:
            out.append(n * n - t)
            t *= 4
    return out


def count_sq(n, mults, form=1):
    """(number of represented values, number of values)"""
    vs = sq_values(n, mults)
    return sum(represented(v, form) for v in vs), len(vs)


def ways_x_le_y(v):
    """Number of (x, y) with 0 <= x <= y and x^2 + y^2 = v, by plain search (small v)."""
    c, x = 0, 0
    while 2 * x * x <= v:
        if is_square(v - x * x):
            c += 1
        x += 1
    return c


def ways_from_factors(v):
    """The same number for v > 0 from the factorization: (r + s) / 2, where r is the product of e + 1
    over the primes == 1 (mod 4), or 0 if some prime == 3 (mod 4) has an odd exponent, and s = 1 if v
    is a square or twice a square, else 0."""
    f = factorize(v)
    if any(e % 2 for q, e in f.items() if q % 4 == 3):
        return 0
    r = 1
    for q, e in f.items():
        if q % 4 == 1:
            r *= e + 1
    return (r + (1 if is_square(v) or (v % 2 == 0 and is_square(v // 2)) else 0)) // 2


def ways_x2_2y2(v):
    """Number of (x, y) with x, y >= 0 and x^2 + 2y^2 = v, by plain search (small v)."""
    c, y = 0, 0
    while 2 * y * y <= v:
        if is_square(v - 2 * y * y):
            c += 1
        y += 1
    return c


def check_cert_squares(n, mults, path, form=1):
    """One row m,k,v,d for every pair (m, k) with m*4^k < n^2: v = n^2 - m*4^k, d odd, d | v,
    gcd(d, v/d) = 1 and d == 3 (mod 4) (form 1) or d == 5, 7 (mod 8) (form 2).  No primality test."""
    need = {}
    for m in mults:
        k = 0
        while m * 4 ** k < n * n:
            need[(m, k)] = n * n - m * 4 ** k
            k += 1
    seen = set()
    for line in open(path):
        parts = line.strip().split(",")
        if len(parts) != 4 or not parts[0].isdigit():
            continue
        m, k, v, d = (int(t) for t in parts)
        assert (m, k) in need and (m, k) not in seen and v == need[(m, k)]
        assert d > 1 and d % 2 == 1 and v % d == 0 and gcd(d, v // d) == 1
        assert bad_prime(d, form)
        seen.add((m, k))
    assert seen == set(need)
    return len(need)


OK = True


def report(label, n, value, expect_zero=True, extra=""):
    global OK
    good = (value == 0) if expect_zero else (value > 0)
    OK &= good
    print("%-4s %-34s n = %-13d representations: %-5d %s" % ("PASS" if good else "FAIL", label, n, value, extra), flush=True)


def main():
    global OK
    print("A157242  odd N = p + 2^x + 11*2^y, p prime, p == 5 (mod 6), x, y >= 1"
          "  (conjectured for all odd N > 28 except 33, 39, 45, 171, 247)")
    oeis = [0] * 14 + [1, 1, 0, 2, 1, 0, 2, 2, 0, 2, 2, 1, 3, 2, 2, 3, 1, 2, 4, 2, 2, 5, 1, 2, 5, 2, 2, 4, 2, 2,
                       3, 2, 3, 4, 4, 3, 6, 2, 3, 6, 5, 1, 7, 4, 2, 6]
    same = [count_a157242(2 * n - 1) for n in range(1, 61)] == oeis
    OK &= same
    print("%-4s the count function reproduces the 60 terms a(1..60) of A157242 given in the OEIS" % ("PASS" if same else "FAIL"))
    top = 2001 if QUICK else 20001
    zeros = [N for N in range(29, top + 1, 2) if count_a157242(N) == 0]
    good = zeros == [33, 39, 45, 171, 247]
    OK &= good
    print("%-4s odd N in (28, %d] with no representation: %s" % ("PASS" if good else "FAIL", top, zeros))
    N = 107576474238337
    pairs, inclass = check_cert_a157242(N, "certificates/A157242_N%d.csv" % N)
    print("PASS certificate (N = %d): all %d pairs (x, y) covered; %d have p == 5 (mod 6), each with an explicit proper divisor"
          % (N, pairs, inclass))
    report("A157242 (recount)", N, count_a157242(N), extra="of %d pairs (x, y)" % pairs)
    for M in (N - 2, N + 2):
        report("  control (neighbour)", M, count_a157242(M), expect_zero=False)
    print("Same form, other cases: k = 11 with p == 1 (mod 6) (A157237, open), and k = 7 (exceptions found by D. S. McNeil in 2009, and one more)")
    zeros = [M for M in range(31, top + 1, 2) if count_a157242(M, 11, 1) == 0]
    good = zeros == [35, 41, 47, 101, 167, 2021][:5 if QUICK else 6]
    OK &= good
    print("%-4s k=11, p == 1 (mod 6): odd N in (30, %d] with no representation: %s" % ("PASS" if good else "FAIL", top, zeros))
    report("k=11, p == 1 (mod 6)", 119171, count_a157242(119171, 11, 1))
    report("k=7,  p == 1 (mod 6)", 445466268731, count_a157242(445466268731, 7, 1))
    for M in (1433987797, 45483661477, 85305551767, 367084269651, 1572650032347):
        report("k=7,  p == 5 (mod 6)", M, count_a157242(M, 7, 5))

    print("A300510 / A300448  for every n > 1 some n^2 - m*4^k (k >= 0) is a sum of two squares:"
          "  m in {2, 5} (Sun's 2-5 conjecture, A300510); m in {1, 2, 5} (A300448, Conjecture 1)")
    d510 = [0, 1, 2, 1, 3, 3, 3, 1, 4, 4, 5, 3, 5, 4, 6, 1, 3, 4, 5, 4, 7, 6, 5, 3, 8, 6, 6, 4, 5, 7, 7, 1, 5, 4, 11, 4,
            7, 5, 6, 4, 6, 8, 5, 6, 12, 5, 5, 3, 6, 9, 8, 6, 7, 6, 10, 4, 7, 7, 6, 7, 5, 9, 9, 1, 8, 5, 10, 4, 9, 11,
            9, 4, 11, 7, 12, 5, 8, 7, 7, 4]
    d448 = [1, 2, 4, 2, 4, 5, 4, 2, 7, 5, 6, 5, 6, 5, 5, 2, 5, 7, 6, 5, 10, 7, 5, 5, 7, 7, 6, 5, 6, 6, 7, 2, 10, 6, 8, 7,
            7, 6, 6, 5, 7, 11, 5, 7, 10, 5, 6, 5, 8, 8, 9, 7, 8, 6, 6, 5, 10, 7, 5, 6, 5, 8, 7, 2, 5, 10, 8, 6, 10, 8,
            7, 7, 11, 7, 7, 6, 10, 7, 5, 5]
    same = ([sum(ways_x_le_y(v) for v in sq_values(n, (2, 5))) for n in range(1, 81)] == d510
            and [sum(1 for v in sq_values(n, (1, 2, 5)) if ways_x_le_y(v)) for n in range(1, 81)] == d448)
    OK &= same
    print("%-4s plain search reproduces the 80 terms a(1..80) of A300510 and of A300448 given in the OEIS" % ("PASS" if same else "FAIL"))
    same = ([sum(ways_from_factors(v) for v in sq_values(n, (2, 5))) for n in range(1, 81)] == d510
            and [count_sq(n, (1, 2, 5))[0] for n in range(1, 81)] == d448)
    OK &= same
    print("%-4s the prime criterion (by factoring) reproduces the same 160 terms" % ("PASS" if same else "FAIL"))
    top = 2001 if QUICK else 20001
    zeros = [n for n in range(2, top + 1) if not any(represented(v) for v in sq_values(n, (2, 5)))]
    good = zeros == []
    OK &= good
    print("%-4s 1 < n <= %d with no representation for m in {2, 5}: %s" % ("PASS" if good else "FAIL", top, zeros))
    n0 = 80397893435
    rows = check_cert_squares(n0, (2, 5), "certificates/A300510_n%d.csv" % n0)
    print("PASS certificate (n = %d, m in {2, 5}): each of the %d numbers n^2 - m*4^k has a divisor d == 3 (mod 4) coprime to its cofactor" % (n0, rows))
    rows = check_cert_squares(n0, (1, 2, 5), "certificates/A300448_conj1_n%d.csv" % n0)
    print("PASS certificate (n = %d, m in {1, 2, 5}): each of the %d numbers n^2 - m*4^k has a divisor d == 3 (mod 4) coprime to its cofactor" % (n0, rows))
    g, tot = count_sq(n0, (2, 5))
    report("A300510 (recount by factoring)", n0, g, extra="of %d numbers n^2 - m*4^k" % tot)
    g, tot = count_sq(n0, (1, 2, 5))
    report("A300448 Conj.1 (recount)", n0, g, extra="of %d numbers n^2 - m*4^k" % tot)
    n2 = 239694783
    rows = check_cert_squares(n2, (1, 5, 17), "certificates/A300448_conj2_n%d.csv" % n2)
    print("PASS certificate (n = %d, m in {1, 5, 17}): each of the %d numbers n^2 - m*4^k has a divisor d == 3 (mod 4) coprime to its cofactor" % (n2, rows))
    g, tot = count_sq(n2, (1, 5, 17))
    report("A300448 Conj.2 (recount)", n2, g, extra="of %d numbers n^2 - m*4^k" % tot)
    for M, e in ((n2 - 2, 17), (n2 + 2, 25)):
        g, tot = count_sq(M, (1, 5, 17))
        good = (g, tot) == (e, 81)
        OK &= good
        print("%-4s   control (neighbour)                n = %-13d sums of two squares: %d of %d" % ("PASS" if good else "FAIL", M, g, tot), flush=True)
    zeros = [n for n in range(1, top + 1) if not any(represented(v) for v in sq_values(n, (1, 5, 17)))]
    good = zeros == []
    OK &= good
    print("%-4s 0 < n <= %d with no representation for m in {1, 5, 17}: %s" % ("PASS" if good else "FAIL", top, zeros))
    for M, e25, e125, a510 in ((n0 - 2, 19, 27, 47), (n0 + 2, 18, 22, 39)):
        got = (count_sq(M, (2, 5))[0], count_sq(M, (1, 2, 5))[0], sum(ways_from_factors(v) for v in sq_values(M, (2, 5))))
        good = got == (e25, e125, a510)
        OK &= good
        print("%-4s   control (neighbour)                n = %-13d sums of two squares: %d of 72 and %d of 109; A300510(n) = %d"
              % (("PASS" if good else "FAIL", M) + got), flush=True)
    for n3, tot3, ctl in ((68136270467, 104, ((-2, 21), (2, 26))), (80570407837, 107, ((-2, 19), (2, 20)))):
        rows = check_cert_squares(n3, (1, 5, 65), "certificates/A300448_conj3_n%d.csv" % n3)
        print("PASS certificate (n = %d, m in {1, 5, 65}): each of the %d numbers n^2 - m*4^k has a divisor d == 3 (mod 4) coprime to its cofactor" % (n3, rows))
        g, tot = count_sq(n3, (1, 5, 65))
        OK &= tot == tot3
        report("A300448 Conj.3 (recount)", n3, g, extra="of %d numbers n^2 - m*4^k" % tot)
        for dlt, e in ctl:
            g, tot = count_sq(n3 + dlt, (1, 5, 65))
            good = (g, tot) == (e, tot3)
            OK &= good
            print("%-4s   control (neighbour)                n = %-13d sums of two squares: %d of %d" % ("PASS" if good else "FAIL", n3 + dlt, g, tot), flush=True)
    zeros = [n for n in range(1, top + 1) if not any(represented(v) for v in sq_values(n, (1, 5, 65)))]
    good = zeros == []
    OK &= good
    print("%-4s 0 < n <= %d with no representation for m in {1, 5, 65}: %s" % ("PASS" if good else "FAIL", top, zeros))

    print("A301452  for every n > 1 some n^2 - 2*4^k or n^2 - 3*4^k (k >= 0) is of the form x^2 + 2y^2  (Sun's 2-3 conjecture)")
    d452 = [0, 2, 2, 2, 2, 5, 3, 2, 4, 4, 4, 5, 5, 5, 6, 2, 4, 6, 5, 4, 9, 5, 4, 5, 5, 7, 10, 5, 6, 7, 8, 2, 6, 6, 7, 6,
            9, 7, 10, 4, 6, 12, 3, 5, 10, 5, 6, 5, 5, 8, 9, 7, 7, 12, 5, 5, 13, 9, 6, 7, 8, 10, 13, 2, 6, 8, 10, 6, 15,
            9, 9, 6, 10, 9, 12, 7, 8, 13, 6, 4]
    same = [sum(ways_x2_2y2(v) for v in sq_values(n, (2, 3))) for n in range(1, 81)] == d452
    OK &= same
    print("%-4s plain search reproduces the 80 terms a(1..80) of A301452 given in the OEIS" % ("PASS" if same else "FAIL"))
    same = all((ways_x2_2y2(v) > 0) == represented(v, 2) for n in range(1, 201) for v in sq_values(n, (2, 3)))
    OK &= same
    print("%-4s the prime criterion (by factoring) agrees with plain search for all n^2 - m*4^k with n <= 200" % ("PASS" if same else "FAIL"))
    zeros = [n for n in range(2, top + 1) if not any(represented(v, 2) for v in sq_values(n, (2, 3)))]
    good = zeros == []
    OK &= good
    print("%-4s 1 < n <= %d with no representation for m in {2, 3}: %s" % ("PASS" if good else "FAIL", top, zeros))
    for n4 in (386749710995, 985224034231, 1241790594939, 1572550146081, 9224999068015):
        rows = check_cert_squares(n4, (2, 3), "certificates/A301452_n%d.csv" % n4, 2)
        print("PASS certificate (n = %d, m in {2, 3}): each of the %d numbers n^2 - m*4^k has a divisor d == 5 or 7 (mod 8) coprime to its cofactor" % (n4, rows))
        if n4 * n4 < 33 * 10 ** 23:             # the range in which is_prime, used by the factoring, is proven
            g, tot = count_sq(n4, (2, 3), 2)
            report("A301452 (recount by factoring)", n4, g, extra="of %d numbers n^2 - m*4^k" % tot)
    for M, e in ((386749710993, 17), (386749710997, 16)):
        g, tot = count_sq(M, (2, 3), 2)
        good = (g, tot) == (e, 76)
        OK &= good
        print("%-4s   control (neighbour)                n = %-13d of the form x^2 + 2y^2: %d of %d" % ("PASS" if good else "FAIL", M, g, tot), flush=True)

    print("Published failures of the same kind, re-checked: m in {2, c} for c = 29, 65, 113 (Q.-H. Hou) and c = 162 (Sun);"
          " m in {1, 5, 11} (Sun); Sun's counterexample to his Square Conjecture A301471 (n^2 - 3*2^k = x^2 + 2y^2)")
    for c, M in ((29, 51883659), (65, 56173837), (113, 65525725), (162, 42211965)):
        g, tot = count_sq(M, (2, c))
        report("m in {2, %d}" % c, M, g, extra="of %d numbers" % tot)
    g, tot = count_sq(13018247, (1, 5, 11))
    report("m in {1, 5, 11}", 13018247, g, extra="of %d numbers" % tot)
    g, tot = count_sq(5884015571, (3, 6), 2)
    report("A301471 (x^2 + 2y^2, 3*2^k)", 5884015571, g, extra="of %d numbers" % tot)

    A302984 = [4375339279, 6285378879, 8525602079, 9731038079, 13523100479, 24920030079, 29312698479,
               30034054079, 32269674879, 40072842719, 51729946479, 58532818879, 67179857279, 71235478079,
               72161342879, 82358435279, 85090270079, 91924062079]
    print("A302984  n = x^2 + 2y^2 + 2^z + 5*2^w  (conjectured for all n > 5)")
    for n in A302984:
        g, pairs = count_a302984(n, has_x2_2y2_criterion)
        report("A302984 (prime criterion)", n, g, extra="of %d pairs (z, w)" % pairs)
    for n in (A302984[0] - 1, A302984[0] + 1):
        g, pairs = count_a302984(n, has_x2_2y2_criterion)
        report("  control (neighbour)", n, g, expect_zero=False)
    if not QUICK:
        g, pairs = count_a302984(A302984[0], lambda m: has_x2_by2(m, 2))
        report("A302984 (plain search)", A302984[0], g, extra="of %d pairs (z, w)" % pairs)

    print("A302982  n = x^2 + 5y^2 + 2^z + 3*2^w  (conjectured for all n > 3)")
    for n in ([558556427] if QUICK else [558556427, 640123571, 683305451]):
        g, pairs = count_a302982(n)
        report("A302982 (plain search)", n, g, extra="of %d pairs (z, w)" % pairs)
    if not QUICK:
        g, pairs = count_a302982(558556426)
        report("  control (neighbour)", 558556426, g, expect_zero=False)

    f31 = lambda c: c * (3 * c + 1) // 2
    f41 = lambda c: c * (4 * c + 1)
    print("A308644 / A308656  n = (smooth number)^2 + two polygonal-type terms  (conjectured for all n > 0)")
    report("A308644 Conj.2 r=3 (3^a 4^b)", 7654105, count_smooth_square(7654105, 3, 4, f31, gen_poly_tester(7, 3)))
    report("A308644 Conj.1(ii) (3^a 5^b)", 108011768, count_smooth_square(108011768, 3, 5, f31, gen_poly_tester(7, 3)))
    report("A308656 Conj.2 f=x(4x+1)", 235644514, count_smooth_square(235644514, 2, 9, f41, gen_poly_tester(3, 1)))
    report("  control (neighbour)", 7654106, count_smooth_square(7654106, 3, 4, f31, gen_poly_tester(7, 3)), expect_zero=False)

    print("A303932  2n = p + 2^k + 3^m with p prime and 11 a quadratic residue (Conj. 1) / non-residue (Conj. 2) mod p")
    for N in (68211369278, 273299806988):
        pairs, nw = check_cert_a303932(N, "certificates/A303932_conj2_N%d.csv" % N, "nr")
        print("PASS certificate (Conj. 2, 2n = %d): all %d pairs (k, m) covered; %d by an explicit square root of 11, the rest by an explicit divisor" % (N, pairs, nw))
        primes, res, nonres = recount_a303932(N)
        report("A303932 Conj.2 (recount)", N, nonres, extra="(%d candidates are prime; 11 is a residue for all of them)" % primes)
    primes, res, nonres = recount_a303932(68211369278 + 2)
    report("  control (neighbour)", 68211369278 + 2, nonres, expect_zero=False)
    for N in (257558146370, 749360555468, 799324864730):
        pairs, nw = check_cert_a303932(N, "certificates/A303932_conj1_N%d.csv" % N, "qr")
        print("PASS certificate (Conj. 1, 2n = %d): all %d pairs (k, m) covered; %d by Euler's criterion (11 a non-residue), the rest by an explicit divisor" % (N, pairs, nw))
        primes, res, nonres = recount_a303932(N)
        report("A303932 Conj.1 (recount)", N, res, extra="(%d candidates are prime; 11 is a non-residue for all of them)" % primes)
    primes, res, nonres = recount_a303932(257558146370 + 2)
    report("  control (neighbour)", 257558146370 + 2, res, expect_zero=False)

    print("\nALL CHECKS PASSED" if OK else "\nSOME CHECK FAILED")
    return 0 if OK else 1


if __name__ == "__main__":
    sys.exit(main())
