#!/usr/bin/env python3
"""Decide whether some n^2 - m*4^k is represented by a binary form and, if none is, write the certificate.

    python make_cert_a300510.py out.csv                      (x^2 + y^2, m in {2, 5}, n = 80397893435)
    python make_cert_a300510.py out.csv FORM M1,M2,... n

FORM 1 is x^2 + y^2 (bad primes q == 3 mod 4), FORM 2 is x^2 + 2y^2 (bad primes q == 5, 7 mod 8).

For every m in the list and every k >= 0 with m*4^k < n^2 the number v = n^2 - m*4^k is factored
completely (trial division, Miller-Rabin, Pollard's rho; nothing is shared with the C search program).
v is represented iff every bad prime divides it to an even power.  If no v is represented, the program
writes one row  m,k,v,d  per pair, where d is
  * the odd part of v, when that is == 3 (mod 4) (FORM 1) or == 5, 7 (mod 8) (FORM 2);
  * otherwise q^e for the least bad prime q whose exponent e in v is odd.
In both cases d | v, gcd(d, v/d) = 1 and d lies in a class that proves v is not represented;
`check_cert_a300510.py` verifies exactly that, without any primality test.
Python 3 standard library.
"""
import math
import random
import sys


def is_prime(n):
    """Deterministic Miller-Rabin, valid for n < 3.3e24."""
    if n < 2:
        return False
    small = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41)
    for q in small:
        if n % q == 0:
            return n == q
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for a in small:
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


def rho(n):
    """A proper factor of the odd composite number n."""
    while True:
        c = random.randrange(1, n)
        x = y = random.randrange(n)
        g = 1
        while g == 1:
            x = (x * x + c) % n
            y = (y * y + c) % n
            y = (y * y + c) % n
            g = math.gcd(abs(x - y), n)
        if g != n:
            return g


def factor(n):
    """dict prime -> exponent"""
    out = {}
    for p in range(2, 2000):
        while n % p == 0:
            out[p] = out.get(p, 0) + 1
            n //= p
    stack = [n] if n > 1 else []
    while stack:
        m = stack.pop()
        if is_prime(m):
            out[m] = out.get(m, 0) + 1
        else:
            d = rho(m)
            stack += [d, m // d]
    return out


def main(out, form=1, mults=(2, 5), n=80397893435):
    bad = (lambda q: q % 4 == 3) if form == 1 else (lambda q: q % 8 in (5, 7))
    rows, reps = [], 0
    for m in mults:
        k = 0
        while m * 4 ** k < n * n:
            v = n * n - m * 4 ** k
            fac = factor(v)
            assert math.prod(p ** e for p, e in fac.items()) == v
            odd_bad = sorted(p for p, e in fac.items() if p != 2 and bad(p) and e % 2 == 1)
            if not odd_bad:
                reps += 1
                print("  m = %d, k = %d: v = %d IS represented" % (m, k, v))
            else:
                u = v
                while u % 2 == 0:
                    u //= 2
                d = u if bad(u) else odd_bad[0] ** fac[odd_bad[0]]
                assert v % d == 0 and math.gcd(d, v // d) == 1 and bad(d)
                rows.append((m, k, v, d))
            k += 1
    print("n = %d: %d pairs (m, k), %d represented" % (n, len(rows) + reps, reps))
    if reps:
        print("NOT a counterexample.")
        return 1
    name = "x^2 + y^2" if form == 1 else "x^2 + 2y^2"
    with open(out, "w", newline="\n") as f:
        f.write("# n=%d: no n^2 - m*4^k with m in %s, k >= 0, m*4^k < n^2 is of the form %s\n"
                % (n, "{" + ", ".join(map(str, mults)) + "}", name))
        f.write("m,k,v,d\n")
        for r in rows:
            f.write("%d,%d,%d,%d\n" % r)
    print("NO REPRESENTATION; certificate written to", out)
    return 0


if __name__ == "__main__":
    a = sys.argv[1:]
    if len(a) == 1:
        sys.exit(main(a[0]))
    if len(a) == 4:
        sys.exit(main(a[0], int(a[1]), tuple(int(t) for t in a[2].split(",")), int(a[3])))
    print(__doc__)
    sys.exit(2)
