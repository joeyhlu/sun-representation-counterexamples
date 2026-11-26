#!/usr/bin/env python3
"""Write the certificate that an odd N is not p + 2^x + k*2^y with p prime, p == c (mod 6), x, y >= 1.

    python make_cert_a157242.py out.csv                 (k = 11, c = 5, N = 107576474238337)
    python make_cert_a157242.py out.csv K C N

For every pair (x, y) with x, y >= 1 and 2^x + k*2^y < N the row gives p = N - 2^x - k*2^y,
p mod 6 and, when p == c (mod 6), the least prime factor of p (trial division up to 10^5, then
Pollard's rho).  The program stops with a message if some p in the class is prime, that is, if N
has a representation after all.  `check_cert_a157242.py` verifies the result without any
primality test.  Python 3 standard library.
"""
import math
import random
import sys


def is_prime(n):
    """Deterministic Miller-Rabin, valid for n < 3.3e24."""
    if n < 2:
        return False
    small = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37)
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
    """A non-trivial factor of the odd composite number n (Pollard's rho, Brent's variant)."""
    while True:
        y, c, m = random.randrange(1, n), random.randrange(1, n), 128
        g = r = q = 1
        while g == 1:
            x = y
            for _ in range(r):
                y = (y * y + c) % n
            k = 0
            while k < r and g == 1:
                ys = y
                for _ in range(min(m, r - k)):
                    y = (y * y + c) % n
                    q = q * abs(x - y) % n
                g = math.gcd(q, n)
                k += m
            r *= 2
        if g == n:
            g = 1
            while g == 1:
                ys = (ys * ys + c) % n
                g = math.gcd(abs(x - ys), n)
        if g != n:
            return g


def least_prime_factor(n):
    for q in (2, 3, 5):
        if n % q == 0:
            return q
    q = 7
    while q < 100000 and q * q <= n:
        if n % q == 0:
            return q
        q += 2
    if q * q > n or is_prime(n):
        return n
    best, stack = n, [n]
    while stack:
        m = stack.pop()
        if m == 1:
            continue
        if is_prime(m):
            best = min(best, m)
            continue
        d = rho(m)
        stack += [d, m // d]
    return best


def main(out, k=11, c=5, N=107576474238337):
    assert N % 2 == 1 and c in (1, 5)
    rows = []
    y = 1
    while k * 2 ** y < N:
        x = 1
        while 2 ** x + k * 2 ** y < N:
            p = N - 2 ** x - k * 2 ** y
            if p % 6 != c:
                rows.append((x, y, p, p % 6, 0))
            else:
                d = least_prime_factor(p) if p > 1 else 0
                if p > 1 and d == p:
                    print("NOT a counterexample: %d = %d + 2^%d + %d*2^%d with %d prime, == %d (mod 6)"
                          % (N, p, x, k, y, p, c))
                    return 1
                rows.append((x, y, p, p % 6, d))
            x += 1
        y += 1
    rows.sort(key=lambda r: (r[1], r[0]))
    with open(out, "w", newline="\n") as f:
        f.write("# N=%d is not p + 2^x + %d*2^y with p prime, p == %d (mod 6), x,y >= 1\n" % (N, k, c))
        f.write("x,y,p,p_mod_6,divisor\n")
        for r in rows:
            f.write(",".join(map(str, r)) + "\n")
    inclass = [r for r in rows if r[3] == c]
    print("N = %d: %d pairs; %d with p == %d (mod 6), all composite; largest least prime factor %d"
          % (N, len(rows), len(inclass), c, max(r[4] for r in inclass)))
    return 0


if __name__ == "__main__":
    a = sys.argv[1:]
    if len(a) == 1:
        sys.exit(main(a[0]))
    if len(a) == 4:
        sys.exit(main(a[0], int(a[1]), int(a[2]), int(a[3])))
    print(__doc__)
    sys.exit(2)
