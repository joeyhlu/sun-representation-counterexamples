#!/usr/bin/env python3
"""Check a certificate that an odd number N is NOT of the form

        N = p + 2^x + k*2^y,     p prime,  p == c (mod 6),  x, y >= 1.

    python check_cert_a157242.py certificates/A157242_N107576474238337.csv
    python check_cert_a157242.py cert.csv K C N          (another k, class c and N)

The default is the counterexample to Sun's conjecture in OEIS A157242: k = 11, c = 5.

Rows of the certificate:  x,y,p,p_mod_6,divisor
    divisor = 0   p mod 6 != c, so this pair needs no further argument
    divisor = d   1 < d < p and d divides p, so p is not prime

The checker rebuilds the list of ALL pairs (x, y) with x, y >= 1 and 2^x + k*2^y < N by
itself, so a missing, extra or altered row is detected.  It uses integer arithmetic only and
trusts no primality test.  Python 3 standard library.
"""
import sys


def pairs(k, N):
    out = []
    y = 1
    while k * 2 ** y < N:
        x = 1
        while 2 ** x + k * 2 ** y < N:
            out.append((x, y))
            x += 1
        y += 1
    return out


def check(path, k=11, c=5, N=107576474238337):
    assert N % 2 == 1 and c in (1, 5)
    cert = {}
    for line in open(path):
        if line.startswith("#") or line.startswith("x,"):
            continue
        x, y, p, pm, d = map(int, line.split(","))
        assert (x, y) not in cert, "duplicate row"
        cert[(x, y)] = (p, pm, d)
    need = pairs(k, N)
    assert set(need) == set(cert) and len(need) == len(cert), "the certificate does not list exactly the admissible pairs"
    other = composite = 0
    for (x, y) in need:
        p, pm, d = cert[(x, y)]
        assert p == N - 2 ** x - k * 2 ** y and p >= 1 and pm == p % 6
        if p % 6 != c:
            other += 1
        elif p == 1:
            composite += 1                     # 1 is not a prime
        else:
            assert 1 < d < p and p % d == 0, "bad divisor at x=%d, y=%d" % (x, y)
            composite += 1
    print("N = %d: %d pairs (x, y) with x, y >= 1 and 2^x + %d*2^y < N." % (N, len(need), k))
    print("  %d of them give p = N - 2^x - %d*2^y with p mod 6 != %d;" % (other, k, c))
    print("  %d give p == %d (mod 6), each with an explicit proper divisor of p." % (composite, c))
    print("CERTIFICATE VALID: N is not p + 2^x + %d*2^y with p prime, p == %d (mod 6), x, y >= 1." % (k, c))
    return 0


if __name__ == "__main__":
    a = sys.argv[1:]
    if len(a) == 1:
        sys.exit(check(a[0]))
    if len(a) == 4:
        sys.exit(check(a[0], int(a[1]), int(a[2]), int(a[3])))
    print(__doc__)
    sys.exit(2)
