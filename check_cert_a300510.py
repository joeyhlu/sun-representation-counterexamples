#!/usr/bin/env python3
"""Check a certificate that no number n^2 - m*4^k is represented by a binary form.  No primality test.

    python check_cert_a300510.py certificate.csv                     (x^2 + y^2, m in {2, 5}, n = 80397893435)
    python check_cert_a300510.py certificate.csv FORM M1,M2,... n

FORM 1 is x^2 + y^2, FORM 2 is x^2 + 2y^2.

The program builds the list of all pairs (m, k) with k >= 0 and m*4^k < n^2 itself and requires one
row  m,k,v,d  for each of them, with

    v = n^2 - m*4^k,    d odd,    d | v,    gcd(d, v/d) = 1,
    d == 3 (mod 4)           for FORM 1,
    d == 5 or 7 (mod 8)      for FORM 2.

Why such a d proves that v is not represented.  Call a prime q bad if q == 3 (mod 4) (FORM 1), resp.
q == 5 or 7 (mod 8) (FORM 2).  The odd numbers all of whose bad prime factors occur to an even power
are == 1 (mod 4), resp. == 1 or 3 (mod 8), so d has a bad prime factor q of odd multiplicity.  Since d
is coprime to v/d, q divides v with the same odd multiplicity.  A number represented by the form has
every bad prime to an even power: if q | x^2 + y^2 and -1 is not a square modulo q, then q | x and
q | y (and likewise with -2 for x^2 + 2y^2), so q^2 can be divided out.

Python 3 standard library.
"""
import csv
import math
import sys


def check(path, form=1, mults=(2, 5), n=80397893435):
    need = {}
    for m in mults:
        k = 0
        while m * 4 ** k < n * n:
            need[(m, k)] = n * n - m * 4 ** k
            k += 1
    seen = set()
    for row in csv.reader(open(path)):
        if not row or not row[0].isdigit():
            continue
        m, k, v, d = (int(t) for t in row)
        assert (m, k) in need and (m, k) not in seen, "unexpected or repeated pair (%d, %d)" % (m, k)
        assert v == need[(m, k)], "wrong v for the pair (%d, %d)" % (m, k)
        assert d > 1 and d % 2 == 1 and v % d == 0, "d does not divide v for the pair (%d, %d)" % (m, k)
        assert math.gcd(d, v // d) == 1, "d is not coprime to v/d for the pair (%d, %d)" % (m, k)
        assert (d % 4 == 3) if form == 1 else (d % 8 in (5, 7)), "d is in an allowed class for the pair (%d, %d)" % (m, k)
        seen.add((m, k))
    assert seen == set(need), "missing pairs: %s" % sorted(set(need) - seen)[:5]
    return len(need)


if __name__ == "__main__":
    a = sys.argv[1:]
    if len(a) == 1:
        form, mults, n = 1, (2, 5), 80397893435
    elif len(a) == 4:
        form, mults, n = int(a[1]), tuple(int(t) for t in a[2].split(",")), int(a[3])
    else:
        print(__doc__)
        sys.exit(2)
    pairs = check(a[0], form, mults, n)
    print("n = %d: %d pairs (m, k) with m in %s, k >= 0 and m*4^k < n^2, each with a divisor d as required."
          % (n, pairs, list(mults)))
    print("CERTIFICATE VALID: no n^2 - m*4^k is of the form %s." % ("x^2 + y^2" if form == 1 else "x^2 + 2y^2"))
