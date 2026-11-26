"""Write the certificate for an A303932 counterexample.

  python make_cert_a303932.py N out.csv [qr|nr]

Conjecture 2 of A303932 (mode nr, default): every even number 2n > 4 is p + 2^k + 3^m with p a
prime such that 11 is a quadratic NON-residue mod p (k, m >= 0).
Conjecture 1 (mode qr): the same with 11 a quadratic residue mod p.

For every pair (k, m) with p = N - 2^k - 3^m >= 1 the certificate gives one of
  small        p is 1 (not a prime)
  divisor d    1 < d < p and d | p                      -> p is not prime
  sqrt r       r*r == 11 (mod p)                        -> 11 is a square mod p   (used in mode nr)
  nonres       p is an odd number with 11^((p-1)/2) == -1 (mod p)
               -> if p is prime, 11 is a non-residue mod p (Euler's criterion)     (used in mode qr)
check_cert_a303932.py verifies these rows with integer arithmetic only (no primality test).
"""
import math
import random
import sys


def is_prime(n):
    if n < 2:
        return False
    for q in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % q == 0:
            return n == q
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
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


def small_factor(n):
    for q in (2, 3, 5):
        if n % q == 0:
            return q
    q = 7
    while q < 200000 and q * q <= n:
        if n % q == 0:
            return q
        q += 2
    while True:                                  # Pollard rho
        c = random.randrange(1, n)
        x = y = random.randrange(2, n)
        d = 1
        while d == 1:
            x = (x * x + c) % n
            y = (y * y + c) % n
            y = (y * y + c) % n
            d = math.gcd(abs(x - y), n)
        if d != n:
            return d


def sqrt_mod(a, p):
    """Tonelli-Shanks: r with r*r == a (mod p), p an odd prime, a a non-zero square mod p."""
    a %= p
    if p % 4 == 3:
        return pow(a, (p + 1) // 4, p)
    q, s = p - 1, 0
    while q % 2 == 0:
        q //= 2
        s += 1
    z = 2
    while pow(z, (p - 1) // 2, p) != p - 1:
        z += 1
    m, c, t, r = s, pow(z, q, p), pow(a, q, p), pow(a, (q + 1) // 2, p)
    while t != 1:
        i, t2 = 0, t
        while t2 != 1:
            t2 = t2 * t2 % p
            i += 1
        b = pow(c, 1 << (m - i - 1), p)
        m, c = i, b * b % p
        t, r = t * c % p, r * b % p
    return r


def main():
    N, out = int(sys.argv[1]), sys.argv[2]
    mode = sys.argv[3] if len(sys.argv) > 3 else "nr"
    rows = []
    k = 0
    while 2 ** k < N:
        m = 0
        while 2 ** k + 3 ** m < N:
            p = N - 2 ** k - 3 ** m
            if p == 1:
                rows.append((k, m, p, "small", 0))
            elif not is_prime(p):
                rows.append((k, m, p, "divisor", small_factor(p)))
            else:
                if p == 2:
                    sys.exit("p = 2 occurs (k=%d, m=%d): 11 is a square mod 2; treat by hand" % (k, m))
                if p == 11:
                    kind = "eleven"
                    rows.append((k, m, p, kind, 0))
                    m += 1
                    continue
                res = pow(11, (p - 1) // 2, p) == 1
                if mode == "nr":
                    if not res:
                        sys.exit("NOT a counterexample: p = %d is prime and 11 is a non-residue (k=%d, m=%d)" % (p, k, m))
                    r = sqrt_mod(11, p)
                    assert r * r % p == 11 % p
                    rows.append((k, m, p, "sqrt", r))
                else:
                    if res:
                        sys.exit("NOT a counterexample: p = %d is prime and 11 is a residue (k=%d, m=%d)" % (p, k, m))
                    rows.append((k, m, p, "nonres", 0))
            m += 1
        k += 1
    with open(out, "w", newline="\n") as f:
        f.write("# N=%d, mode=%s: p = N - 2^k - 3^m for all k, m >= 0 with p >= 1\n" % (N, mode))
        f.write("k,m,p,kind,witness\n")
        for r in rows:
            f.write(",".join(map(str, r)) + "\n")
    kinds = {}
    for r in rows:
        kinds[r[3]] = kinds.get(r[3], 0) + 1
    print("N = %d: %d pairs; %s" % (N, len(rows), ", ".join("%s: %d" % kv for kv in sorted(kinds.items()))))


if __name__ == "__main__":
    main()
