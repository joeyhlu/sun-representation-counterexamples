"""Check a certificate written by make_cert_a303932.py.  Integer arithmetic only.

  python check_cert_a303932.py N cert.csv [qr|nr]

The checker rebuilds the complete list of pairs (k, m), k, m >= 0, with 2^k + 3^m < N, and for
each p = N - 2^k - 3^m accepts exactly one of:
  small    p == 1
  divisor  1 < d < p and p % d == 0                        (p is composite)
  sqrt     r*r % p == 11 % p                               (11 is a square mod p; mode nr)
  nonres   p odd and pow(11, (p-1)//2, p) == p - 1         (11 is not a square mod p if p is prime; mode qr)
  eleven   p == 11                                         (11 == 0 mod 11: neither a non-zero residue nor a non-residue)
In mode nr a valid certificate proves: there is no prime p with N = p + 2^k + 3^m and 11 a
quadratic non-residue mod p.   (If p is prime and 11 == r^2 mod p then 11 is a residue or p = 11.)
"""
import sys

N, path = int(sys.argv[1]), sys.argv[2]
mode = sys.argv[3] if len(sys.argv) > 3 else "nr"
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
assert set(need) == set(cert) and len(need) == len(cert), "certificate does not list exactly the pairs (k, m)"
count = {}
for (k, m) in need:
    p, kind, w = cert[(k, m)]
    assert p == N - 2 ** k - 3 ** m and p >= 1
    if kind == "small":
        assert p == 1
    elif kind == "divisor":
        assert 1 < w < p and p % w == 0
    elif kind == "sqrt":
        assert mode == "nr" and w * w % p == 11 % p
    elif kind == "nonres":
        assert mode == "qr" and p % 2 == 1 and p > 2 and pow(11, (p - 1) // 2, p) == p - 1
    elif kind == "eleven":
        assert p == 11
    else:
        raise AssertionError("unknown row kind")
    count[kind] = count.get(kind, 0) + 1
print("N = %d: %d pairs (k, m): %s" % (N, len(need), ", ".join("%s %d" % kv for kv in sorted(count.items()))))
print("CERTIFICATE VALID (%s): N is not p + 2^k + 3^m with p a prime for which 11 is a quadratic %s."
      % (mode, "non-residue" if mode == "nr" else "residue (p = 2 and p = 11 do not occur)"))
