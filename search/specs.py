"""Term lists (spec files) and a brute-force reference for mixsum2.c.

  python specs.py write outdir          write the .spec files
  python specs.py zeros NAME N          brute force: all n <= N (of the spec's parity) with no representation
  python specs.py count NAME N          first 60 values of the representation count
  python specs.py check NAME n          list the representations of one n

Conjectures (all by Zhi-Wei Sun, OEIS, 2018):
  a303932qr : every even 2n > 2 is p + 2^k + 3^m, p prime with 11 a quadratic residue mod p     (A303932, Conj. 1)
  a303932nr : every even 2n > 4 is p + 2^k + 3^m, p prime with 11 a quadratic non-residue mod p (A303932, Conj. 2)
  a302984   : every n > 5 is x^2 + 2y^2 + 2^z + 5*2^w                                            (A302984)
  a302983   : every n > 3 is x^2 + 2y^2 + 2^z + 3*2^w                                            (A302983; used as a control)
"""
import sys
from math import isqrt

LIMIT = 1 << 62


def pw(base, mult=1, limit=LIMIT):
    out, v = [], 1
    while mult * v < limit:
        out.append(mult * v)
        v *= base
    return out


def legendre_11(p):
    """1 if 11 is a non-zero square mod the odd prime p, -1 if a non-square, 0 if p = 11."""
    if p == 11:
        return 0
    return 1 if pow(11, (p - 1) // 2, p) == 1 else -1


def classes_mod_44(sign):
    """Residues r mod 44 (coprime to 44) such that primes p == r have (11/p) = sign."""
    qr11 = {1, 3, 4, 5, 9}
    out = []
    for r in range(44):
        if r % 2 == 0 or r % 11 == 0:
            continue
        s = (1 if (r % 11) in qr11 else -1) * (1 if r % 4 == 1 else -1)     # quadratic reciprocity
        if s == sign:
            out.append(r)
    return out


SPECS = {
    # name: (base, pmod, allowed residues, extra values, parity, A, B, text)
    "a303932qr": ("prime", 44, classes_mod_44(1), [2, 11], "even", pw(2), pw(3),
                  "2n = p + 2^k + 3^m, (11/p) = 1; p = 2 and p = 11 are also admitted, so a number found here fails under every reading"),
    "a303932nr": ("prime", 44, classes_mod_44(-1), [], "even", pw(2), pw(3),
                  "2n = p + 2^k + 3^m, (11/p) = -1"),
    "a302984": ("x2p2y2", 0, [], [], "all", pw(2), pw(2, 5), "n = x^2 + 2y^2 + 2^z + 5*2^w"),
    "a302983": ("x2p2y2", 0, [], [], "all", pw(2), pw(2, 3), "n = x^2 + 2y^2 + 2^z + 3*2^w"),
}


def write(outdir):
    import os
    os.makedirs(outdir, exist_ok=True)
    for name, (base, pmod, allow, extra, parity, A, B, text) in SPECS.items():
        with open(os.path.join(outdir, name + ".spec"), "w", newline="\n") as f:
            f.write("name %s\nbase %s\nparity %s\n" % (name, base, parity))
            if pmod:
                f.write("pmod %d %d %s\n" % (pmod, len(allow), " ".join(map(str, allow))))
            if extra:
                f.write("extra %d %s\n" % (len(extra), " ".join(map(str, extra))))
            f.write("A %d %s\n" % (len(A), " ".join(map(str, A))))
            f.write("B %d %s\n" % (len(B), " ".join(map(str, B))))
        print("%-10s %s" % (name, text))


def base_table(name, N):
    base, pmod, allow, extra, parity, A, B, text = SPECS[name]
    t = bytearray(N + 1)
    if base == "prime":
        s = bytearray([1]) * (N + 1)
        s[0:2] = b"\x00\x00"
        for i in range(2, isqrt(N) + 1):
            if s[i]:
                s[i * i::i] = bytearray(len(s[i * i::i]))
        al = set(allow)
        for p in range(3, N + 1, 2):
            if s[p] and (not pmod or p % pmod in al):
                t[p] = 1
        for v in extra:
            if v <= N:
                t[v] = 1
    else:
        x = 0
        while x * x <= N:
            y = 0
            while x * x + 2 * y * y <= N:
                t[x * x + 2 * y * y] = 1
                y += 1
            x += 1
    return t


def counts(name, N):
    base, pmod, allow, extra, parity, A, B, text = SPECS[name]
    t = base_table(name, N)
    cnt = [0] * (N + 1)
    for a in A:
        if a > N:
            break
        for b in B:
            if a + b > N:
                break
            s = a + b
            for m in range(0, N - s + 1):
                if t[m]:
                    cnt[m + s] += 1
    return cnt


def main():
    mode = sys.argv[1]
    if mode == "write":
        write(sys.argv[2] if len(sys.argv) > 2 else "specs2")
        return
    name = sys.argv[2]
    base, pmod, allow, extra, parity, A, B, text = SPECS[name]
    if mode in ("zeros", "count"):
        N = int(float(sys.argv[3]))
        cnt = counts(name, N)
        if mode == "count":
            print(name, text)
            print(",".join(str(c) for c in cnt[1:61]))
        par = {"all": (0, 1), "odd": (1,), "even": (0,)}[parity]
        print("no representation, n <= %d: %s" % (N, [n for n in range(1, N + 1) if cnt[n] == 0 and n % 2 in par]))
    elif mode == "check":
        n = int(sys.argv[3])
        # exact membership without tables
        def member(m):
            if base == "prime":
                if m in extra:
                    return True
                if m < 3 or m % 2 == 0 or (pmod and m % pmod not in allow):
                    return False
                return all(m % q for q in range(3, isqrt(m) + 1, 2))
            y = 0
            while 2 * y * y <= m:
                r = isqrt(m - 2 * y * y)
                if r * r == m - 2 * y * y:
                    return True
                y += 1
            return False
        good = 0
        for a in A:
            if a > n:
                break
            for b in B:
                if a + b > n:
                    break
                if member(n - a - b):
                    good += 1
                    if good <= 5:
                        print("  REPRESENTATION: %d = %d + %d + %d" % (n, n - a - b, a, b))
        print("%s n=%d: %d representations" % (name, n, good))


if __name__ == "__main__":
    main()
