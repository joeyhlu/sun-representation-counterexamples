# Search programs

Both programs need gcc with OpenMP on a 64-bit machine.

```
gcc -O3 -march=native -fopenmp -o mixsum2 mixsum2.c -lm
gcc -O3 -march=native -fopenmp -o direct_scan direct_scan.c -lm
python specs.py write specs          # writes the term lists used by mixsum2
```

## mixsum2: bit-parallel exhaustive search

```
./mixsum2 specs/a302984.spec   1 1e11 out_a302984      # x^2 + 2y^2 + 2^z + 5*2^w
./mixsum2 specs/a303932nr.spec 1 1e12 out_a303932nr    # p + 2^k + 3^m, 11 a non-residue mod p
./mixsum2 specs/a303932qr.spec 1 1e12 out_a303932qr    # p + 2^k + 3^m, 11 a residue mod p
```

Each run writes `<out>.found`, one line per n in the range with no representation (the first
few lines are the trivial small n). A run can be interrupted and restarted; it resumes from
`<out>.state`.

How it works. Bit i of a bitset stands for the integer i. With P the indicator of the base
set (primes in the allowed residue classes, or numbers x² + 2y²), the program forms
Q = OR of (P shifted by a) over the terms a of the first list, then C = OR of (Q shifted by b)
over the second list, using only terms up to 2^28 at this stage. Every n with C[n] = 0 is then
tested exactly against all pairs (a, b) of any size, by a deterministic Miller-Rabin test for
primes or by factoring for x² + 2y². The bit stage only ever sets bits for true
representations, so it cannot hide a number that has none.

For the base set x² + 2y² the sieve marks the numbers in which some prime ≡ 5, 7 (mod 8)
occurs to an odd power; a number is in the set when it is unmarked and its odd part is ≡ 1 or
3 (mod 8).

In `a303932qr.spec` the primes 2 and 11 are also admitted, so a number reported by that run
fails Conjecture 1 whichever way "quadratic residue" is read for those two primes.

## direct_scan: independent check by table look-up

```
./direct_scan 2 1e10      # A302984
./direct_scan 1 1e9       # A302982
./direct_scan 3 3e8       # A308644 Conjecture 2, r = 3
./direct_scan 4 3e8       # A308644 Conjecture 1, second form
./direct_scan 5 3e8       # A308656 Conjecture 2, f(x) = x(4x+1)
```

It marks every sum of two quadratic values up to N in a bit table with a double loop and
then looks up n − t for every offset t. It shares no code and no number-theoretic criterion
with `mixsum2`. Memory use is N/8 bytes.

## specs.py

Defines the term lists and has a slow brute-force mode that was used to test `mixsum2` on
small ranges:

```
python specs.py zeros a302984 40000
python specs.py count a303932qr 200
```
