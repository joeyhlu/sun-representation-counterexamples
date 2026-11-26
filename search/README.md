# Search programs

All programs need gcc with OpenMP on a 64-bit machine.

```
gcc -O3 -march=native -fopenmp -o mixsum2 mixsum2.c -lm
gcc -O3 -march=native -fopenmp -o direct_scan direct_scan.c -lm
gcc -O3 -march=native -fopenmp -o sunmix_multi sunmix_multi.c -lm
gcc -O3 -march=native -fopenmp -o suntarget2 suntarget2.c -lm
gcc -O2 -march=native -fopenmp -o sqsparse sqsparse.c -lm
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

## sunmix_multi: exhaustive search for odd N ≠ p + 2^x + k·2^y

```
./sunmix_multi 7,11 5 1 1e13 out_c5                                 # k = 7 and 11, p ≡ 5 (mod 6), all odd N < 10^13
./sunmix_multi 7,11 1 1 1e13 out_c1                                 # the same for p ≡ 1 (mod 6)
./sunmix_multi 11 5 107570000000000 107580000000000 out_window      # a window around 107576474238337
```

Arguments: the list of k, the class of p modulo 6 (1 or 5; 0 for no condition), the range
[Nlo, Nhi), an output name, and optionally the number of threads. `<out>.found` gets one line
"k N" for every odd N in the range with no representation, including the trivial small N. A
run can be interrupted and restarted; it resumes from `<out>.state`.

How it works. Bit j stands for the odd number 2j + 1. With P the primes in the class,
Q = OR of (P shifted by 2^(x−1)) for x ≤ 26 and C = OR of (Q shifted by k·2^(y−1)) for
k·2^(y−1) ≤ 2^25. Every N with C = 0 is then tested exactly against all pairs (x, y), with a
deterministic Miller-Rabin test. As in `mixsum2`, the bit stage can only clear numbers that
do have a representation.

## suntarget2: targeted search for the same form

```
./suntarget2 11 5 1e13 1.024e16 3e11 out_t11_c5 16 1e-11 1024       # found 107576474238337
./suntarget2 11 1 1e13 1.024e16 3e11 out_t11_c1 16 1e-11 1024       # nothing
```

Arguments: k, the class of p modulo 6, the range [Xlo, Xhi), the number of N to look at per
octave, an output name, threads, the smallest estimated probability of being an exception
for which an N is still tested exactly, and the largest prime used for that estimate. A
negative number of N per octave prints the plan without running it. The same two commands
with the range 1.024e16 to 4.5e18 complete the runs reported in the main README.

This search is not exhaustive. For every octave it ranks the residue classes of N modulo
9699690 by how many pairs (x, y) can still give a p that is in the right class modulo 6 and
is not divisible by a prime up to 19, and it only looks at N in the best classes. Among those
it tests exactly the N for which so few candidates survive the primes 23 to 1021 that the
estimated probability of all of them being composite reaches the threshold. A number is
reported only after the complete test over all pairs (x, y), so what it reports is correct,
but a number it does not report may still be an exception.

## sqsparse: odd n with no n² − m·4^k represented by x² + y² or x² + 2y²

```
./sqsparse scan 1 2,5 6e9 1e11 12                   # A300510: found 80397893435 (about 15 minutes on 16 threads)
./sqsparse one  1 2,5 80397893435                   # every pair (m, k) with a witness
./sqsparse one  1 1,2,5 80397893435                 # the same for A300448, Conjecture 1
./sqsparse scan 1 2,29 3 7e7 200                    # every odd n decided: only 51883659
./sqsparse scan 1 1,5,17 3 239694785 200            # A300448 Conjecture 2, every odd n decided: only 239694783
./sqsparse scan 2 3,6 3 6e9 14                      # Square Conjecture (3*2^k = 3*4^j or 6*4^j): only 5884015571
./sqsparse scan 1 1,5,65 2e7 1e11 14                # A300448 Conjecture 3: found 68136270467 and 80570407837
./sqsparse cls  2 2,3 1e11 1e12 13 7,13,23,31,37 -14  # A301452 in the best residue classes: found 386749710995 and 985224034231
./sqsparse one  2 2,3 386749710995                  # every pair (m, k) with a witness
```

Arguments of `scan`: the form (1 for x² + y², 2 for x² + 2y²), the list of multipliers m, the
range [NLO, NHI) of n, the limit T, and optionally the number of threads and the sieve bound
(default 65536). Only odd n are examined: if an even n has no representation, then neither
has n/2.

How it works. Call a prime bad if it must divide a represented number to an even power
(q ≡ 3 mod 4 for x² + y², q ≡ 5, 7 mod 8 for x² + 2y²). Stage 1 sieves blocks of odd n with the
bad primes up to the bound: q divides n² − m·4^k exactly when n ≡ ±s·2^k (mod q) with s² ≡ m,
so each residue class of n modulo q strikes a fixed set of pairs (m, k), kept as a bit mask.
Stage 2 takes the n with at most T pairs not struck and decides those pairs exactly (trial
division, strong probable prime tests, Pollard's rho) until one is represented. If none is,
all pairs of that n are decided exactly, because the sieve also strikes a pair when q² divides
the number. Only an n that fails this complete test is printed (`NOREP n`).

With T at least the number of pairs (for example 200) every odd n is decided, and the run is
exhaustive. With a smaller T the run is not exhaustive: an n with more than T pairs left is
skipped. `cls` restricts the sieve further, to the residue classes of n modulo a product of
small bad primes in which those primes alone leave the fewest pairs; `stats` prints, for the
n with at most T pairs left, how many of those pairs are represented.
