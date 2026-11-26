# Counterexamples to thirteen conjectures of Zhi-Wei Sun on representations of integers

Zhi-Wei Sun has posted in the OEIS many conjectures of the form "every positive integer is a
sum of a few sparse terms". Thirteen of those statements are false: one from 2009 on
p + 2^x + 11·2^y, five from March 2018 on n² minus a small multiple of a power of 4 (his "2-5
conjecture", his "2-3 conjecture" and three related statements), and seven from 2018 and 2019.
Each counterexample is one integer, and `python verify.py --quick` re-checks all of them in
about a minute using only the Python standard library.

| OEIS | Conjecture | Checked before up to | Counterexample |
|---|---|---|---|
| [A157242](https://oeis.org/A157242) | every odd N > 28 other than 33, 39, 45, 171, 247 is p + 2^x + 11·2^y with x, y ≥ 1 and p a prime ≡ 5 (mod 6) | 10^13 (D. S. McNeil, 2009, as reported by Sun) | **107576474238337** |
| [A300510](https://oeis.org/A300510), the "2-5 conjecture" | for every n > 1 some n² − 2·4^k or n² − 5·4^k with k ≥ 0 is a sum of two squares | n < 6·10^9 (G. Resta, 2019) | **80397893435** |
| [A300448](https://oeis.org/A300448), Conjecture 1 | for every n > 0 some n² − m·4^k with m ∈ {1, 2, 5} and k ≥ 0 is a sum of two squares | n ≤ 2·10^7 (Sun, 2018) | **80397893435** |
| [A300448](https://oeis.org/A300448), Conjecture 2 | for every n > 0 some n² − m·4^k with m ∈ {1, 5, 17} and k ≥ 0 is a sum of two squares | n ≤ 2·10^7 (Sun, 2018) | **239694783** |
| [A300448](https://oeis.org/A300448), Conjecture 3 | for every n > 0 some n² − m·4^k with m ∈ {1, 5, 65} and k ≥ 0 is a sum of two squares | n ≤ 2·10^7 (Sun, 2018) | **68136270467** |
| [A301452](https://oeis.org/A301452), the "2-3 conjecture" | for every n > 1 some n² − 2·4^k or n² − 3·4^k with k ≥ 0 is of the form x² + 2y² | n ≤ 5·10^7 (Sun, 2018) | **386749710995** |
| [A302984](https://oeis.org/A302984) | every n > 5 is x² + 2y² + 2^z + 5·2^w | 10^9 | **4375339279** |
| [A302982](https://oeis.org/A302982) | every n > 3 is x² + 5y² + 2^z + 3·2^w | 2·10^8 | **558556427** |
| [A303932](https://oeis.org/A303932), Conjecture 1 | every even number 2n > 2 is p + 2^k + 3^m with p a prime such that 11 is a quadratic residue mod p | n ≤ 5·10^8 | **2n = 257558146370** |
| [A303932](https://oeis.org/A303932), Conjecture 2 | every even number 2n > 4 is p + 2^k + 3^m with p a prime such that 11 is a quadratic non-residue mod p | n ≤ 5·10^8 | **2n = 68211369278** |
| [A308644](https://oeis.org/A308644), Conjecture 2 (r = 3) | every n > 0 is (3^a·4^b)² + c(3c+1)/2 + d(7d+3)/2 | 10^6 | **7654105** |
| [A308644](https://oeis.org/A308644), Conjecture 1 (second assertion) | every n > 0 is (3^a·5^b)² + c(3c+1)/2 + d(7d+3)/2 | 10^6 | **108011768** |
| [A308656](https://oeis.org/A308656), Conjecture 2 (f(x) = x(4x+1)) | every n > 0 is (2^a·9^b)² + c(4c+1) + d(3d+1)/2 | 10^6 | **235644514** |

In the last seven rows, x, y, z, w, a, b, k, m are non-negative integers and c, d are
arbitrary integers, as in the OEIS entries. Those seven counterexamples are the smallest ones,
and so are 107576474238337 and 239694783. For the 2-5 conjecture, the 2-3 conjecture and
Conjectures 1 and 3 of A300448 the smallest counterexample is not known. See the next two
sections.

## The 2009 conjecture on p + 2^x + 11·2^y (A157242)

On 24 February 2009 Sun conjectured that, except for 33, 39, 45, 171 and 247, every odd
integer greater than 28 can be written as p + 2^x + 11·2^y, where x and y are positive
integers and p is a prime congruent to 5 modulo 6
([his post to the Number Theory List](https://listserv.nodak.edu/cgi-bin/wa.exe?A2=NMBRTHRY;20d0ba9.0902);
OEIS A157242, where the entry carries the keyword "nice"). Sun checked this up to 5·10^7 and
Qing-Hu Hou up to 1.5·10^8. D. S. McNeil checked it up to 10^12, which is the bound the OEIS
entry gives, and then up to 10^13
([Sun's post of 28 February 2009](https://listserv.nodak.edu/cgi-bin/wa.exe?A2=NMBRTHRY;e19bbf22.0902)).

**The odd number N = 107576474238337 = 13 · 43 · 2129 · 7411 · 12197 has no such
representation.** In the indexing of the OEIS entry, a(53788237119169) = 0.

* There are 1974 pairs (x, y) with x, y ≥ 1 and 2^x + 11·2^y < N.
* N ≡ 1 (mod 3), so p = N − 2^x − 11·2^y is congruent to 5 modulo 6 exactly when x is even
  and y is odd. That leaves 504 pairs.
* All 504 values of p are composite. 363 of them have a prime factor at most 19, 90 have
  their least prime factor between 23 and 1000, and 51 have no prime factor below 1000. The
  largest least prime factor is 6704531, for x = 16 and y = 17.

The counterexample does not come from a covering congruence: 51 of the candidates have no
prime factor below 1000 and are composite for no visible reason. Bjorn Poonen had predicted
on heuristic grounds in March 2009 that for every k infinitely many odd numbers are not of
the form p + 2^x + k·2^y, even with no condition on p modulo 6
([his post](https://listserv.nodak.edu/cgi-bin/wa.exe?A2=NMBRTHRY;bb043890.0903)), so a
failure was to be expected somewhere. We found no record of an actual counterexample to the
conjecture; see "Prior work".

How N relates to smaller numbers:

* Every odd N below 10^13 was tested (`search/sunmix_multi.c`). The only odd numbers in
  (28, 10^13) with no representation are Sun's five. This repeats McNeil's verification.
* The number 107576474238337 came from a targeted search (`search/suntarget2.c`) that looks
  only at the residue classes in which a counterexample is most likely. That search cannot
  show that there is none in between.
* The exhaustive program was then run on every odd number from 10^13 to 1.0758·10^14. The
  only one with no representation is N itself. So **107576474238337 is the smallest
  counterexample**: the odd numbers in (28, 1.0758·10^14) that are not of the form
  p + 2^x + 11·2^y are 33, 39, 45, 171, 247 and 107576474238337.

The companion conjecture with p ≡ 1 (mod 6) ([A157237](https://oeis.org/A157237): every odd
N > 30 other than 35, 41, 47, 101, 167, 2021, 119171) is still open. Below 10^13 its only
exceptions are those seven, and the targeted search, run up to 4.5·10^18, found no
counterexample to it.

## The 2018 conjectures on n² minus a multiple of a power of 4 (A300510, A300448, A301452)

### The "2-5 conjecture" (A300510) and Conjecture 1 of A300448

On 7 March 2018 Sun conjectured that for every integer n > 1 there is a k ≥ 0 such that
n² − 2·4^k or n² − 5·4^k is a sum of two squares (OEIS A300510; in A301452 he calls this his
"2-5 conjecture"). It says that n² = (2^k)² + (2^k·m)² + x² + y² with m = 1 or 2, a special way
of writing n² as a sum of four squares. The OEIS entry records three checks: by Sun for
2 ≤ n ≤ 5·10^7, by Qing-Hu Hou up to 10^9, and by Giovanni Resta for 1 < n < 6·10^9 (June 2019).
In April 2020 Sun asked on MathOverflow whether a counterexample can be found
([question 356905](https://mathoverflow.net/questions/356905)). The one answer there listed
failures for the two coefficients separately, which is not a counterexample, as Sun pointed out
in a comment.

**For n = 80397893435 = 5 · 1163 · 1237 · 11177 there is no such k.** In the notation of the OEIS
entry, a(80397893435) = 0.

* There are 72 numbers to consider: n² − 2·4^k and n² − 5·4^k for k = 0, …, 35.
* A positive integer is a sum of two squares exactly when every prime ≡ 3 (mod 4) divides it to
  an even power. Each of the 72 numbers has a prime ≡ 3 (mod 4) to an odd power.
* For comparison, 19 of the 72 numbers are sums of two squares for n − 2, and 18 for n + 2. The
  OEIS sequence, which counts the pairs x ≤ y as well, has a(n − 2) = 47 and a(n + 2) = 39.

The same n refutes the weaker Conjecture 1 of [A300448](https://oeis.org/A300448), which also
allows n² − 4^k. That adds the 37 numbers n² − 4^k with k = 0, …, 36, and none of the 109 numbers
is a sum of two squares. Sun reports that he checked the conjectures of A300448 for n ≤ 2·10^7.

How n relates to smaller numbers:

* If an even n has no representation, then neither has n/2, because
  n² − m·4^(k+1) = 4·((n/2)² − m·4^k). So the smallest counterexample is odd.
* Below 6·10^9 there is no counterexample to A300510, by Resta's check.
* Between 6·10^9 and 10^11 the search here sieved every odd n, but it decided exactly only the n
  with at most 12 pairs (m, k) left after the sieve (see "How the numbers were found"). Among
  those, 80397893435 is the only one with no representation. The other odd n were not decided,
  so **80397893435 is not known to be the smallest counterexample.**

### Conjectures 2 and 3 of A300448

**Conjecture 2 of A300448** uses m ∈ {1, 5, 17}: every positive square n² should be
4^k·m + x² + y² for one of these m. (Sun first stated it for {1, 5, 11} as well and withdrew
that case the same day, 6 March 2018, because it fails at n = 13018247; the edit history of the
entry shows this.) It fails at n = 239694783 = 3 · 47 · 179 · 9497. None of the 81 numbers
n² − 4^k (k ≤ 27), n² − 5·4^k (k ≤ 26) and n² − 17·4^k (k ≤ 25) is a sum of two squares; for
n − 2 and n + 2, 17 and 25 of the 81 are. Here the search was complete: every odd n below
239694783 was decided exactly and has a representation, and by the first remark above the
smallest counterexample is odd. So **239694783 is the smallest counterexample to Conjecture 2.**

**Conjecture 3 of A300448** uses m ∈ {1, 5, 65}. It fails at
n = 68136270467 = 103 · 661517189: none of the 104 numbers n² − 4^k (k ≤ 35), n² − 5·4^k
(k ≤ 34) and n² − 65·4^k (k ≤ 32) is a sum of two squares, while 21 and 26 of the 104 are for
n − 2 and n + 2. It fails again at the prime n = 80570407837 (107 numbers; 19 and 20 of them
are sums of two squares for n − 2 and n + 2). The search sieved every odd n from 2·10^7 to
10^11 and decided exactly those with at most 14 pairs (m, k) left after the sieve; these two
are the only failures among them. A second run, on the odd n below 68136270467 only, raised
that limit to 20 and found no smaller failure. The other odd n were not decided, so
**68136270467 is not known to be the smallest counterexample to Conjecture 3.** Sun remarks
that Conjecture 3 would imply the conjecture of [A300441](https://oeis.org/A300441), which
allows infinitely many values of m. That conjecture is not settled by these numbers.

### The "2-3 conjecture" (A301452)

On 21 March 2018 Sun conjectured that for every integer n > 1 there is a k ≥ 0 such that
n² − 2·4^k or n² − 3·4^k is of the form x² + 2y² (OEIS A301452, where he calls this the "2-3
conjecture" and reports a check for 2 ≤ n ≤ 5·10^7). It is Question 2 of his MathOverflow
question cited above; Question 1 there is the 2-5 conjecture.

**For n = 386749710995 = 5 · 77349942199 there is no such k.** In the notation of the OEIS
entry, a(386749710995) = 0.

* There are 76 numbers to consider: n² − 2·4^k and n² − 3·4^k for k = 0, …, 37.
* A positive integer is of the form x² + 2y² exactly when every prime ≡ 5 or 7 (mod 8) divides
  it to an even power. Each of the 76 numbers has a prime ≡ 5 or 7 (mod 8) to an odd power.
* For comparison, 17 of the 76 numbers are of the form x² + 2y² for n − 2, and 16 for n + 2.

Four more counterexamples came out of the same search: 985224034231 = 73 · 37181 · 362987,
1241790594939 = 3³ · 181 · 2663 · 95419, 1572550146081 = 3³ · 41057 · 1418579 and
9224999068015 = 5 · 1487 · 1240753069.

How these numbers relate to smaller ones:

* As for the 2-5 conjecture, the smallest counterexample is odd.
* Every odd n from 5·10^7 to 2·10^10 was sieved, and those with at most 16 pairs (m, k) left
  were decided exactly; from 2·10^10 to 10^11 the limit was 12. No failure turned up.
* From 10^11 to 10^12 the search looked only at the n in 210592 of the 2400671 residue classes
  modulo 7 · 13 · 23 · 31 · 37 (the classes in which those five primes strike the most pairs),
  with the limit 13. That gave 386749710995 and 985224034231. A run of the same kind from
  10^12 to 10^13, with the limit 14, was stopped after about 44 % of its 280672 classes; it had
  found the other three.
* So most n below 386749710995 were never decided, and **386749710995 is not known to be the
  smallest counterexample.**

Sun's Square Conjecture ([A301471](https://oeis.org/A301471): some n² − 3·2^k is of the form
x² + 2y²) is false as well, as Sun found himself in April 2018; the counterexample is
5884015571.

## Verify it yourself

```
python verify.py --quick     # about half a minute
python verify.py             # adds slow brute-force cross-checks, a few minutes
```

Every check enumerates all admissible choices of the sparse terms (the powers of 2, 3, 4, 5
and 9) and shows that what is left over never has the required shape. That is a finite
computation, so a passing check is a proof.

* **x² + 5y² and the polygonal forms.** Plain search: loop over one variable and test the
  rest with an exact integer square root.
* **x² + 2y².** Two independent tests: the same plain search, and the classical criterion
  that m = x² + 2y² exactly when every prime ≡ 5 or 7 (mod 8) divides m to an even power.
  The two agree on all 983 leftovers for n = 4375339279.
* **Primes (A303932).** The files in [`certificates/`](certificates/) list every pair (k, m)
  for each counterexample N. For most pairs they give a proper divisor of p = N − 2^k − 3^m.
  For the others they give r with r² ≡ 11 (mod p), which shows 11 is a quadratic residue
  (Conjecture 2), or p satisfies 11^((p−1)/2) ≡ −1 (mod p), which by Euler's criterion shows
  11 is a non-residue if p is prime (Conjecture 1). `check_cert_a303932.py` checks this with
  integer arithmetic alone and trusts no primality test.
* **Primes (A157242).** [`certificates/A157242_N107576474238337.csv`](certificates/A157242_N107576474238337.csv)
  lists all 1974 pairs (x, y) with p = N − 2^x − 11·2^y, p mod 6 and, for the 504 values in
  the class 5 mod 6, a proper divisor of p. `check_cert_a157242.py` rebuilds the list of
  pairs itself and checks every row with integer arithmetic alone. `verify.py` also recounts
  with a deterministic Miller-Rabin test (0 representations for N, 51 and 43 for N − 2 and
  N + 2) and checks its counting function against the 60 terms of A157242 printed in the OEIS.
* **Sums of two squares (A300510, A300448) and numbers x² + 2y² (A301452).**
  [`certificates/A300510_n80397893435.csv`](certificates/A300510_n80397893435.csv) lists the 72
  numbers v = n² − m·4^k and, for each, a divisor d ≡ 3 (mod 4) of v that is coprime to v/d.
  Such a d has a prime factor ≡ 3 (mod 4) of odd multiplicity, and then so has v, so no
  primality test is needed. `check_cert_a300510.py` rebuilds the list of pairs (m, k) itself and
  checks every row.
  [`certificates/A300448_conj1_n80397893435.csv`](certificates/A300448_conj1_n80397893435.csv)
  does the same for the 109 numbers with m ∈ {1, 2, 5}, and
  [`certificates/A300448_conj2_n239694783.csv`](certificates/A300448_conj2_n239694783.csv) for
  the 81 numbers of Conjecture 2. The two files `certificates/A300448_conj3_n*.csv` cover
  Conjecture 3. The five files `certificates/A301452_n*.csv` do the same for the 2-3
  conjecture, with a divisor d ≡ 5 or 7 (mod 8) coprime to v/d; such a d has a prime factor
  ≡ 5 or 7 (mod 8) of odd multiplicity. `verify.py` also factors every number itself and
  recounts (0 of 72 and 0 of 109 for 80397893435, 19 and 18 of 72 for its neighbours; 0 of 81
  for 239694783, 17 and 25 for its neighbours; 0 of 104 for 68136270467 and 0 of 76 for
  386749710995, with their neighbours as controls), checks its counting functions against the
  80 terms of A300510, of A300448 and of A301452 printed in the OEIS, and re-checks the
  published failures of related statements: 51883659, 56173837, 65525725 and 42211965 for the
  variants with 29, 65, 113 and 162 in place of 5, 13018247 for {1, 5, 11}, and Sun's
  5884015571 for his Square Conjecture. For 9224999068015 it checks the certificate only: the
  squares involved are beyond the range in which its primality test is proven.
* **Lean.** [`lean/SunA157242.lean`](lean/SunA157242.lean) proves in Lean 4, using only the
  core library, that 107576474238337 has no representation, and hence that the conjecture of
  A157242 is false. [`lean/SunA303932.lean`](lean/SunA303932.lean) does the same for A303932,
  Conjecture 2, at 2n = 68211369278. The repository's GitHub Actions workflow runs the Lean
  kernel on both files (Lean 4.36.0 at the time of writing). Each final theorem depends on the
  axioms `propext` and `Quot.sound` only, with no `sorry` and no `native_decide`.

## How the numbers were found

* `search/mixsum2.c` is an exhaustive search for n ≠ m + a + b with a, b from two short lists
  and m from a base set (primes in given residue classes, or numbers x² + 2y²). It builds the
  base set as a bitset, ORs shifted copies of it, and then tests every surviving n exactly
  against all pairs (a, b). It scans about 10^9 integers per second on an 8-core desktop.
* `search/direct_scan.c` is a second, independent program: it marks every sum of two
  quadratic values in a bit table by a double loop and then checks each n by table look-up.
* `search/sunmix_multi.c` is the exhaustive search for odd N ≠ p + 2^x + k·2^y with p in a
  residue class modulo 6. It works like `mixsum2`, on a bitset of the odd numbers.
* `search/suntarget2.c` is a targeted search for the same problem at sizes where testing every
  number is out of reach. For N in certain residue classes modulo 9699690 = 2·3·5·7·11·13·17·19
  most candidates p are forced to be divisible by a prime up to 19 or to lie in the wrong class
  modulo 6. The program ranks the classes by the number of candidates that survive, scans the
  best ones (about 3·10^11 numbers per octave [X, 2X) in the runs here), counts for each N
  how many candidates also survive the primes from 23 to 1021, and tests the most promising N
  exactly, over all pairs (x, y). It reports only numbers that pass that complete test, so a
  reported number is a true exception, but it does not look at most numbers. As a check, for
  k = 7 and p ≡ 5 (mod 6) it finds all four exceptions between 4·10^10 and 3.2·10^12 while
  looking at 1.3·10^11 of the 1.6·10^12 odd numbers in that range.
* `search/sqsparse.c` searches for odd n such that no n² − m·4^k is a sum of two squares (or
  is of the form x² + 2y²). A prime q ≡ 3 (mod 4) divides n² − m·4^k only if n ≡ ±s·2^k (mod q)
  with s² ≡ m, so a sieve over n with those primes up to 65536 gives, for every n, the pairs
  (m, k) that none of them divides. The n with few pairs left are then decided exactly, by
  factoring, and a number is reported only if every one of its pairs fails that exact test. Run
  on all odd n from 6·10^9 to 10^11 with the limit "at most 12 pairs left", it reported
  80397893435 and nothing else. As a check, when it decides every odd n below 7·10^7 it finds
  exactly one failure for each variant with 29, 65, 113 or 162 in place of 5, namely the number
  recorded in the OEIS entry, and with the limit 14 it finds exactly Sun's 5884015571 among
  the odd n below 6·10^9 for the Square Conjecture. For m ∈ {1, 5, 17} a run on all odd n from
  2·10^7 to 10^11 with the limit 14 reported 239694783 and nothing else, and a second run that
  decided every odd n below 239694785 exactly found no other failure. The same kind of run for
  m ∈ {1, 5, 65} reported 68136270467 and 80570407837. For the 2-3 conjecture the program was
  also used in a second mode, which sieves only the residue classes of n modulo a product of
  small bad primes in which those primes already strike the most pairs; the ranges are listed
  in the section on that conjecture.

| Statement | Range searched exhaustively | Counterexamples in that range |
|---|---|---|
| A157242 (k = 11, p ≡ 5 mod 6) | odd N < 1.0758·10^14 (`sunmix_multi`) | 107576474238337, and none other beyond Sun's 33, 39, 45, 171, 247 |
| A157237 (k = 11, p ≡ 1 mod 6) | odd N < 10^13 (`sunmix_multi`) | none beyond Sun's 35, 41, 47, 101, 167, 2021, 119171 |
| A300448, Conjecture 2 | odd n < 239694785, each decided exactly (`sqsparse`); an even counterexample would give a smaller odd one | 239694783 |
| A302984 | n < 10^11 (`mixsum2`); n ≤ 10^10 again with `direct_scan` | 18, listed in [`data/`](data/), all ≡ 79 (mod 80); both programs give the same four below 10^10 |
| A302982 | n ≤ 10^9 (`direct_scan`) | 558556427, 640123571, 683305451 |
| A303932, Conjecture 1 | even numbers below 10^12 (`mixsum2`) | 257558146370, 749360555468, 799324864730 |
| A303932, Conjecture 2 | even numbers below 10^12 (`mixsum2`) | 68211369278, 273299806988 |
| A308644 Conj. 2 (r = 3) | n ≤ 3·10^8 (`direct_scan`) | 7654105 |
| A308644 Conj. 1 (second assertion) | n ≤ 3·10^8 (`direct_scan`) | 108011768 |
| A308656 Conj. 2 (f = x(4x+1)) | n ≤ 3·10^8 (`direct_scan`) | 235644514 |
| A308644 Conj. 1, first assertion (d(7d+1)/2) | n ≤ 3·10^8 (`direct_scan`) | none |

For A157242, for A300448, Conjecture 2, and for the seven statements of 2018 and 2019 the
counterexample in the first table is the smallest one, and that rests on these runs. For
A302984 two independent programs agree; for the other statements one program made the
exhaustive run.

The same exhaustive run below 10^13 covered k = 7, Sun's earlier guesses of the same shape,
which D. S. McNeil refuted in 2009 (see Sun's post of 25 February 2009). For p ≡ 1 (mod 6) the
only exception above 33 is McNeil's 445466268731. For p ≡ 5 (mod 6) the exceptions above 1983
are McNeil's 1433987797, 45483661477, 85305551767, 367084269651, and one more below 10^13:
1572650032347.

The targeted search for k = 11 covered the octaves from 10^13 to 4.5·10^18 for both classes.
Its only hit was 107576474238337.

In A303932, Conjecture 1, the search also admitted p = 11 (where 11 ≡ 0) as if it were a
residue. The three numbers found have no representation either way, and no further even
number below 10^12 fails when p = 11 is excluded.

## Prior work

We looked for earlier reports of these counterexamples on the OEIS pages themselves, on
GitHub (repositories and code), on arXiv and by web search of the numbers, on
26 November 2026, and found none. Mathematical Reviews, zbMATH and Chinese-language sources
were not searched. If you know of an earlier source, please open an issue.

For A300510, A300448 and A301452 we read the OEIS entries, searched the OEIS for the numbers
given here, read Sun's MathOverflow question of April 2020 with its answer and comments, and
searched GitHub and the web. No counterexample to any of these conjectures is recorded in any
of them. None of the three sequences is among the OEIS conjectures in Google DeepMind's
formal-conjectures.

For A157242 we also read the OEIS edit history and draft queue of the entry, Sun's posts to
the Number Theory List of February and March 2009, and his web page
[Mixed Sums of Primes and Other Terms](http://maths.nju.edu.cn/~zwsun/MSPT.htm) (last
modified in November 2010). The largest verification bound we found is McNeil's 10^13, and no
counterexample is mentioned anywhere. Sun's
[Notes on conjectures](http://maths.nju.edu.cn/~zwsun/conj.pdf) (last update 30 August 2009)
list his mixed-sum conjectures with Fibonacci, Pell and Catalan numbers but not this one. The
conjecture is not among the OEIS conjectures
formalized in Google DeepMind's
[formal-conjectures](https://github.com/google-deepmind/formal-conjectures); its companion
A157237 is, and is marked open there. The k = 7 analogue with p ≡ 5 (mod 6), A157225, is
recorded there as refuted at 1433987797, a number that Sun's post of 25 February 2009 already
credits to McNeil.

Close relatives of these conjectures were already known to fail:

* A303338 (x² + 2y² + 3·2^z + 4^w) at 12558941213 and A303637 (x² + y² + 2^z + 5·2^w) at
  18836421387, found by Jiao-Min Lin in 2022 and recorded in the OEIS by Sun.
* A308656, Conjecture 1, at 2109982225 (Lin, 2022).
* A308640, Conjecture 1, at 1764630 (Scott Moore, 2026).
* A301534 (every prime p > 7 with p ≡ 7 mod 12 is x² + 3y² + 15·2^z) at p = 28662666283,
  reported on 26 November 2026 in the GitHub repository `senamakel/conjectures-lean`. We
  re-checked that number.

## Provenance

The searches and checks were run on 25 and 26 November 2026 on one desktop PC by the owner
of this repository, working with Claude (Anthropic's AI model) in Claude Code. The
counterexamples to A302982, A308644 and A308656 were first found by a quick scan and then
re-derived with the programs here. The counterexample to A157242 was found by
`search/suntarget2.c` and then confirmed by the exact test of `search/sunmix_multi.c`, by
`verify.py`, and by factoring all 504 candidates with the `factor` program of uutils
coreutils 0.10.0. A second AI system (OpenAI's GPT-6.2) was asked to referee the claim; it
recomputed the counts quoted above for 107576474238337 with a program of its own and obtained
the same numbers. The counterexample to A300510 and A300448 was found by `search/sqsparse.c`
and confirmed by `make_cert_a300510.py`, which shares no code with it, and by `verify.py`.
GPT-6.2 refereed that claim too: with a program of its own it factored all 109 numbers for n
and for n ± 2, obtained the same counts, and checked every row of the first certificate; it
also searched the edit histories and drafts of the two OEIS entries, the versions of Sun's
paper arXiv:1701.05868 and the papers citing it, and found the MathOverflow thread cited
above, but no earlier counterexample. The counterexamples to Conjectures 2 and 3 of A300448
and to A301452 were found and confirmed in the same way; they were not part of that referee
report. Every counterexample to the thirteen statements that is named in this file, and every
exception listed for k = 7 and k = 11, is covered by `verify.py`.

## Files

| Path | Contents |
|---|---|
| `verify.py` | re-checks every counterexample |
| `certificates/` | primality-free certificates for the A157242 counterexample, the A300510, A300448 and A301452 counterexamples and the five A303932 counterexamples |
| `make_cert_a157242.py`, `check_cert_a157242.py` | writer and checker for the A157242 certificate |
| `make_cert_a300510.py`, `check_cert_a300510.py` | writer and checker for the A300510, A300448 and A301452 certificates |
| `make_cert_a303932.py`, `check_cert_a303932.py` | writer and checker for the A303932 certificates |
| `lean/` | Lean 4 proofs for A157242 and for A303932, Conjecture 2, with their generators |
| `data/` | the complete lists found for A302984, A302982 and A303932 |
| `search/` | the search programs and the definitions of the term lists |
| `OEIS_COMMENTS.md` | suggested comments for the OEIS entries |

## License

MIT, see [`LICENSE`](LICENSE).
