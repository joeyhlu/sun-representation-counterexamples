# Counterexamples to eight conjectures of Zhi-Wei Sun on representations of integers

Zhi-Wei Sun has posted in the OEIS many conjectures of the form "every positive integer is a
sum of a few sparse terms". Eight of those statements are false: one from 2009 on
p + 2^x + 11·2^y, and seven from 2018 and 2019. Each counterexample is one integer, and
`python verify.py --quick` re-checks all of them in about 20 seconds using only the Python
standard library.

| OEIS | Conjecture | Checked before up to | Counterexample |
|---|---|---|---|
| [A157242](https://oeis.org/A157242) | every odd N > 28 other than 33, 39, 45, 171, 247 is p + 2^x + 11·2^y with x, y ≥ 1 and p a prime ≡ 5 (mod 6) | 10^13 (D. S. McNeil, 2009, as reported by Sun) | **107576474238337** |
| [A302984](https://oeis.org/A302984) | every n > 5 is x² + 2y² + 2^z + 5·2^w | 10^9 | **4375339279** |
| [A302982](https://oeis.org/A302982) | every n > 3 is x² + 5y² + 2^z + 3·2^w | 2·10^8 | **558556427** |
| [A303932](https://oeis.org/A303932), Conjecture 1 | every even number 2n > 2 is p + 2^k + 3^m with p a prime such that 11 is a quadratic residue mod p | n ≤ 5·10^8 | **2n = 257558146370** |
| [A303932](https://oeis.org/A303932), Conjecture 2 | every even number 2n > 4 is p + 2^k + 3^m with p a prime such that 11 is a quadratic non-residue mod p | n ≤ 5·10^8 | **2n = 68211369278** |
| [A308644](https://oeis.org/A308644), Conjecture 2 (r = 3) | every n > 0 is (3^a·4^b)² + c(3c+1)/2 + d(7d+3)/2 | 10^6 | **7654105** |
| [A308644](https://oeis.org/A308644), Conjecture 1 (second assertion) | every n > 0 is (3^a·5^b)² + c(3c+1)/2 + d(7d+3)/2 | 10^6 | **108011768** |
| [A308656](https://oeis.org/A308656), Conjecture 2 (f(x) = x(4x+1)) | every n > 0 is (2^a·9^b)² + c(4c+1) + d(3d+1)/2 | 10^6 | **235644514** |

In the rows after the first, x, y, z, w, a, b, k, m are non-negative integers and c, d are
arbitrary integers, as in the OEIS entries. Those seven counterexamples are the smallest ones;
for A157242 see the next section.

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
* A test of every odd number from 10^13 up to N was still running when this was written, so
  **107576474238337 is not yet known to be the smallest counterexample.**

The companion conjecture with p ≡ 1 (mod 6) ([A157237](https://oeis.org/A157237): every odd
N > 30 other than 35, 41, 47, 101, 167, 2021, 119171) is still open. Below 10^13 its only
exceptions are those seven, and the targeted search, run up to 4.5·10^18, found no
counterexample to it.

## Verify it yourself

```
python verify.py --quick     # about 20 seconds
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

| Statement | Range searched exhaustively | Counterexamples in that range |
|---|---|---|
| A157242 (k = 11, p ≡ 5 mod 6) | odd N < 10^13 (`sunmix_multi`) | none beyond Sun's 33, 39, 45, 171, 247 |
| A157237 (k = 11, p ≡ 1 mod 6) | odd N < 10^13 (`sunmix_multi`) | none beyond Sun's 35, 41, 47, 101, 167, 2021, 119171 |
| A302984 | n < 10^11 (`mixsum2`); n ≤ 10^10 again with `direct_scan` | 18, listed in [`data/`](data/), all ≡ 79 (mod 80); both programs give the same four below 10^10 |
| A302982 | n ≤ 10^9 (`direct_scan`) | 558556427, 640123571, 683305451 |
| A303932, Conjecture 1 | even numbers below 10^12 (`mixsum2`) | 257558146370, 749360555468, 799324864730 |
| A303932, Conjecture 2 | even numbers below 10^12 (`mixsum2`) | 68211369278, 273299806988 |
| A308644 Conj. 2 (r = 3) | n ≤ 3·10^8 (`direct_scan`) | 7654105 |
| A308644 Conj. 1 (second assertion) | n ≤ 3·10^8 (`direct_scan`) | 108011768 |
| A308656 Conj. 2 (f = x(4x+1)) | n ≤ 3·10^8 (`direct_scan`) | 235644514 |
| A308644 Conj. 1, first assertion (d(7d+1)/2) | n ≤ 3·10^8 (`direct_scan`) | none |

For the seven statements of 2018 and 2019 the counterexample in the first table is the
smallest one, and that rests on these runs. For A302984 two independent programs agree; for
the other statements one program made the exhaustive run.

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

## Provenance

The searches and checks were run on 25 and 26 November 2026 on one desktop PC by the owner
of this repository, working with Claude (Anthropic's AI model) in Claude Code. The
counterexamples to A302982, A308644 and A308656 were first found by a quick scan and then
re-derived with the programs here. The counterexample to A157242 was found by
`search/suntarget2.c` and then confirmed by the exact test of `search/sunmix_multi.c`, by
`verify.py`, and by factoring all 504 candidates with the `factor` program of uutils
coreutils 0.10.0. A second AI system (OpenAI's GPT-6.2) was asked to referee the claim; it
recomputed the counts quoted above for 107576474238337 with a program of its own and obtained
the same numbers. Every counterexample to the eight statements that is named in this file,
and every exception listed for k = 7 and k = 11, is covered by `verify.py`.

## Files

| Path | Contents |
|---|---|
| `verify.py` | re-checks every counterexample |
| `certificates/` | primality-free certificates for the A157242 counterexample and the five A303932 counterexamples |
| `make_cert_a157242.py`, `check_cert_a157242.py` | writer and checker for the A157242 certificate |
| `make_cert_a303932.py`, `check_cert_a303932.py` | writer and checker for the A303932 certificates |
| `lean/` | Lean 4 proofs for A157242 and for A303932, Conjecture 2, with their generators |
| `data/` | the complete lists found for A302984, A302982 and A303932 |
| `search/` | the search programs and the definitions of the term lists |
| `OEIS_COMMENTS.md` | suggested comments for the OEIS entries |

## License

MIT, see [`LICENSE`](LICENSE).
