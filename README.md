# Counterexamples to seven conjectures of Zhi-Wei Sun on representations of integers

In 2018 and 2019 Zhi-Wei Sun posted in the OEIS a family of conjectures of the form
"every positive integer is a sum of a few sparse terms". Seven of those statements are false.
Each counterexample is one integer, and `python verify.py --quick` re-checks all of them in
about 20 seconds using only the Python standard library.

| OEIS | Conjecture | Checked by its author up to | Smallest counterexample |
|---|---|---|---|
| [A302984](https://oeis.org/A302984) | every n > 5 is x² + 2y² + 2^z + 5·2^w | 10^9 | **4375339279** |
| [A302982](https://oeis.org/A302982) | every n > 3 is x² + 5y² + 2^z + 3·2^w | 2·10^8 | **558556427** |
| [A303932](https://oeis.org/A303932), Conjecture 1 | every even number 2n > 2 is p + 2^k + 3^m with p a prime such that 11 is a quadratic residue mod p | n ≤ 5·10^8 | **2n = 257558146370** |
| [A303932](https://oeis.org/A303932), Conjecture 2 | every even number 2n > 4 is p + 2^k + 3^m with p a prime such that 11 is a quadratic non-residue mod p | n ≤ 5·10^8 | **2n = 68211369278** |
| [A308644](https://oeis.org/A308644), Conjecture 2 (r = 3) | every n > 0 is (3^a·4^b)² + c(3c+1)/2 + d(7d+3)/2 | 10^6 | **7654105** |
| [A308644](https://oeis.org/A308644), Conjecture 1 (second assertion) | every n > 0 is (3^a·5^b)² + c(3c+1)/2 + d(7d+3)/2 | 10^6 | **108011768** |
| [A308656](https://oeis.org/A308656), Conjecture 2 (f(x) = x(4x+1)) | every n > 0 is (2^a·9^b)² + c(4c+1) + d(3d+1)/2 | 10^6 | **235644514** |

Here x, y, z, w, a, b, k, m are non-negative integers and c, d are arbitrary integers, as in
the OEIS entries.

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
* **Lean.** [`lean/SunA303932.lean`](lean/SunA303932.lean) proves the failure of A303932,
  Conjecture 2, at 2n = 68211369278 in Lean 4, using only the core library. The repository's
  GitHub Actions workflow runs the Lean kernel on it (Lean 4.36.0 at the time of writing);
  the theorem depends on the axioms `propext` and `Quot.sound` only, with no `sorry` and no
  `native_decide`.

## How the numbers were found

* `search/mixsum2.c` is an exhaustive search for n ≠ m + a + b with a, b from two short lists
  and m from a base set (primes in given residue classes, or numbers x² + 2y²). It builds the
  base set as a bitset, ORs shifted copies of it, and then tests every surviving n exactly
  against all pairs (a, b). It scans about 10^9 integers per second on an 8-core desktop.
* `search/direct_scan.c` is a second, independent program: it marks every sum of two
  quadratic values in a bit table by a double loop and then checks each n by table look-up.

| Statement | Range searched exhaustively | Counterexamples in that range |
|---|---|---|
| A302984 | n < 10^11 (`mixsum2`); n ≤ 10^10 again with `direct_scan` | 18, listed in [`data/`](data/), all ≡ 79 (mod 80); both programs give the same four below 10^10 |
| A302982 | n ≤ 10^9 (`direct_scan`) | 558556427, 640123571, 683305451 |
| A303932, Conjecture 1 | even numbers below 10^12 (`mixsum2`) | 257558146370, 749360555468, 799324864730 |
| A303932, Conjecture 2 | even numbers below 10^12 (`mixsum2`) | 68211369278, 273299806988 |
| A308644 Conj. 2 (r = 3) | n ≤ 3·10^8 (`direct_scan`) | 7654105 |
| A308644 Conj. 1 (second assertion) | n ≤ 3·10^8 (`direct_scan`) | 108011768 |
| A308656 Conj. 2 (f = x(4x+1)) | n ≤ 3·10^8 (`direct_scan`) | 235644514 |
| A308644 Conj. 1, first assertion (d(7d+1)/2) | n ≤ 3·10^8 (`direct_scan`) | none |

The "smallest counterexample" column of the first table rests on these runs. For A302984 two
independent programs agree; for the other statements one program made the exhaustive run.

In A303932, Conjecture 1, the search also admitted p = 11 (where 11 ≡ 0) as if it were a
residue. The three numbers found have no representation either way, and no further even
number below 10^12 fails when p = 11 is excluded.

## Prior work

We looked for earlier reports of these counterexamples on the OEIS pages themselves, on
GitHub (repositories and code), on arXiv and by web search of the numbers, on
26 November 2026, and found none. Mathematical Reviews, zbMATH and Chinese-language sources
were not searched. If you know of an earlier source, please open an issue.

Close relatives of these conjectures were already known to fail:

* A303338 (x² + 2y² + 3·2^z + 4^w) at 12558941213 and A303637 (x² + y² + 2^z + 5·2^w) at
  18836421387, found by Jiao-Min Lin in 2022 and recorded in the OEIS by Sun.
* A308656, Conjecture 1, at 2109982225 (Lin, 2022).
* A308640, Conjecture 1, at 1764630 (Scott Moore, 2026).

## Provenance

The searches and checks were run on 25 and 26 November 2026 on one desktop PC by the owner
of this repository, working with Claude (Anthropic's AI model) in Claude Code. The
counterexamples to A302982, A308644 and A308656 were first found by a quick scan and then
re-derived with the programs here. Every counterexample named in this file is covered by
`verify.py`.

## Files

| Path | Contents |
|---|---|
| `verify.py` | re-checks every counterexample |
| `certificates/` | primality-free certificates for the five A303932 counterexamples |
| `make_cert_a303932.py`, `check_cert_a303932.py` | writer and checker for those certificates |
| `lean/` | Lean 4 proof for A303932, Conjecture 2, and its generator |
| `data/` | the complete lists found for A302984, A302982 and A303932 |
| `search/` | the search programs and the definitions of the term lists |
| `OEIS_COMMENTS.md` | suggested comments for the OEIS entries |

## License

MIT, see [`LICENSE`](LICENSE).
