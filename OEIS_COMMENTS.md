# Suggested OEIS comments

Draft comments for the COMMENTS sections of the entries concerned. The OEIS appends the
contributor's signature and date (`~~~~`).

## A157242

```
The conjecture is false: a(53788237119169) = 0, that is, the odd number 107576474238337 = 13*43*2129*7411*12197 is not of the form p + 2^x + 11*2^y with p a prime congruent to 5 mod 6 and x, y positive integers. Among the 1974 pairs (x,y) with 2^x + 11*2^y < 107576474238337, p = 107576474238337 - 2^x - 11*2^y is congruent to 5 mod 6 exactly for the 504 pairs with x even and y odd, and all 504 values are composite. This is the least counterexample: a test of every odd number below 107576474238337 found no other exception beyond 33, 39, 45, 171, 247 (up to 10^13 this agrees with D. S. McNeil's verification reported in Sun's message to the Number Theory List of Feb 28 2009). - ~~~~
```

Suggested LINKS line:

```
Zhi-Wei Sun, <a href="https://listserv.nodak.edu/cgi-bin/wa.exe?A2=NMBRTHRY;20d0ba9.0902">A curious conjecture about p+2^x+11*2^y</a>, Number Theory List, Feb 25 2009.
```

The test of all odd numbers between 10^13 and 107576474238337 finished on 26 November 2026
and found no other exception (see the README), so the comment says "least".

## A300510

```
The conjecture is false: a(80397893435) = 0. For n = 80397893435 = 5*1163*1237*11177 and each k = 0..35, both n^2 - 2*4^k and n^2 - 5*4^k have a prime factor congruent to 3 mod 4 with odd exponent, so neither is a sum of two squares. This n was found by a search that is not exhaustive between 6*10^9 and 10^11, so it need not be the least counterexample. - ~~~~
```

Suggested LINKS line:

```
Zhi-Wei Sun, <a href="https://mathoverflow.net/questions/356905">Write n^2 as x^2+y^2+2*4^z or x^2+y^2+5*4^z</a>, MathOverflow, April 2020.
```

## A300448

```
Conjecture 1 is false: a(80397893435) = 0. For n = 80397893435 and each admissible k, the numbers n^2 - 4^k (k = 0..36), n^2 - 2*4^k and n^2 - 5*4^k (k = 0..35) all have a prime factor congruent to 3 mod 4 with odd exponent. - ~~~~
```

```
Conjecture 2 is false: the square of 239694783 = 3*47*179*9497 cannot be written as 4^k*m + x^2 + y^2 with k, x, y nonnegative integers and m among 1, 5, 17. This is the least positive integer with that property. - ~~~~
```

```
Conjecture 3 is false: the squares of 68136270467 = 103*661517189 and of the prime 80570407837 cannot be written as 4^k*m + x^2 + y^2 with k, x, y nonnegative integers and m among 1, 5, 65. These numbers were found by a search that is not exhaustive, so 68136270467 need not be the least such number. - ~~~~
```

## A301452

```
The 2-3 conjecture is false: a(386749710995) = 0. For n = 386749710995 = 5*77349942199 and each k = 0..37, both n^2 - 2*4^k and n^2 - 3*4^k have a prime factor congruent to 5 or 7 mod 8 with odd exponent, so neither is of the form x^2 + 2*y^2. Also a(n) = 0 for n = 985224034231, 1241790594939, 1572550146081, 9224999068015. These n were found by a search that is not exhaustive, so 386749710995 need not be the least counterexample. - ~~~~
```

Suggested LINKS line:

```
Zhi-Wei Sun, <a href="https://mathoverflow.net/questions/356905">Write n^2 as x^2+y^2+2*4^z or x^2+y^2+5*4^z</a>, MathOverflow, April 2020.
```

## A157237

```
There is no counterexample to the conjecture below 10^13. - ~~~~
```

## A157218 and A157225 (k = 7), for the record

Neither entry mentions that the statements Sun first considered for k = 7 fail. Sun's message
of Feb 25 2009 credits the counterexamples to D. S. McNeil.

```
A157218: The odd number 445466268731 is not of the form p + 2^x + 7*2^y with p a prime congruent to 1 mod 6 and x, y positive integers (D. S. McNeil, 2009, see Sun's message to the Number Theory List of Feb 25 2009); it is the only odd number in (33, 10^13) with this property. - ~~~~
```

```
A157225: The odd numbers greater than 1983 and less than 10^13 that are not of the form p + 2^x + 7*2^y with p a prime congruent to 5 mod 6 and x, y positive integers are 1433987797, 45483661477, 85305551767, 367084269651 (D. S. McNeil, 2009, see Sun's message to the Number Theory List of Feb 25 2009) and 1572650032347. - ~~~~
```

## A302984

```
The conjecture is false: a(4375339279) = 0, and 4375339279 is the least n > 5 with a(n) = 0. The numbers n with 5 < n < 10^11 and a(n) = 0 are 4375339279, 6285378879, 8525602079, 9731038079, 13523100479, 24920030079, 29312698479, 30034054079, 32269674879, 40072842719, 51729946479, 58532818879, 67179857279, 71235478079, 72161342879, 82358435279, 85090270079, 91924062079; all of them are congruent to 79 modulo 80. - ~~~~
```

## A302982

```
The conjecture is false: a(558556427) = 0. The only n with 3 < n <= 10^9 and a(n) = 0 are 558556427, 640123571 and 683305451. - ~~~~
```

## A303932

```
Both conjectures are false. Conjecture 1 fails for n = 128779073185: among the numbers 257558146370 - 2^k - 3^m (k, m >= 0) exactly 19 are prime, and 11 is a quadratic nonresidue modulo each of them. Conjecture 2 fails for n = 34105684639: among the numbers 68211369278 - 2^k - 3^m exactly 20 are prime, and 11 is a quadratic residue modulo each of them. For 2n < 10^12 the counterexamples to Conjecture 1 are 2n = 257558146370, 749360555468, 799324864730, and those to Conjecture 2 are 2n = 68211369278, 273299806988. - ~~~~
```

## A308644

```
Conjecture 2 is false for r = 3: 7654105 cannot be written as (3^a*4^b)^2 + c*(3c+1)/2 + d*(7d+3)/2 with a, b nonnegative integers and c, d integers. The second assertion of Conjecture 1 is also false: 108011768 cannot be written as (3^a*5^b)^2 + c*(3c+1)/2 + d*(7d+3)/2. In both cases this is the only counterexample up to 3*10^8. The first assertion of Conjecture 1, a(n) > 0, holds for all n = 1..3*10^8. - ~~~~
```

## A308656

```
Conjecture 2 is false for f(x) = x*(4x+1): 235644514 cannot be written as (2^a*9^b)^2 + c*(4c+1) + d*(3d+1)/2 with a, b nonnegative integers and c, d integers; this is the only counterexample up to 3*10^8. - ~~~~
```
