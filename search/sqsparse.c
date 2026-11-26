/* sqsparse.c - search for odd n such that NO number n^2 - m*4^k (m in a given list, k >= 0,
 * m*4^k < n^2) is represented by a fixed binary quadratic form.
 *
 *   form 1:  x^2 + y^2     v > 0 is represented iff every prime q = 3 (mod 4) divides v to an even power
 *   form 2:  x^2 + 2 y^2   v > 0 is represented iff every prime q = 5, 7 (mod 8) divides v to an even power
 *
 * Examples (Zhi-Wei Sun, OEIS):
 *   form 1, m = 2,5   "2-5 conjecture" (A300510): for every n > 1 some n^2 - 2*4^k or n^2 - 5*4^k is x^2 + y^2
 *   form 2, m = 2,3   "2-3 conjecture" (A301452): for every n > 1 some n^2 - 2*4^k or n^2 - 3*4^k is x^2 + 2y^2
 *   form 2, m = 3,6   Square Conjecture (A301471): some n^2 - 3*2^k is x^2 + 2y^2   (3*2^k = 3*4^j or 6*4^j)
 *
 * Only odd n are scanned: n^2 - m*4^(k+1) = 4 ((n/2)^2 - m*4^k), so if an even n has no representation
 * then neither has n/2, and the least n without a representation is odd.
 *
 * Method.  Call a prime "bad" if it must divide a represented number to an even power.  For a bad prime
 * q not dividing m, q | n^2 - m*4^k iff n = +-s*2^k (mod q) with s^2 = m (mod q).  Stage 1 sieves blocks
 * of odd n with the bad primes up to B and records, for every n, the set of pairs (m, k) not hit
 * ("alive").  Stage 2 takes the n with at most T alive pairs and decides alive pairs exactly (trial
 * division, strong probable prime tests, Pollard rho) until one is represented.  If no alive pair is
 * represented, ALL pairs of that n are decided exactly, because stage 1 also strikes a pair when
 * q^2 | n^2 - m*4^k.  An n is reported only after that exact check.  With T below the number of pairs
 * the search is not exhaustive: an n with more than T alive pairs is skipped.
 *
 * Usage:
 *   ./sqsparse scan  FORM M1,M2,.. NLO NHI T [threads] [B]       every odd n in [NLO, NHI) is sieved
 *   ./sqsparse stats FORM M1,M2,.. NLO NHI T [threads] [B]       table: alive pairs against represented pairs
 *   ./sqsparse cls   FORM M1,M2,.. NLO NHI T Q1,Q2,.. A0 [threads] [B]
 *         only the residue classes of n modulo Q1*Q2*... in which at most A0 pairs survive the primes
 *         Q1, Q2, ... are sieved.  A0 = 0: print the distribution of that number and stop.
 *         A0 = -p: keep the classes with the fewest surviving pairs, at most p percent of all classes.
 *   ./sqsparse one   FORM M1,M2,.. n [n ...]                     every pair of n, with a witness
 *
 * Build: gcc -O2 -march=native -fopenmp -o sqsparse sqsparse.c -lm
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <omp.h>

typedef unsigned __int128 u128;
typedef uint64_t u64;
typedef uint32_t u32;

#define MAXM 4
static int FORM, NM;
static u64 MUL[MAXM];

static inline int bad_prime(u64 q) { return FORM == 1 ? (q & 3) == 3 : ((q & 7) == 5 || (q & 7) == 7); }
static inline int good_odd(u128 v) { u64 w = (u64)v; return FORM == 1 ? (w & 3) == 1 : ((w & 7) == 1 || (w & 7) == 3); }

/* ---------------------------------------------------------------- 128-bit modular arithmetic */
static inline void mul128(u128 a, u128 b, u128 *hi, u128 *lo) {
    u64 a0 = (u64)a, a1 = (u64)(a >> 64), b0 = (u64)b, b1 = (u64)(b >> 64);
    u128 p00 = (u128)a0 * b0, p01 = (u128)a0 * b1, p10 = (u128)a1 * b0, p11 = (u128)a1 * b1;
    u128 mid = (p00 >> 64) + (u64)p01 + (u64)p10;
    *lo = (u128)(u64)p00 | (mid << 64);
    *hi = p11 + (p01 >> 64) + (p10 >> 64) + (mid >> 64);
}
static inline u128 inv128(u128 m) {            /* m odd: inverse modulo 2^128 */
    u128 x = m;                                 /* correct to 3 bits */
    for (int i = 0; i < 7; i++) x *= 2 - m * x;
    return x;
}
/* Montgomery product a*b/2^128 mod m, for odd m < 2^127 and a, b < m */
static inline u128 mmul(u128 a, u128 b, u128 m, u128 minv) {
    u128 hi, lo, qh, ql;
    mul128(a, b, &hi, &lo);
    mul128(lo * minv, m, &qh, &ql);
    u128 r = hi - qh;
    if (hi < qh) r += m;
    return r;
}
static inline u128 addmod(u128 a, u128 b, u128 m) { u128 s = a + b; if (s >= m) s -= m; return s; }

static int sprp(u128 n, u128 ninv, u128 one, unsigned base) {      /* strong probable prime test */
    u128 b = 0;
    for (unsigned i = 0; i < base; i++) b = addmod(b, one, n);     /* base in Montgomery form */
    if (b == 0) return 1;
    u128 d = n - 1; int s = 0;
    while (!(d & 1)) { d >>= 1; s++; }
    u128 x = one;
    int top = 127;
    while (!((d >> top) & 1)) top--;
    for (int i = top; i >= 0; i--) {
        x = mmul(x, x, n, ninv);
        if ((d >> i) & 1) x = mmul(x, b, n, ninv);
    }
    u128 mone = n - one;
    if (x == one || x == mone) return 1;
    for (int i = 1; i < s; i++) {
        x = mmul(x, x, n, ninv);
        if (x == mone) return 1;
        if (x == one) return 0;
    }
    return 0;
}
static const unsigned BASES[12] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37};
/* n odd, n > 37, n < 2^127.  Deterministic for n < 3.18e23 (Sorenson and Webster, 2017). */
static int is_prime128(u128 n) {
    u128 ninv = inv128(n), one = (0 - n) % n;
    for (int i = 0; i < 12; i++) if (!sprp(n, ninv, one, BASES[i])) return 0;
    return 1;
}
static inline int ctz128(u128 x) { u64 lo = (u64)x; return lo ? __builtin_ctzll(lo) : 64 + __builtin_ctzll((u64)(x >> 64)); }
static u128 gcd128(u128 a, u128 b) {
    if (a == 0) return b;
    if (b == 0) return a;
    int sh = ctz128(a | b);
    a >>= ctz128(a);
    do { b >>= ctz128(b); if (a > b) { u128 t = a; a = b; b = t; } b -= a; } while (b);
    return a << sh;
}
static u128 rho(u128 n) {                       /* a proper factor of the odd composite n */
    u128 ninv = inv128(n), one = (0 - n) % n;
    for (unsigned c = 1;; c++) {
        u128 cc = 0;
        for (unsigned i = 0; i < c; i++) cc = addmod(cc, one, n);
        u128 y = addmod(one, one, n), x = y, ys = y, q = one, g = 1;
        u64 r = 1;
        const u64 m = 128;
        do {
            x = y;
            for (u64 i = 0; i < r; i++) y = addmod(mmul(y, y, n, ninv), cc, n);
            u64 k = 0;
            do {
                ys = y;
                u64 lim = (m < r - k) ? m : r - k;
                for (u64 i = 0; i < lim; i++) {
                    y = addmod(mmul(y, y, n, ninv), cc, n);
                    q = mmul(q, x > y ? x - y : y - x, n, ninv);
                }
                g = gcd128(q, n);
                k += m;
            } while (k < r && g == 1);
            r <<= 1;
        } while (g == 1);
        if (g == n) {
            g = 1;
            do {
                ys = addmod(mmul(ys, ys, n, ninv), cc, n);
                g = gcd128(x > ys ? x - ys : ys - x, n);
            } while (g == 1);
        }
        if (g != n) return g;
    }
}
static void factor_rec(u128 v, u128 *fac, int *nf) {     /* v odd, no prime factor below 256 */
    if (v == 1) return;
    if (v < 65536 || is_prime128(v)) { fac[(*nf)++] = v; return; }
    u128 d = rho(v);
    factor_rec(d, fac, nf);
    factor_rec(v / d, fac, nf);
}

static u32 SMALLP[64]; static int NSMALL;      /* odd primes below 256 */

/* Is v represented by the form?  If not and wit != NULL: *wit = a bad prime dividing v to an odd power,
 * or 0 if the odd part of v is in a forbidden class (3 mod 4, resp. 5 or 7 mod 8). */
static int is_rep(u128 v, u128 *wit) {
    if (v == 0) return 1;
    v >>= ctz128(v);
    if (!good_odd(v)) { if (wit) *wit = 0; return 0; }
    u64 hi = (u64)(v >> 64), lo = (u64)v;
    for (int i = 0; i < NSMALL; i++) {
        u64 p = SMALLP[i];
        u64 r = hi ? (u64)((((u128)(hi % p)) << 64 | lo) % p) : lo % p;
        if (r == 0) {
            int e = 0;
            do { v /= p; e++; } while (v % p == 0);
            if ((e & 1) && bad_prime(p)) { if (wit) *wit = p; return 0; }
            hi = (u64)(v >> 64); lo = (u64)v;
        }
    }
    if (v == 1) return 1;
    if (v < 65536) { if (bad_prime((u64)v)) { if (wit) *wit = v; return 0; } return 1; }
    u128 ninv = inv128(v), one = (0 - v) % v;
    if (sprp(v, ninv, one, 2)) {
        int pr = 1;
        for (int i = 1; i < 12 && pr; i++) pr = sprp(v, ninv, one, BASES[i]);
        if (pr) { if (bad_prime((u64)v)) { if (wit) *wit = v; return 0; } return 1; }
    }
    u128 fac[128]; int nf = 0;
    factor_rec(v, fac, &nf);
    for (int i = 0; i < nf; i++) {
        if (!bad_prime((u64)fac[i])) continue;
        int e = 0;
        for (int j = 0; j < nf; j++) if (fac[j] == fac[i]) e++;
        if (e & 1) { if (wit) *wit = fac[i]; return 0; }
    }
    return 1;
}

/* ---------------------------------------------------------------- tables */
typedef struct { u32 q, r, j, inv; u64 mask; } ent_t;   /* class r (mod q) kills the pairs `mask` of family j */
static ent_t *ENT; static long NENT;
static int KTOP[MAXM];                          /* largest k used anywhere in the run */

/* trial primes for the quick test of stage 2: the odd primes below 2048 */
#define MAXTP 320
static u32 TP[MAXTP]; static int NTP;
static uint16_t TM[MAXM][64][MAXTP];            /* m*4^k mod p */
static int STRICT;                              /* 1: 12 bases also in the quick test */
static int STATS;                               /* 1: for every flagged n count the represented alive pairs */
static u64 STAB[64][64];                        /* STAB[A][r]: flagged n with A alive pairs, r of them represented */

/* Quick decision for v = n^2 - m*4^k.  n2r[i] = n^2 mod TP[i], tm = TM[j][k].
 * Returns 1 (represented), 0 (not represented) or -1 (cofactor composite: undecided). */
static int quick_rep(u128 v, const u32 *n2r, const uint16_t *tm) {
    v >>= ctz128(v);
    if (!good_odd(v)) return 0;
    for (int i = 0; i < NTP; i++) {
        if (n2r[i] != tm[i]) continue;
        u64 p = TP[i]; int e = 0;
        while (v % p == 0) { v /= p; e++; }
        if ((e & 1) && bad_prime(p)) return 0;
    }
    if (v == 1) return 1;
    if (v < 2048ULL * 2048ULL) return !bad_prime((u64)v);
    u128 ninv = inv128(v), one = (0 - v) % v;
    if (!sprp(v, ninv, one, 2)) return -1;
    int nb = STRICT ? 12 : 2;
    for (int i = 1; i < nb; i++) if (!sprp(v, ninv, one, BASES[i])) return -1;
    return !bad_prime((u64)v);
}

static u64 powmod64(u64 a, u64 e, u64 m) { u64 r = 1; a %= m; while (e) { if (e & 1) r = (u128)r * a % m; a = (u128)a * a % m; e >>= 1; } return r; }
static u64 sqrtmod(u64 a, u64 p) {              /* p odd prime, a a non-zero square mod p */
    if ((p & 3) == 3) return powmod64(a, (p + 1) / 4, p);
    u64 q = p - 1; int s = 0;
    while (!(q & 1)) { q >>= 1; s++; }
    u64 z = 2;
    while (powmod64(z, (p - 1) / 2, p) != p - 1) z++;
    u64 c = powmod64(z, q, p), x = powmod64(a, (q + 1) / 2, p), t = powmod64(a, q, p);
    int m = s;
    while (t != 1) {
        int i = 0; u64 tt = t;
        while (tt != 1) { tt = (u128)tt * tt % p; i++; }
        u64 b = c;
        for (int j = 0; j < m - i - 1; j++) b = (u128)b * b % p;
        x = (u128)x * b % p; c = (u128)b * b % p; t = (u128)t * c % p; m = i;
    }
    return x;
}
static int kmax_of(u64 n, u64 m) {              /* largest k with m*4^k < n^2, or -1 */
    u128 n2 = (u128)n * n; int k = -1;
    while (k + 1 < 62 && ((u128)m << (2 * (k + 1))) < n2) k++;
    return k;
}
static void build_tables(u64 nhi, u32 B) {
    char *comp = calloc(B + 1, 1);
    for (u32 i = 2; (u64)i * i <= B; i++) if (!comp[i]) for (u32 j = i * i; j <= B; j += i) comp[j] = 1;
    NSMALL = 0;
    for (u32 p = 3; p < 256; p += 2) if (!comp[p]) SMALLP[NSMALL++] = p;
    NTP = 0;
    for (u32 p = 3; p < 2048 && p <= B; p += 2) if (!comp[p]) TP[NTP++] = p;
    for (int j = 0; j < NM; j++)
        for (int i = 0; i < NTP; i++) {
            u64 t = MUL[j] % TP[i];
            for (int k = 0; k < 64; k++) { TM[j][k][i] = (uint16_t)t; t = t * 4 % TP[i]; }
        }
    for (int j = 0; j < NM; j++) KTOP[j] = kmax_of(nhi, MUL[j]);
    long cap = 1 << 16; ENT = malloc(cap * sizeof(ent_t)); NENT = 0;
    for (u32 q = 3; q <= B; q += 2) {
        if (comp[q] || !bad_prime(q)) continue;
        for (int j = 0; j < NM; j++) {
            u64 a = MUL[j] % q;
            if (a == 0 || powmod64(a, (q - 1) / 2, q) != 1) continue;
            u64 s = sqrtmod(a, q);
            u32 res[130]; u64 msk[130]; int nc = 0;
            u64 r = s;
            for (int k = 0; k <= KTOP[j]; k++) {
                u32 two[2] = {(u32)r, (u32)(q - r)};
                for (int t = 0; t < 2; t++) {
                    int f = -1;
                    for (int c = 0; c < nc; c++) if (res[c] == two[t]) { f = c; break; }
                    if (f < 0) { f = nc++; res[f] = two[t]; msk[f] = 0; }
                    msk[f] |= 1ULL << k;
                }
                r = r * 2 % q;
            }
            for (int c = 0; c < nc; c++) {
                if (NENT == cap) { cap *= 2; ENT = realloc(ENT, cap * sizeof(ent_t)); }
                ENT[NENT].q = q; ENT[NENT].r = res[c]; ENT[NENT].mask = msk[c]; ENT[NENT].j = j; ENT[NENT].inv = 0; NENT++;
            }
        }
    }
    free(comp);
}

/* ---------------------------------------------------------------- one n in detail */
static void explain(u64 n) {
    u128 n2 = (u128)n * n; int reps = 0, tot = 0;
    printf("n = %llu\n", (unsigned long long)n);
    for (int j = 0; j < NM; j++) {
        int km = kmax_of(n, MUL[j]);
        for (int k = 0; k <= km; k++) {
            u128 v = n2 - ((u128)MUL[j] << (2 * k)), w = 0;
            int rp = is_rep(v, &w);
            tot++; reps += rp;
            printf("  m=%llu k=%d : %s", (unsigned long long)MUL[j], k, rp ? "REPRESENTED" : "no");
            if (!rp) {
                if (w == 0) printf("  (odd part in a forbidden class)");
                else if (w >> 64) printf("  (bad prime %llu*2^64+%llu to an odd power)", (unsigned long long)(w >> 64), (unsigned long long)(u64)w);
                else printf("  (bad prime %llu to an odd power)", (unsigned long long)(u64)w);
            }
            printf("\n");
        }
    }
    printf("  => %d of %d pairs represented%s\n", reps, tot, reps ? "" : "   *** NO REPRESENTATION ***");
}

/* ---------------------------------------------------------------- stage 2 for one n */
typedef struct { u64 hist[160], flag, over, test, full, found; } tstat;
static int TT;                                  /* threshold on the number of alive pairs */

static void handle_n(u64 n, u64 *alive, const int *km, int A, tstat *st) {
    if (A > TT + 2 * NM) { st->hist[A < 159 ? A : 159]++; return; }
    u128 n2 = (u128)n * n;
    /* the pairs with k <= 1 may fail at the prime 2: decide that part exactly */
    for (int j = 0; j < NM; j++)
        for (int k = 0; k <= 1 && k <= km[j]; k++)
            if ((alive[j] >> k) & 1) {
                u128 v = n2 - ((u128)MUL[j] << (2 * k));
                v >>= ctz128(v);
                if (!good_odd(v)) { alive[j] &= ~(1ULL << k); A--; }
            }
    st->hist[A < 159 ? A : 159]++;
    if (A > TT) return;
    st->flag++;
    if (STATS) {
        int r = 0;
        for (int j = 0; j < NM; j++) {
            u64 a = alive[j];
            while (a) { int k = __builtin_ctzll(a); a &= a - 1; r += is_rep(n2 - ((u128)MUL[j] << (2 * k)), NULL); }
        }
        if (A < 64 && r < 64) {
            #pragma omp atomic
            STAB[A][r]++;
        }
        return;
    }
    int rep = 0;
    u32 n2r[MAXTP];
    for (int t = 0; t < NTP; t++) { u64 r = n % TP[t]; n2r[t] = (u32)(r * r % TP[t]); }
    u64 undecided[MAXM];
    for (int j = 0; j < NM; j++) undecided[j] = 0;
    for (int j = 0; j < NM && !rep; j++) {      /* pass 1: quick test */
        u64 a = alive[j];
        while (a) {
            int k = __builtin_ctzll(a); a &= a - 1;
            st->test++;
            int qr = quick_rep(n2 - ((u128)MUL[j] << (2 * k)), n2r, TM[j][k]);
            if (qr == 1) { rep = 1; break; }
            if (qr < 0) undecided[j] |= 1ULL << k;
        }
    }
    for (int j = 0; j < NM && !rep; j++) {      /* pass 2: full factorization */
        u64 a = undecided[j];
        while (a) {
            int k = __builtin_ctzll(a); a &= a - 1;
            st->full++;
            if (is_rep(n2 - ((u128)MUL[j] << (2 * k)), NULL)) { rep = 1; break; }
        }
    }
    if (rep) return;
    int cnt = 0;                                /* no alive pair is represented: decide ALL pairs exactly */
    for (int j = 0; j < NM; j++)
        for (int k = 0; k <= km[j]; k++) cnt += is_rep(n2 - ((u128)MUL[j] << (2 * k)), NULL);
    if (cnt) { st->over++; return; }
    st->found++;
    #pragma omp critical
    {
        printf("NOREP %llu alive=%d\n", (unsigned long long)n, A);
        fflush(stdout);
    }
}

#define LBLK 16384                              /* sieve positions per block */

static void report(const tstat *tot, double secs, const char *what) {
    u64 total = 0;
    for (int a = 0; a < 160; a++) total += tot->hist[a];
    printf("DONE %s: %llu, flagged (alive <= %d): %llu, quick pair tests: %llu, full factorizations: %llu, sieve-only strikes undone: %llu, NOREP found: %llu, %.0fs\n",
           what, (unsigned long long)total, TT, (unsigned long long)tot->flag, (unsigned long long)tot->test,
           (unsigned long long)tot->full, (unsigned long long)tot->over, (unsigned long long)tot->found, secs);
    if (STATS) {
        printf("STATS rows: alive pairs A; columns: number r of represented alive pairs (r = 0, 1, 2, ...); last: mean r / A\n");
        for (int a = 0; a < 64; a++) {
            u64 t = 0, sr = 0; int rmax = 0;
            for (int r = 0; r < 64; r++) if (STAB[a][r]) { t += STAB[a][r]; sr += (u64)r * STAB[a][r]; rmax = r; }
            if (!t) continue;
            printf("A=%2d n=%9llu :", a, (unsigned long long)t);
            for (int r = 0; r <= rmax; r++) printf(" %llu", (unsigned long long)STAB[a][r]);
            printf("  | rho=%.4f\n", a ? (double)sr / ((double)t * a) : 0.0);
        }
    }
    printf("HIST alive:count");
    for (int a = 0; a < 160; a++) if (tot->hist[a]) printf(" %d:%llu", a, (unsigned long long)tot->hist[a]);
    printf("\n");
}
static void merge(tstat *tot, const tstat *l) {
    for (int a = 0; a < 160; a++) tot->hist[a] += l->hist[a];
    tot->flag += l->flag; tot->over += l->over; tot->test += l->test; tot->full += l->full; tot->found += l->found;
}

/* ---------------------------------------------------------------- scan: every odd n */
static void run_scan(u64 nlo, u64 nhi, int threads) {
    long super = 4096;                          /* blocks per work unit */
    {
        u64 blocks = (nhi - nlo) / (2ULL * LBLK) + 1;
        while (super > 1 && blocks / super < (u64)threads * 8) super /= 2;
    }
    u64 span = 2ULL * LBLK * super;
    long nunits = (long)((nhi - nlo + span - 1) / span), done_units = 0;
    tstat tot; memset(&tot, 0, sizeof tot);
    double t0 = omp_get_wtime();
    #pragma omp parallel num_threads(threads)
    {
        u64 *kill[MAXM];
        for (int j = 0; j < NM; j++) kill[j] = malloc(LBLK * sizeof(u64));
        u32 *nxt = malloc(NENT * sizeof(u32));
        tstat st; memset(&st, 0, sizeof st);
        #pragma omp for schedule(dynamic, 1)
        for (long u = 0; u < nunits; u++) {
            u64 ustart = nlo + (u64)u * span;   /* odd */
            for (long e = 0; e < NENT; e++) {
                u64 q = ENT[e].q;
                u64 d = (ENT[e].r + q - ustart % q) % q;
                nxt[e] = (u32)(d * ((q + 1) / 2) % q);          /* first i with ustart + 2i = r (mod q) */
            }
            for (long b = 0; b < super; b++) {
                u64 n0 = ustart + 2ULL * LBLK * b;
                if (n0 >= nhi) break;
                for (int j = 0; j < NM; j++) memset(kill[j], 0, LBLK * sizeof(u64));
                for (long e = 0; e < NENT; e++) {
                    u32 i = nxt[e], q = ENT[e].q; u64 mk = ENT[e].mask; u64 *kj = kill[ENT[e].j];
                    while (i < LBLK) { kj[i] |= mk; i += q; }
                    nxt[e] = i - LBLK;
                }
                u64 nlast = n0 + 2ULL * (LBLK - 1);
                int k0[MAXM], k1[MAXM], same = 1;
                for (int j = 0; j < NM; j++) { k0[j] = kmax_of(n0, MUL[j]); k1[j] = kmax_of(nlast, MUL[j]); if (k0[j] != k1[j]) same = 0; }
                for (int i = 0; i < LBLK; i++) {
                    u64 n = n0 + 2ULL * i;
                    if (n >= nhi) break;
                    u64 alive[MAXM]; int A = 0, km[MAXM];
                    for (int j = 0; j < NM; j++) {
                        km[j] = same ? k0[j] : kmax_of(n, MUL[j]);
                        u64 cand = km[j] < 0 ? 0 : ((2ULL << km[j]) - 1);
                        alive[j] = cand & ~kill[j][i];
                        A += __builtin_popcountll(alive[j]);
                    }
                    handle_n(n, alive, km, A, &st);
                }
            }
            #pragma omp critical
            {
                done_units++;
                if (done_units % 16 == 0 || done_units == nunits)
                    fprintf(stderr, "units %ld/%ld  (n < about %.4e)  %.0fs\n", done_units, nunits,
                            (double)nlo + (double)done_units * (double)span, omp_get_wtime() - t0);
            }
        }
        #pragma omp critical
        merge(&tot, &st);
    }
    report(&tot, omp_get_wtime() - t0, "odd n scanned");
}

/* ---------------------------------------------------------------- cls: selected residue classes */
static void run_cls(u64 nlo, u64 nhi, int threads, char *qlist, int a0max) {
    u32 qs[12]; int nqs = 0; u64 M0 = 1;
    for (char *t = strtok(qlist, ","); t && nqs < 12; t = strtok(NULL, ",")) { qs[nqs] = (u32)atoi(t); M0 *= qs[nqs]; nqs++; }
    u64 S = 2 * M0;
    /* move the classes of the primes Q_i out of the sieve table into per-residue tables */
    u64 *km0[12];
    for (int i = 0; i < nqs; i++) km0[i] = calloc((size_t)qs[i] * NM, sizeof(u64));
    long w = 0;
    for (long e = 0; e < NENT; e++) {
        int hit = -1;
        for (int i = 0; i < nqs; i++) if (ENT[e].q == qs[i]) hit = i;
        if (hit >= 0) km0[hit][(size_t)ENT[e].r * NM + ENT[e].j] |= ENT[e].mask;
        else { ENT[w] = ENT[e]; ENT[w].inv = (u32)powmod64(S % ENT[e].q, ENT[e].q - 2, ENT[e].q); w++; }
    }
    NENT = w;
    u64 top[MAXM];
    for (int j = 0; j < NM; j++) top[j] = KTOP[j] < 0 ? 0 : ((2ULL << KTOP[j]) - 1);
    /* pass over all residues r mod M0 */
    u64 h0[160] = {0};
    long ncls = 0, cap = 1 << 20;
    u64 *cls = malloc(cap * sizeof(u64));
    int dry = (a0max == 0);
    for (int pass = 0; pass < 2; pass++) {
        u32 res[12] = {0};
        for (u64 r = 0; r < M0; r++) {
            int a0 = 0;
            for (int j = 0; j < NM; j++) {
                u64 mk = 0;
                for (int i = 0; i < nqs; i++) mk |= km0[i][(size_t)res[i] * NM + j];
                a0 += __builtin_popcountll(top[j] & ~mk);
            }
            if (pass == 0) h0[a0 < 159 ? a0 : 159]++;
            else if (a0 <= a0max) {
                if (ncls == cap) { cap *= 2; cls = realloc(cls, cap * sizeof(u64)); }
                cls[ncls++] = r;
            }
            for (int i = 0; i < nqs; i++) if (++res[i] == qs[i]) res[i] = 0;
        }
        if (pass == 0 && a0max < 0) {           /* A0 = -p: keep the best classes, at most p percent of them */
            double target = -a0max / 100.0 * (double)M0, cs = 0; int sel = -1;
            for (int a = 0; a < 160; a++) {
                if (!h0[a]) continue;
                if (sel >= 0 && cs + (double)h0[a] > target) break;
                cs += (double)h0[a]; sel = a;
            }
            a0max = sel;
        }
        if (dry) break;
    }
    {
        /* model weight of a class with a0 surviving pairs: (1 - f (1 - x))^a0, f = share of pairs that
           survive the other primes, x = chance that a surviving pair is not represented */
        double f = 0.6, x = 0.27, wt = 0, ws = 0, cs = 0;
        for (int a = 0; a < 160; a++) wt += (double)h0[a] * pow(1 - f * (1 - x), a);
        printf("CLASSES modulus %llu; pairs surviving the primes of the modulus: value:classes cumulative-share-of-classes cumulative-share-of-model-weight\n", (unsigned long long)M0);
        for (int a = 0; a < 160; a++) {
            if (!h0[a]) continue;
            cs += (double)h0[a]; ws += (double)h0[a] * pow(1 - f * (1 - x), a);
            printf("  %3d:%llu %.5f %.4f%s\n", a, (unsigned long long)h0[a], cs / (double)M0, ws / wt, a == a0max ? "   <-- kept up to here" : "");
        }
        fflush(stdout);
    }
    if (dry) return;
    fprintf(stderr, "%ld of %llu classes kept, %ld sieve classes of other primes\n", ncls, (unsigned long long)M0, NENT);
    tstat tot; memset(&tot, 0, sizeof tot);
    double t0 = omp_get_wtime();
    long done = 0;
    #pragma omp parallel num_threads(threads)
    {
        u64 *kill[MAXM];
        for (int j = 0; j < NM; j++) kill[j] = malloc(LBLK * sizeof(u64));
        u32 *nxt = malloc(NENT * sizeof(u32));
        tstat st; memset(&st, 0, sizeof st);
        #pragma omp for schedule(dynamic, 4)
        for (long c = 0; c < ncls; c++) {
            u64 r = cls[c];
            u64 rp = (r & 1) ? r : r + M0;      /* the odd representative modulo 2*M0 */
            u64 mk0[MAXM];
            for (int j = 0; j < NM; j++) {
                mk0[j] = 0;
                for (int i = 0; i < nqs; i++) mk0[j] |= km0[i][(size_t)(r % qs[i]) * NM + j];
            }
            u64 tlo = nlo > rp ? (nlo - rp + S - 1) / S : 0;
            u64 thi = nhi > rp ? (nhi - rp + S - 1) / S : 0;      /* n = rp + S*t < nhi  iff  t < thi */
            if (tlo >= thi) continue;
            {
                u64 nb = rp + S * tlo, lastq = 0, nbq = 0;
                for (long e = 0; e < NENT; e++) {
                    u64 q = ENT[e].q;
                    if (q != lastq) { nbq = nb % q; lastq = q; }
                    u64 d = (ENT[e].r + q - nbq) % q;
                    nxt[e] = (u32)(d * ENT[e].inv % q);          /* first i with nb + S*i = r (mod q) */
                }
            }
            for (u64 tb = tlo; tb < thi; tb += LBLK) {
                u64 n0 = rp + S * tb;
                int len = (thi - tb < LBLK) ? (int)(thi - tb) : LBLK;
                for (int j = 0; j < NM; j++) { u64 *kj = kill[j], v = mk0[j]; for (int i = 0; i < LBLK; i++) kj[i] = v; }
                for (long e = 0; e < NENT; e++) {
                    u32 i = nxt[e], q = ENT[e].q; u64 mk = ENT[e].mask; u64 *kj = kill[ENT[e].j];
                    while (i < LBLK) { kj[i] |= mk; i += q; }
                    nxt[e] = i - LBLK;
                }
                u64 nlast = n0 + S * (u64)(len - 1);
                int k0[MAXM], k1[MAXM], same = 1;
                for (int j = 0; j < NM; j++) { k0[j] = kmax_of(n0, MUL[j]); k1[j] = kmax_of(nlast, MUL[j]); if (k0[j] != k1[j]) same = 0; }
                for (int i = 0; i < len; i++) {
                    u64 n = n0 + S * (u64)i;
                    u64 alive[MAXM]; int A = 0, km[MAXM];
                    for (int j = 0; j < NM; j++) {
                        km[j] = same ? k0[j] : kmax_of(n, MUL[j]);
                        u64 cand = km[j] < 0 ? 0 : ((2ULL << km[j]) - 1);
                        alive[j] = cand & ~kill[j][i];
                        A += __builtin_popcountll(alive[j]);
                    }
                    handle_n(n, alive, km, A, &st);
                }
            }
            #pragma omp critical
            {
                done++;
                if (done % 512 == 0 || done == ncls)
                    fprintf(stderr, "classes %ld/%ld  %.0fs\n", done, ncls, omp_get_wtime() - t0);
            }
        }
        #pragma omp critical
        merge(&tot, &st);
    }
    report(&tot, omp_get_wtime() - t0, "n sieved in the kept classes");
}

int main(int argc, char **argv) {
    if (argc < 5) {
        fprintf(stderr, "usage: %s scan|stats FORM M1,M2,.. NLO NHI T [threads] [B]\n"
                        "       %s cls FORM M1,M2,.. NLO NHI T Q1,Q2,.. A0 [threads] [B]\n"
                        "       %s one FORM M1,M2,.. n [n ...]\n", argv[0], argv[0], argv[0]);
        return 2;
    }
    FORM = atoi(argv[2]);
    NM = 0;
    for (char *t = strtok(argv[3], ","); t && NM < MAXM; t = strtok(NULL, ",")) MUL[NM++] = strtoull(t, NULL, 10);
    STATS = !strcmp(argv[1], "stats");
    if (!strcmp(argv[1], "one")) {
        build_tables(1ULL << 61, 1000);
        for (int i = 4; i < argc; i++) explain(strtoull(argv[i], NULL, 10));
        return 0;
    }
    int iscls = !strcmp(argv[1], "cls");
    if (argc < (iscls ? 9 : 7)) return 2;
    u64 nlo = (u64)strtod(argv[4], NULL), nhi = (u64)strtod(argv[5], NULL);
    TT = atoi(argv[6]);
    int ap = iscls ? 9 : 7;
    int threads = argc > ap ? atoi(argv[ap]) : omp_get_max_threads();
    u32 B = argc > ap + 1 ? (u32)atoi(argv[ap + 1]) : 65536;
    if (nlo < 3) nlo = 3;
    nlo |= 1;
    build_tables(nhi, B);
    STRICT = TT >= 128;
    fprintf(stderr, "form %d, m =", FORM);
    for (int j = 0; j < NM; j++) fprintf(stderr, " %llu", (unsigned long long)MUL[j]);
    fprintf(stderr, ", odd n in [%llu, %llu), T = %d, B = %u, %ld sieve classes, %d threads\n",
            (unsigned long long)nlo, (unsigned long long)nhi, TT, B, NENT, threads);
    if (iscls) run_cls(nlo, nhi, threads, argv[7], atoi(argv[8]));
    else run_scan(nlo, nhi, threads);
    return 0;
}
