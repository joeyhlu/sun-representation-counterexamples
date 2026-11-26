// suntarget2.c
//
// Targeted (non-exhaustive) search for odd N with NO representation
//
//        N = p + 2^x + k*2^y,     p prime,  p == c (mod 6),  x, y >= 1.
//
// (Zhi-Wei Sun's conjectures of Feb 2009: OEIS A157218/A157225 for k = 7, A157237/A157242 for
// k = 11.)  An exhaustive scan is hopeless beyond about 1e14, but exceptions are not spread
// evenly: for N in certain residue classes modulo M = 2*3*5*7*11*13*17*19 almost every
// candidate p = N - 2^x - k*2^y is forced to be divisible by a prime <= 19 or to lie in the wrong
// class mod 6.  Such N are millions of times more likely to be exceptions than a typical N.
//
// For every octave [X, 2X) the program
//   1. walks all odd residue classes r mod M and counts the "alive" pairs (x, y) of each class,
//   2. keeps the classes with the fewest alive pairs, as many as the test budget allows,
//   3. for every N = r + M*i of a kept class adds up, over the primes 23 <= q <= QMAX, how many
//      alive candidates q kills (a sum of periodic patterns, a few machine operations per N),
//   4. sends the N whose estimated number of unkilled candidates, C*exp(-kills/C), is at most U0
//      to an exact test: Miller-Rabin on the alive candidates, and, if none is prime, a
//      complete check of ALL pairs (x, y).
// Only step 4's complete check decides what is reported, so every reported N is a true
// exception; the search simply does not look at most N.
//
// Build:  gcc -O3 -march=native -fopenmp -o suntarget2 suntarget2.c -lm
// Usage:  ./suntarget2 k c Xlo Xhi scanned_per_octave outprefix [threads] [P0] [QMAX]
//
#define _GNU_SOURCE
#include <math.h>
#include <omp.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

typedef uint64_t u64;
typedef uint32_t u32;
typedef uint16_t u16;
typedef unsigned __int128 u128;

// ---------------------------------------------------------------- primality for n < 2^63
static inline u64 mont_ninv(u64 n) {              // n * x == 1 (mod 2^64), n odd
    u64 x = n;
    for (int i = 0; i < 6; i++) x *= 2 - n * x;
    return x;
}
static inline u64 mont_mul(u64 a, u64 b, u64 n, u64 ninv) {
    u128 t = (u128)a * b;
    u64 m = (u64)t * ninv;
    u64 th = (u64)(t >> 64), mh = (u64)(((u128)m * n) >> 64);
    return th >= mh ? th - mh : th - mh + n;
}
// strong probable prime test to base a (n odd, n > 2, a reduced mod n, a != 0)
static int sprp(u64 n, u64 a) {
    u64 ninv = mont_ninv(n);
    u64 one = (0 - n) % n;                         // 2^64 mod n
    u64 mone = n - one;
    u64 r2 = (u64)(((u128)one * one) % n);         // 2^128 mod n
    u64 x = mont_mul(a % n, r2, n, ninv);          // a in Montgomery form
    if (x == 0) return 1;
    u64 d = n - 1;
    int s = __builtin_ctzll(d);
    d >>= s;
    u64 y = one;
    while (d) {
        if (d & 1) y = mont_mul(y, x, n, ninv);
        x = mont_mul(x, x, n, ninv);
        d >>= 1;
    }
    if (y == one || y == mone) return 1;
    for (int r = 1; r < s; r++) {
        y = mont_mul(y, y, n, ninv);
        if (y == mone) return 1;
        if (y == one) return 0;
    }
    return 0;
}
// deterministic for all n < 2^64 (Sinclair's seven bases)
static int is_prime_det(u64 n) {
    static const u64 sp[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47};
    if (n < 2) return 0;
    for (int i = 0; i < 15; i++) {
        if (n == sp[i]) return 1;
        if (n % sp[i] == 0) return 0;
    }
    if (n < 2209) return 1;
    static const u64 bases[] = {2, 325, 9375, 28178, 450775, 9780504, 1795265022};
    for (int i = 0; i < 7; i++) {
        u64 a = bases[i] % n;
        if (a == 0) continue;
        if (!sprp(n, a)) return 0;
    }
    return 1;
}

// ---------------------------------------------------------------- problem data
static int K, CLS;
#define M1 9699690ULL                      // 2*3*5*7*11*13*17*19
static const int CP[6] = {5, 7, 11, 13, 17, 19};
#define MAXPAIRS 4200
static int NP_all, NP_base;                // pairs with s < 2X (all) and s < X (base); base pairs come first
static u64 PS[MAXPAIRS];
static unsigned char PS6[MAXPAIRS], PSC[6][MAXPAIRS];
#define MAXQ 700
static int NQ, QP[MAXQ];                   // extra primes 23..QMAX
static u16 *SMQ[MAXQ];                     // SMQ[t][j] = s_j mod q_t
static u32 INVL[MAXQ];                     // M1^{-1} mod q_t

// exact, complete: does N have a representation?
static int has_rep_full(u64 N) {
    for (int y = 1; y < 62; y++) {
        u64 ky = (u64)K << y;
        if (ky >= N) break;
        u64 rest = N - ky;
        for (int x = 1; x < 63; x++) {
            u64 px = 1ULL << x;
            if (px >= rest) break;
            u64 p = rest - px;
            if (p % 6 != (u64)CLS) continue;
            if (is_prime_det(p)) return 1;
        }
    }
    return 0;
}

static void build_pairs(u64 X) {           // base pairs: s < X; then extra pairs: X <= s < 2X
    NP_all = 0;
    for (int pass = 0; pass < 2; pass++) {
        for (int x = 1; x < 63; x++) {
            for (int y = 1; y < 62; y++) {
                u128 s = ((u128)1 << x) + ((u128)K << y);
                if (s >= (u128)2 * X) break;
                int base = s < X;
                if ((pass == 0) != base) continue;
                if (NP_all >= MAXPAIRS) { fputs("too many pairs\n", stderr); exit(1); }
                PS[NP_all++] = (u64)s;
            }
        }
        if (pass == 0) NP_base = NP_all;
    }
    for (int j = 0; j < NP_all; j++) {
        PS6[j] = (unsigned char)(PS[j] % 6);
        for (int t = 0; t < 6; t++) PSC[t][j] = (unsigned char)(PS[j] % CP[t]);
        for (int t = 0; t < NQ; t++) SMQ[t][j] = (u16)(PS[j] % QP[t]);
    }
}

// number of alive base pairs of the class r mod M1 (0 if the class is not odd)
static inline int alive_list(u64 r, int *out, int npairs) {
    int r6 = (int)(r % 6), rc[6], n = 0;
    for (int t = 0; t < 6; t++) rc[t] = (int)(r % CP[t]);
    for (int j = 0; j < npairs; j++) {
        if ((r6 + 6 - PS6[j]) % 6 != CLS) continue;
        int dead = 0;
        for (int t = 0; t < 6; t++) if (PSC[t][j] == rc[t]) { dead = 1; break; }
        if (!dead) out[n++] = j;
    }
    return n;
}

// ---------------------------------------------------------------- class enumeration (fast DFS over bitmasks)
#define MW ((MAXPAIRS + 63) / 64)
static u64 KM[6][19][MW];                  // KM[t][a]: base pairs killed when N == a mod CP[t]
static u64 A6[6][MW];                      // A6[r6]: base pairs whose p is in the right class mod 6
static long HIST[MAXPAIRS + 1];
static u64 *SEL; static long NSEL, SELCAP; static int CMAX;
static int NWB;

static void dfs(int lvl, u64 r, u64 mod, const u64 *alive, int mode) {
    if (lvl == 6) {
        int cnt = 0;
        for (int w = 0; w < NWB; w++) cnt += __builtin_popcountll(alive[w]);
        if (mode == 0) HIST[cnt]++;
        else if (cnt <= CMAX) {
            if (NSEL == SELCAP) { SELCAP = SELCAP ? 2 * SELCAP : 1 << 16; SEL = realloc(SEL, SELCAP * sizeof(u64)); }
            SEL[NSEL++] = r;
        }
        return;
    }
    int q = CP[lvl];
    u64 inv = 0;                            // mod^{-1} mod q
    for (u64 t = 1; t < (u64)q; t++) if ((mod % q) * t % q == 1) { inv = t; break; }
    for (int a = 0; a < q; a++) {
        u64 nxt[MW];
        for (int w = 0; w < NWB; w++) nxt[w] = alive[w] & ~KM[lvl][a][w];
        u64 r2 = r + mod * ((((u64)a + q - r % q) % q) * inv % q);
        dfs(lvl + 1, r2, mod * q, nxt, mode);
    }
}

static void enumerate_classes(int mode) {
    for (int r6 = 1; r6 < 6; r6 += 2) dfs(0, (u64)r6, 6, A6[r6], mode);
}

// ---------------------------------------------------------------- per-thread work
#define BL (1 << 13)
typedef struct {
    u16 acc[BL];
    u16 *pat;                               // concatenated kill patterns
    int alive[MAXPAIRS], extra[MAXPAIRS];
    u64 tested, passed, found;
    u64 res[256]; int nres;
} tls_t;
static long PATOFF[MAXQ + 1];
static int QL[MAXQ];
static double U0;                           // pass N to the exact stage if C*exp(-K/C) <= U0
typedef struct { u64 r, i0, i1; } task_t;

static inline int cand_is_prime(u64 p) {    // p odd, coprime to 2..19
    if (p < 100) return is_prime_det(p);
    if (p % 23 == 0 || p % 29 == 0 || p % 31 == 0 || p % 37 == 0 || p % 41 == 0 || p % 43 == 0 || p % 47 == 0 ||
        p % 53 == 0 || p % 59 == 0 || p % 61 == 0 || p % 67 == 0 || p % 71 == 0 || p % 73 == 0 || p % 79 == 0 ||
        p % 83 == 0 || p % 89 == 0 || p % 97 == 0) return 0;
    if (!sprp(p, 2)) return 0;
    return is_prime_det(p);
}

static void do_task(tls_t *S, const task_t *T) {
    u64 r = T->r;
    int na = alive_list(r, S->alive, NP_base);
    int ne = 0;                              // extra pairs (X <= s < 2X) that are alive for this class
    {
        int tmp[MAXPAIRS];
        int n2 = alive_list(r, tmp, NP_all);
        for (int i = 0; i < n2; i++) if (tmp[i] >= NP_base) S->extra[ne++] = tmp[i];
    }
    // kill patterns
    memset(S->pat, 0, (size_t)PATOFF[NQ] * sizeof(u16));
    for (int t = 0; t < NQ; t++) {
        u32 q = (u32)QP[t], rq = (u32)(r % q), inv = INVL[t];
        u16 *pt = S->pat + PATOFF[t];
        const u16 *sm = SMQ[t];
        for (int a = 0; a < na; a++) {
            u32 d = (u32)sm[S->alive[a]] + q - rq;
            if (d >= q) d -= q;
            pt[(u64)d * inv % q]++;
        }
        for (int rep = 1; rep < QL[t] / (int)q; rep++) memcpy(pt + (size_t)rep * q, pt, q * sizeof(u16));   // tile
    }
    // threshold on the kill count: C*exp(-K/C) <= U0  <=>  K >= C*ln(C/U0)
    int k0 = 0;
    if ((double)na > U0) k0 = (int)ceil((double)na * log((double)na / U0));
    for (u64 i0 = T->i0; i0 <= T->i1; i0 += BL) {
        long len = (T->i1 - i0 + 1 < (u64)BL) ? (long)(T->i1 - i0 + 1) : BL;
        memset(S->acc, 0, (size_t)len * sizeof(u16));
        for (int t = 0; t < NQ; t++) {
            int q = QP[t];
            const u16 *pt = S->pat + PATOFF[t];
            long pos = 0;
            int idx = (int)(i0 % (u64)q);
            while (pos < len) {
                long n = QL[t] - idx;        // the pattern is tiled to length QL[t] (a multiple of q)
                if (n > len - pos) n = len - pos;
                u16 *restrict dst = S->acc + pos;
                const u16 *restrict src = pt + idx;
                for (long u = 0; u < n; u++) dst[u] += src[u];
                pos += n;
                idx = 0;
            }
        }
        S->tested += (u64)len;
        for (long u = 0; u < len; u++) {
            if (S->acc[u] < k0) continue;
            S->passed++;
            u64 N = r + (i0 + (u64)u) * M1;
            int rep = 0;
            for (int a = 0; a < na; a++) {
                u64 p = N - PS[S->alive[a]];
                if (cand_is_prime(p)) { rep = 1; break; }
            }
            if (rep) continue;
            for (int a = 0; a < ne && !rep; a++) {
                u64 s = PS[S->extra[a]];
                if (s < N && cand_is_prime(N - s)) rep = 1;
            }
            if (rep) continue;
            if (!has_rep_full(N)) {
                S->found++;
                if (S->nres < 256) S->res[S->nres++] = N;
            }
        }
    }
}

static double now(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + 1e-9 * ts.tv_nsec;
}

int main(int argc, char **argv) {
    if (argc < 7) { fprintf(stderr, "usage: %s k c Xlo Xhi tests_per_octave outprefix [threads] [passfrac] [QMAX]\n", argv[0]); return 1; }
    K = atoi(argv[1]); CLS = atoi(argv[2]);
    u64 XLO = (u64)strtold(argv[3], NULL), XHI = (u64)strtold(argv[4], NULL);
    double budget = strtod(argv[5], NULL);
    int dry = budget < 0;                   // negative budget: only print the plan (classes, tests, model yield)
    if (dry) budget = -budget;
    const char *prefix = argv[6];
    int nth = argc > 7 ? atoi(argv[7]) : omp_get_max_threads();
    double P0 = argc > 8 ? atof(argv[8]) : 1e-11;   // smallest estimated exception probability worth an exact test
    int qmax = argc > 9 ? atoi(argv[9]) : 1024;
    if ((CLS != 1 && CLS != 5) || K < 1 || !(K & 1) || XHI > (1ULL << 62) || XLO < 1000000 || XLO >= XHI) { fputs("bad args\n", stderr); return 1; }
    omp_set_num_threads(nth);

    NQ = 0;
    PATOFF[0] = 0;
    for (int q = 23; q <= qmax && NQ < MAXQ; q += 2) {
        int pr = 1;
        for (int d = 3; d * d <= q; d += 2) if (q % d == 0) { pr = 0; break; }
        if (!pr) continue;
        QP[NQ] = q;
        SMQ[NQ] = malloc(MAXPAIRS * sizeof(u16));
        u32 inv = 0;
        for (u32 t = 1; t < (u32)q; t++) if ((M1 % q) * t % q == 1) { inv = t; break; }
        INVL[NQ] = inv;
        QL[NQ] = q * ((256 + q - 1) / q);           // tiled pattern length (kept small: all patterns must fit in L2)
        PATOFF[NQ + 1] = PATOFF[NQ] + QL[NQ];
        NQ++;
    }
    char fn_out[512], fn_log[512];
    snprintf(fn_out, sizeof fn_out, "%s.found", prefix);
    snprintf(fn_log, sizeof fn_log, "%s.octaves", prefix);
    tls_t **TL = calloc((size_t)nth, sizeof(tls_t *));
    for (int i = 0; i < nth; i++) {
        TL[i] = aligned_alloc(64, (sizeof(tls_t) + 63) / 64 * 64);
        memset(TL[i], 0, sizeof(tls_t));
        TL[i]->pat = malloc((size_t)PATOFF[NQ] * sizeof(u16));
    }
    double rho = 1.0;
    { int ps[] = {2, 3, 5, 7, 11, 13, 17, 19}; for (int i = 0; i < 8; i++) rho *= ps[i] / (ps[i] - 1.0); }
    double t0 = now(), yield_tot = 0;
    u64 grand_tested = 0, grand_found = 0;
    fprintf(stderr, "k=%d c=%d  N in [%llu, %llu)  scan budget %.3g/octave  P0=%.3g  primes 23..%d (%d)  threads=%d\n",
            K, CLS, (unsigned long long)XLO, (unsigned long long)XHI, budget, P0, QP[NQ - 1], NQ, nth);
    for (u64 X = XLO; X < XHI; X *= 2) {
        u64 Xend = 2 * X < XHI ? 2 * X : XHI;
        build_pairs(X);
        NWB = (NP_base + 63) / 64;
        {   // a candidate that survives the primes <= QMAX is prime with probability about pi_s
            double pis = 1.7810724 * log((double)QP[NQ - 1]) / log(1.5 * (double)X);
            if (pis > 0.95) pis = 0.95;
            U0 = log(1.0 / P0) / -log1p(-pis);
        }
        memset(KM, 0, sizeof KM); memset(A6, 0, sizeof A6);
        for (int j = 0; j < NP_base; j++) {
            for (int t = 0; t < 6; t++) KM[t][PSC[t][j]][j >> 6] |= 1ULL << (j & 63);
            for (int r6 = 1; r6 < 6; r6 += 2) if ((r6 + 6 - PS6[j]) % 6 == CLS) A6[r6][j >> 6] |= 1ULL << (j & 63);
        }
        memset(HIST, 0, sizeof HIST);
        enumerate_classes(0);
        // model: P(class) = (1 - rho/ln(1.5X))^C ; tests per class = (Xend - X)/M1
        double per = (double)(Xend - X) / (double)M1, lam = -log1p(-rho / log(1.5 * (double)X));
        double tot = 0, got = 0, tests = 0;
        long ncls = 0;
        for (int cc = 0; cc <= NP_base; cc++) tot += (double)HIST[cc] * per * exp(-lam * cc);
        CMAX = -1;
        for (int cc = 0; cc <= NP_base; cc++) {
            if (!HIST[cc]) continue;
            if (tests + (double)HIST[cc] * per > budget && CMAX >= 0) break;
            CMAX = cc; tests += (double)HIST[cc] * per; ncls += HIST[cc];
            got += (double)HIST[cc] * per * exp(-lam * cc);
        }
        NSEL = 0;
        enumerate_classes(1);
        u64 before_found = grand_found;
        if (dry) {
            fprintf(stderr, "[%.3g, %.3g): pairs=%d Cmax=%d classes=%ld tests=%.3g model yield %.4g of %.4g\n",
                    (double)X, (double)Xend, NP_base, CMAX, (long)NSEL, tests, got, tot);
            yield_tot += got;
            grand_tested += (u64)tests;
            continue;
        }
        // tasks: (class, chunk of at most 2^24 values of i)
        long ntask = 0, taskcap = NSEL + 1024;
        task_t *TK = malloc((size_t)taskcap * sizeof(task_t));
        for (long ci = 0; ci < NSEL; ci++) {
            u64 r = SEL[ci];
            u64 i_lo = X > r ? (X - r + M1 - 1) / M1 : 0;
            if (r + i_lo * M1 >= Xend) continue;
            u64 i_hi = (Xend - 1 - r) / M1;
            for (u64 a = i_lo; a <= i_hi; a += (1ULL << 24)) {
                u64 b = a + (1ULL << 24) - 1;
                if (b > i_hi) b = i_hi;
                if (ntask == taskcap) { taskcap *= 2; TK = realloc(TK, (size_t)taskcap * sizeof(task_t)); }
                TK[ntask].r = r; TK[ntask].i0 = a; TK[ntask].i1 = b; ntask++;
            }
        }
        #pragma omp parallel for schedule(dynamic, 1)
        for (long ti = 0; ti < ntask; ti++) do_task(TL[omp_get_thread_num()], &TK[ti]);
        free(TK);
        u64 tested = 0, passed = 0;
        FILE *fo = fopen(fn_out, "a");
        for (int i = 0; i < nth; i++) {
            tls_t *S = TL[i];
            tested += S->tested; passed += S->passed;
            for (int j = 0; j < S->nres; j++) {
                fprintf(fo, "%llu\n", (unsigned long long)S->res[j]);
                fprintf(stderr, "  NO REPRESENTATION: k=%d c=%d N = %llu\n", K, CLS, (unsigned long long)S->res[j]);
            }
            grand_found += S->found;
            S->tested = S->passed = S->found = 0; S->nres = 0;
        }
        fclose(fo);
        grand_tested += tested;
        yield_tot += got;
        FILE *fl = fopen(fn_log, "a");
        fprintf(fl, "%llu %llu %d %ld %llu %llu %.4g %.4g %llu %.0f\n", (unsigned long long)X, (unsigned long long)Xend, CMAX, (long)NSEL,
                (unsigned long long)tested, (unsigned long long)passed, got, tot, (unsigned long long)(grand_found - before_found), now() - t0);
        fclose(fl);
        fprintf(stderr, "[%.3g, %.3g): pairs=%d Cmax=%d classes=%ld scanned=%.3g exact=%.3g model yield %.3g of %.3g  found=%llu  %.0fs\n",
                (double)X, (double)Xend, NP_base, CMAX, (long)NSEL, (double)tested, (double)passed, got, tot,
                (unsigned long long)(grand_found - before_found), now() - t0);
    }
    fprintf(stderr, "DONE k=%d c=%d: tested %.4g numbers, model yield %.3g, found %llu (see %s)\n", K, CLS, (double)grand_tested, yield_tot,
            (unsigned long long)grand_found, fn_out);
    return 0;
}
