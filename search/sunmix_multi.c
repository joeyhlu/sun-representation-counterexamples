// sunmix.c
//
// Search for counterexamples to Zhi-Wei Sun's "mixed sum" conjectures (2009):
//
//     every odd integer N (beyond a small bound, with finitely many listed exceptions)
//     can be written as  N = p + 2^x + k*2^y,   p prime,  x >= 1, y >= 1,
//     optionally with p == c (mod 6), c in {1, 5}.
//
//   OEIS A157218: k=7,  p==1 (mod 6), all odd N > 34            (checked < 1.5e8 in 2009)
//   OEIS A157225: k=7,  p==5 (mod 6), exceptions 25,31,1983     (REFUTED 2026: 1433987797)
//   OEIS A157237: k=11, p==1 (mod 6), exceptions 35,41,47,101,167,2021,119171
//   OEIS A157242: k=11, p==5 (mod 6), exceptions 33,39,45,171,247 (checked < 1e12)
//   OEIS A155860/A155904 and Sun's project k=3,5,...,61 without the congruence.
//
// METHOD (bit-parallel sumset, then exact check of survivors).
//   Odd numbers 2j+1 are indexed by j.  P = bitset of primes in the class.
//       Q[j] = OR_{x=1..X} P[j - 2^(x-1)]        <=>  2j+1 = p + 2^x          (2^x <= 2^X)
//       C[j] = OR_{y=1..Y} Q[j - k*2^(y-1)]      <=>  2j+1 = p + 2^x + k*2^y  (small x,y)
//   Every j with C[j]=0 ("survivor") is then checked EXACTLY with all pairs (x,y), using a
//   deterministic 64-bit Miller-Rabin test.  N is reported iff no representation exists.
//   The bit-parallel phase can only produce false survivors, never false negatives, so
//   correctness rests solely on has_rep().
//
// Build:  gcc -O3 -march=native -fopenmp -o sunmix sunmix.c -lm
// Usage:  ./sunmix_multi k1,k2,... cls Nlo Nhi outprefix [threads] [wlog]      (cls = 0, 1 or 5)
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
typedef unsigned __int128 u128;

#define MAXK 40
static int KS[MAXK], NK, CLS;
static u64 NLO, NHI;

// ------------------------------------------------------------ exact primality (n < 2^64)
static inline u64 mulmod(u64 a, u64 b, u64 m) { return (u64)((u128)a * b % m); }
static u64 powmod(u64 a, u64 e, u64 m) {
    u64 r = 1;
    a %= m;
    while (e) {
        if (e & 1) r = mulmod(r, a, m);
        a = mulmod(a, a, m);
        e >>= 1;
    }
    return r;
}
static int is_prime64(u64 n) {
    static const u64 sp[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97};
    if (n < 2) return 0;
    for (int i = 0; i < 25; i++) {
        if (n == sp[i]) return 1;
        if (n % sp[i] == 0) return 0;
    }
    if (n < 97 * 97) return 1;
    u64 d = n - 1;
    int s = 0;
    while (!(d & 1)) { d >>= 1; s++; }
    for (int i = 0; i < 12; i++) {             // bases 2..37: deterministic for n < 3.3e24
        u64 x = powmod(sp[i], d, n);
        if (x == 1 || x == n - 1) continue;
        int comp = 1;
        for (int r = 1; r < s; r++) {
            x = mulmod(x, x, n);
            if (x == n - 1) { comp = 0; break; }
        }
        if (comp) return 0;
    }
    return 1;
}

// exact: does N = p + 2^x + K*2^y have a solution (x,y >= 1, p prime in class)?
static int has_rep(u64 N, int K) {
    for (int y = 1; y < 57; y++) {
        u64 ky = (u64)K << y;
        if (ky >= N) break;
        u64 rest = N - ky;
        for (int x = 1; x < 63; x++) {
            u64 px = 1ULL << x;
            if (px >= rest) break;
            u64 p = rest - px;
            if (CLS && p % 6 != (u64)CLS) continue;
            if (is_prime64(p)) return 1;
        }
    }
    return 0;
}

// ------------------------------------------------------------ base primes
static u32 *BP;
static size_t NBP;
static void make_base_primes(u64 limit) {      // odd primes <= limit
    u64 n = limit + 1;
    unsigned char *c = calloc(n + 1, 1);
    size_t cnt = 0;
    for (u64 i = 3; i <= limit; i += 2) {
        if (!c[i]) {
            cnt++;
            for (u64 j = i * i; j <= limit; j += 2 * i) c[j] = 1;
        }
    }
    BP = malloc((cnt + 1) * sizeof(u32));
    NBP = 0;
    for (u64 i = 3; i <= limit; i += 2) if (!c[i]) BP[NBP++] = (u32)i;
    free(c);
}

// ------------------------------------------------------------ geometry
#define BLKW (1L << 15)                 // block: 2^15 words = 2^21 bits = 256 KB
#define CHUNKW (1L << 24)               // chunk: 2^24 words = 2^30 bits (odd numbers) = 2^31 integers
static int X, YS[MAXK];                 // 2^x for x<=X, K*2^y for y<=Y handled bit-parallel
static u64 SH_S[64], SH_T[MAXK][64];    // shifts in index space
static long PREW;                       // words of history needed before a chunk
static long TW;                         // words of Q history needed

typedef struct {
    u64 *P, *Q;                         // cover global words [gwb, gwb + PREW + CHUNKW)
    u64 *nxt;                           // next multiple index per base prime
    u64 C[BLKW];
    u64 survivors, found;
    u64 *res; int *resk; size_t nres, cap;
} tls_t;

// dst[0..n) |= src bits shifted up by s;  dst word i is global word gw0+i;  src word 0 is global word sgw
static inline void or_shift(u64 *restrict dst, long gw0, long n, const u64 *restrict src, long sgw, u64 s) {
    long ws = (long)(s >> 6);
    int bs = (int)(s & 63);
    const u64 *a = src + (gw0 - ws - sgw);
    if (bs == 0) {
        for (long i = 0; i < n; i++) dst[i] |= a[i];
    } else {
        const u64 *b = a - 1;
        int rs = 64 - bs;
        for (long i = 0; i < n; i++) dst[i] |= (a[i] << bs) | (b[i] >> rs);
    }
}

static void do_chunk(tls_t *S, long cw0) {     // chunk = global words [cw0, cw0 + CHUNKW)
    long gwb = cw0 - PREW;                     // arrays start here (may be negative: virtual zeros)
    long nw = PREW + CHUNKW;
    long qstart = cw0 - TW;                    // Q needed from this global word on
    memset(S->P, 0, (size_t)nw * 8);
    memset(S->Q, 0, (size_t)nw * 8);
    long sv0 = gwb < 0 ? 0 : gwb;              // first real global word to sieve
    u64 jlo = (u64)sv0 * 64, jend = (u64)(cw0 + CHUNKW) * 64;
    u64 nmax = 2 * jend + 1;
    size_t nbp = 0;
    while (nbp < NBP && (u64)BP[nbp] * BP[nbp] <= nmax) nbp++;
    for (size_t i = 0; i < nbp; i++) {
        u64 p = BP[i];
        u64 r = (p - 1) / 2;                   // 2j+1 == 0 mod p  <=>  j == (p-1)/2 mod p
        u64 j0 = jlo + ((r + p - jlo % p) % p);
        u64 jsq = (p * p - 1) / 2;
        if (j0 < jsq) j0 = jsq;
        S->nxt[i] = j0;
    }
    // class mask: N=2j+1 == 1 mod 6 <=> j == 0 mod 3 ; == 5 mod 6 <=> j == 2 mod 3 ; word w, bit b: j mod 3 = (w+b) mod 3
    u64 cm[3];
    for (int wm = 0; wm < 3; wm++) {
        u64 m = 0;
        for (int b = 0; b < 64; b++) {
            int jm = (wm + b) % 3;
            if (CLS == 0 || (CLS == 1 && jm == 0) || (CLS == 5 && jm == 2)) m |= 1ULL << b;
        }
        cm[wm] = m;
    }
    for (long bw = sv0; bw < cw0 + CHUNKW; bw += BLKW) {
        long bn = BLKW;
        if (bw + bn > cw0 + CHUNKW) bn = cw0 + CHUNKW - bw;
        u64 *Pb = S->P + (bw - gwb);
        // --- sieve block: mark composites
        u64 bj0 = (u64)bw * 64, bj1 = (u64)(bw + bn) * 64;
        for (size_t i = 0; i < nbp; i++) {
            u64 j = S->nxt[i];
            if (j >= bj1) continue;
            u64 p = BP[i];
            for (; j < bj1; j += p) Pb[(j - bj0) >> 6] |= 1ULL << ((j - bj0) & 63);
            S->nxt[i] = j;
        }
        for (long i = 0; i < bn; i++) Pb[i] = ~Pb[i] & cm[(bw + i) % 3];
        if (bw == 0) Pb[0] &= ~1ULL;           // j=0 is N=1: not prime
        // --- Q block
        if (bw + bn > qstart) {
            u64 *Qb = S->Q + (bw - gwb);
            for (int x = 1; x <= X; x++) or_shift(Qb, bw, bn, S->P, gwb, SH_S[x]);
        }
        // --- C block and survivors, for every k
        if (bw >= cw0) {
            for (int ki = 0; ki < NK; ki++) {
                memset(S->C, 0, (size_t)bn * 8);
                for (int y = 1; y <= YS[ki]; y++) or_shift(S->C, bw, bn, S->Q, gwb, SH_T[ki][y]);
                for (long i = 0; i < bn; i++) {
                    u64 u = ~S->C[i];
                    while (u) {
                        int b = __builtin_ctzll(u);
                        u &= u - 1;
                        u64 N = 2 * ((u64)(bw + i) * 64 + (u64)b) + 1;
                        if (N < NLO || N >= NHI) continue;
                        S->survivors++;
                        if (!has_rep(N, KS[ki])) {
                            S->found++;
                            if (S->nres == S->cap) {
                                S->cap = S->cap ? 2 * S->cap : 64;
                                S->res = realloc(S->res, S->cap * sizeof(u64));
                                S->resk = realloc(S->resk, S->cap * sizeof(int));
                            }
                            S->resk[S->nres] = KS[ki];
                            S->res[S->nres++] = N;
                        }
                    }
                }
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
    if (argc < 6) { fprintf(stderr, "usage: %s k cls Nlo Nhi outprefix [threads]\n", argv[0]); return 1; }
    NK = 0;
    for (char *tok = strtok(argv[1], ","); tok && NK < MAXK; tok = strtok(NULL, ",")) KS[NK++] = atoi(tok);
    char ksig[256] = "";
    for (int i = 0; i < NK; i++) { char t[16]; snprintf(t, sizeof t, i ? "_%d" : "%d", KS[i]); strcat(ksig, t); }
    CLS = atoi(argv[2]);
    NLO = (u64)strtold(argv[3], NULL);
    NHI = (u64)strtold(argv[4], NULL);
    const char *prefix = argv[5];
    int nth = argc > 6 ? atoi(argv[6]) : omp_get_max_threads();
    if (NK < 1 || (CLS != 0 && CLS != 1 && CLS != 5) || NHI > (1ULL << 62)) { fprintf(stderr, "bad args\n"); return 1; }
    for (int i = 0; i < NK; i++) if (KS[i] < 1 || !(KS[i] & 1) || KS[i] > 127) { fprintf(stderr, "bad k\n"); return 1; }
    omp_set_num_threads(nth);

    // bit-parallel window: index shifts <= 2^wlog  (2^x <= 2^(wlog+1), K*2^y <= 2^(wlog+1)); default wlog = 25
    int wlog = argc > 7 ? atoi(argv[7]) : 25;
    if (wlog < 8 || wlog > 25) { fputs("bad wlog", stderr); return 1; }
    X = wlog + 1;
    for (int x = 1; x <= X; x++) SH_S[x] = 1ULL << (x - 1);
    long SW = (long)(SH_S[X] >> 6) + 2;
    TW = 0;
    for (int ki = 0; ki < NK; ki++) {
        YS[ki] = 0;
        for (int y = 1; ((u64)KS[ki] << (y - 1)) <= (1ULL << wlog); y++) { SH_T[ki][y] = (u64)KS[ki] << (y - 1); YS[ki] = y; }
        long tw = (long)(SH_T[ki][YS[ki]] >> 6) + 2;
        if (tw > TW) TW = tw;
    }
    PREW = ((SW + TW + BLKW - 1) / BLKW) * BLKW;   // multiple of the block size
    make_base_primes((u64)sqrtl((long double)NHI) + 2);
    fprintf(stderr, "k=%s cls=%d N in [%llu, %llu) X=%d base primes=%zu threads=%d\n", ksig, CLS,
            (unsigned long long)NLO, (unsigned long long)NHI, X, NBP, nth);

    char fn_state[512], fn_out[512];
    snprintf(fn_state, sizeof fn_state, "%s.state", prefix);
    snprintf(fn_out, sizeof fn_out, "%s.found", prefix);
    long c_first = (long)((NLO / 2) / 64 / CHUNKW);
    long c_last = (long)(((NHI + 1) / 2 + 63) / 64 / CHUNKW);       // inclusive bound on chunk index
    long c_next = c_first;
    u64 tot_surv = 0, tot_found = 0;
    double el0 = 0;
    FILE *fs = fopen(fn_state, "r");
    if (fs) {
        long a; unsigned long long b, c; double e; int cc; unsigned long long lo, hi; char kk[256];
        if (fscanf(fs, "%255s %d %llu %llu %ld %llu %llu %lf", kk, &cc, &lo, &hi, &a, &b, &c, &e) == 8 &&
            strcmp(kk, ksig) == 0 && cc == CLS && lo == NLO && hi == NHI) {
            c_next = a; tot_surv = b; tot_found = c; el0 = e;
            fprintf(stderr, "resuming at chunk %ld\n", a);
        }
        fclose(fs);
    } else {
        FILE *fo = fopen(fn_out, "w");
        if (fo) fclose(fo);
    }

    tls_t **TL = calloc((size_t)nth, sizeof(tls_t *));
    double t0 = now();
    while (c_next <= c_last) {
        long batch = nth;
        if (c_next + batch > c_last + 1) batch = c_last + 1 - c_next;
        #pragma omp parallel
        {
            int id = omp_get_thread_num();
            if (!TL[id]) {
                tls_t *S = calloc(1, sizeof(tls_t));
                S->P = aligned_alloc(64, (size_t)(PREW + CHUNKW) * 8);
                S->Q = aligned_alloc(64, (size_t)(PREW + CHUNKW) * 8);
                S->nxt = malloc((NBP + 1) * sizeof(u64));
                TL[id] = S;
            }
            #pragma omp for schedule(dynamic, 1)
            for (long c = 0; c < batch; c++) do_chunk(TL[id], (c_next + c) * CHUNKW);
        }
        FILE *fo = fopen(fn_out, "a");
        for (int i = 0; i < nth; i++) {
            tls_t *S = TL[i];
            if (!S) continue;
            for (size_t j = 0; j < S->nres; j++) {
                fprintf(fo, "%d %llu\n", S->resk[j], (unsigned long long)S->res[j]);
                if (S->res[j] > 1000)
                    fprintf(stderr, "  NO REPRESENTATION: k=%d cls=%d N = %llu\n", S->resk[j], CLS, (unsigned long long)S->res[j]);
            }
            S->nres = 0;
            tot_surv += S->survivors; tot_found += S->found;
            S->survivors = S->found = 0;
        }
        fclose(fo);
        c_next += batch;
        double el = el0 + now() - t0;
        fs = fopen(fn_state, "w");
        fprintf(fs, "%s %d %llu %llu %ld %llu %llu %.1f\n", ksig, CLS, (unsigned long long)NLO, (unsigned long long)NHI,
                c_next, (unsigned long long)tot_surv, (unsigned long long)tot_found, el);
        fclose(fs);
        long double reached = (long double)c_next * CHUNKW * 128.0L;
        fprintf(stderr, "N < %.4Le done: survivors=%llu, no-rep found=%llu, %.0fs\n", reached,
                (unsigned long long)tot_surv, (unsigned long long)tot_found, el);
    }
    fprintf(stderr, "DONE k=%s cls=%d [%llu,%llu): %llu numbers with no representation (see %s)\n", ksig, CLS,
            (unsigned long long)NLO, (unsigned long long)NHI, (unsigned long long)tot_found, fn_out);
    return 0;
}
