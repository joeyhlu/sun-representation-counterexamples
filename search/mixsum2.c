// mixsum2.c
//
// Exhaustive search for positive integers n that are NOT of the form
//
//         n = m + a + b,      m in S,  a in A,  b in B,
//
// where A and B are finite increasing lists of positive integers and the base set S is either
//   base prime  : the primes, optionally only those in given residue classes modulo M
//                 (plus a short list of extra admitted values such as 2), or
//   base x2p2y2 : the integers of the form x^2 + 2*y^2 (x, y >= 0; zero included).
// This covers several of Zhi-Wei Sun's 2018 conjectures (OEIS A303932, A302984, ...).
//
// METHOD.  Bit i of a bitset stands for the integer i.  With P = indicator of S,
//       Q = OR_{a in A, a <= W} (P << a),      C = OR_{b in B, b <= W} (Q << b).
// Every n with C[n] = 0 is then tested EXACTLY against all pairs (a, b) of any size by
// has_rep(), which decides membership in S by Miller-Rabin (primes) or by factoring (x^2+2y^2).
// The bit stage only ever sets bits for true representations, so correctness of a reported n
// rests on has_rep() alone.
//
// Build:  gcc -O3 -march=native -fopenmp -o mixsum2 mixsum2.c -lm
// Usage:  ./mixsum2 spec.txt Nlo Nhi outprefix [threads] [wlogA] [wlogB] [chunklog]
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

#define MAXT 200
static u64 A[MAXT], B[MAXT];
static int NA, NB;
static char NAME[128];
static u64 NLO, NHI;
enum { BASE_PRIME = 0, BASE_X2P2Y2 = 1 };
static int BASE = BASE_PRIME;
static int PMOD = 0;                    // residue-class restriction for primes (0 = none)
static unsigned char PALLOW[8192];
static u64 EXTRA[16];
static int NEXTRA = 0;
static int PARITY = 0;                  // 0 all n, 1 odd n only, 2 even n only
static u64 AW[MAXT], BW[MAXT];
static int NAW, NBW;

// ------------------------------------------------------------ arithmetic
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
static u64 gcd64(u64 a, u64 b) { while (b) { u64 t = a % b; a = b; b = t; } return a; }
static u64 isqrt64(u64 n) {
    u64 r = (u64)sqrtl((long double)n);
    while (r * r > n) r--;
    while ((r + 1) * (r + 1) <= n) r++;
    return r;
}
static u64 rho64(u64 n) {                      // a non-trivial factor of an odd composite n
    for (u64 c = 1;; c++) {
        u64 x = 2, y = 2, d = 1;
        while (d == 1) {
            x = (mulmod(x, x, n) + c) % n;
            y = (mulmod(y, y, n) + c) % n;
            y = (mulmod(y, y, n) + c) % n;
            d = gcd64(x > y ? x - y : y - x, n);
        }
        if (d != n) return d;
    }
}

// small primes for the exact x^2 + 2y^2 test
static u32 *SP;
static size_t NSP;
static inline int good8(u64 q) { return (q & 7) == 1 || (q & 7) == 3; }

// exact: is m of the form x^2 + 2y^2 ?   (iff every prime == 5, 7 mod 8 divides m to an even power)
static int is_x2p2y2(u64 m) {
    if (m == 0) return 1;
    u64 u = m >> __builtin_ctzll(m);
    if (!good8(u)) return 0;                   // odd part must be == 1, 3 (mod 8)
    for (size_t i = 0; i < NSP; i++) {
        u64 p = SP[i];
        if (p * p > u) return u == 1 || good8(u);          // u is 1 or prime
        if (u % p == 0) {
            int e = 0;
            do { u /= p; e++; } while (u % p == 0);
            if (!good8(p) && (e & 1)) return 0;
        }
    }
    if (u == 1) return 1;
    // u > 1 has no prime factor below SP limit, and SP limit^3 > Nhi: u is prime, p^2 or p*q
    if (is_prime64(u)) return good8(u);
    u64 r = isqrt64(u);
    if (r * r == u) return 1;
    u64 p = rho64(u), q = u / p;
    return good8(p) && good8(q);
}

static inline int in_base(u64 m) {
    if (BASE == BASE_X2P2Y2) return is_x2p2y2(m);
    for (int i = 0; i < NEXTRA; i++) if (m == EXTRA[i]) return 1;
    if (m < 3 || !(m & 1)) return 0;
    if (PMOD && !PALLOW[m % (u64)PMOD]) return 0;
    return is_prime64(m);
}

// exact: does N have a representation m + a + b ?
static int has_rep(u64 N) {
    for (int j = 0; j < NB && B[j] <= N; j++) {
        u64 rest = N - B[j];
        for (int i = 0; i < NA && A[i] <= rest; i++)
            if (in_base(rest - A[i])) return 1;
    }
    return 0;
}

// ------------------------------------------------------------ base primes
static u32 *BP;
static size_t NBP;
static void make_base_primes(u64 limit, int kind) {    // kind 0: odd primes; kind 1: primes == 5, 7 mod 8
    unsigned char *c = calloc(limit + 2, 1);
    size_t cnt = 0;
    for (u64 i = 3; i <= limit; i += 2) {
        if (!c[i]) {
            cnt++;
            for (u64 j = i * i; j <= limit; j += 2 * i) c[j] = 1;
        }
    }
    BP = malloc((cnt + 1) * sizeof(u32));
    NBP = 0;
    for (u64 i = 3; i <= limit; i += 2) if (!c[i] && (kind == 0 || (i & 7) == 5 || (i & 7) == 7)) BP[NBP++] = (u32)i;
    free(c);
}
static void make_small_primes(u64 limit) {
    unsigned char *c = calloc(limit + 2, 1);
    SP = malloc((limit / 2 + 16) * sizeof(u32));
    NSP = 0;
    for (u64 i = 2; i <= limit; i++) {
        if (!c[i]) {
            SP[NSP++] = (u32)i;
            for (u64 j = i * i; j <= limit; j += i) c[j] = 1;
        }
    }
    free(c);
}

// ------------------------------------------------------------ geometry
#define BLKW (1L << 15)                 // block: 2^15 words = 2^21 integers
static long CHW;
static long HA, HB;
static u64 *P, *Q;
static u64 *CMASK; static long CMW;     // periodic class mask for primes

typedef struct {
    u64 *nxt; u32 *tm;
    u64 C[BLKW];
    u64 survivors;
    u64 *res; size_t nres, cap;
} tls_t;

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

// word w of the indicator "odd part of i is == 1, 3 (mod 8)" (i = 64w + bit); 0 counts as good
static inline u64 gmask(u64 w) {
    u64 m = 0x0A0A0A0A0A0A0A0AULL;                                  // odd i == 1, 3 mod 8
    m |= (1ULL << 2) | (1ULL << 6) | (1ULL << 18) | (1ULL << 22) | (1ULL << 34) | (1ULL << 38) | (1ULL << 50) | (1ULL << 54);   // i == 2, 6 mod 16
    m |= (1ULL << 4) | (1ULL << 12) | (1ULL << 36) | (1ULL << 44);  // i == 4, 12 mod 32
    m |= (1ULL << 8) | (1ULL << 24);                                // i == 8, 24 mod 64
    if ((w & 1) == 0) m |= (1ULL << 16) | (1ULL << 48);             // i == 16, 48 mod 128
    if ((w & 3) <= 1) m |= 1ULL << 32;                              // i == 32, 96 mod 256
    if (w == 0) m |= 1ULL;                                          // i = 0
    else { u64 u = w >> __builtin_ctzll(w); if (good8(u)) m |= 1ULL; }   // i = 64w: odd part of w
    return m;
}

static void sieve_slab(long bw0, long bw1, long pg0, size_t nbp, u64 *nxt, u32 *tm) {
    u64 jlo = (u64)bw0 * 64;
    if (BASE == BASE_PRIME) {
        for (size_t i = 0; i < nbp; i++) {
            u64 p = BP[i];
            u64 j = ((jlo + p - 1) / p) * p;
            if (!(j & 1)) j += p;
            if (j < p * p) j = p * p;
            nxt[i] = j;
        }
    } else {
        for (size_t i = 0; i < nbp; i++) {
            u64 p = BP[i];
            u64 j = ((jlo + p - 1) / p) * p;
            if (j == 0) j = p;
            nxt[i] = j;
            tm[i] = (u32)((j / p) % p);
        }
    }
    for (long bw = bw0; bw < bw1; bw += BLKW) {
        u64 *Pb = P + (bw - pg0);
        memset(Pb, 0, (size_t)BLKW * 8);
        u64 bj0 = (u64)bw * 64, bj1 = bj0 + (u64)BLKW * 64;
        if (BASE == BASE_PRIME) {
            for (size_t i = 0; i < nbp; i++) {
                u64 j = nxt[i];
                if (j >= bj1) continue;
                u64 p2 = 2 * (u64)BP[i];
                for (; j < bj1; j += p2) Pb[(j - bj0) >> 6] |= 1ULL << ((j - bj0) & 63);
                nxt[i] = j;
            }
            for (long i = 0; i < BLKW; i++) {
                u64 v = ~Pb[i] & 0xAAAAAAAAAAAAAAAAULL;
                if (CMASK) v &= CMASK[(bw + i) % CMW];
                Pb[i] = v;
            }
            if (bw == 0) Pb[0] &= ~2ULL;
            for (int e = 0; e < NEXTRA; e++) {
                u64 v = EXTRA[e];
                if (v >= bj0 && v < bj1) Pb[(v - bj0) >> 6] |= 1ULL << ((v - bj0) & 63);
            }
        } else {
            for (size_t i = 0; i < nbp; i++) {
                u64 j = nxt[i];
                if (j >= bj1) continue;
                u64 q = BP[i];
                u32 t = tm[i];
                for (; j < bj1; j += q) {
                    int mark;
                    if (t != 0) mark = 1;                      // q divides j exactly once
                    else {
                        u64 r = j / q;
                        int v = 1;
                        while (r % q == 0) { r /= q; v++; }
                        mark = v & 1;
                    }
                    if (mark) Pb[(j - bj0) >> 6] |= 1ULL << ((j - bj0) & 63);
                    if (++t == q) t = 0;
                }
                nxt[i] = j;
                tm[i] = t;
            }
            for (long i = 0; i < BLKW; i++) Pb[i] = ~Pb[i] & gmask((u64)(bw + i));
        }
    }
}

static double now(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + 1e-9 * ts.tv_nsec;
}

static int read_spec(const char *fn) {
    FILE *f = fopen(fn, "r");
    if (!f) return 0;
    char key[64], val[64];
    int ok = 0;
    while (fscanf(f, "%63s", key) == 1) {
        if (!strcmp(key, "name")) { if (fscanf(f, "%127s", NAME) != 1) return 0; }
        else if (!strcmp(key, "base")) {
            if (fscanf(f, "%63s", val) != 1) return 0;
            if (!strcmp(val, "prime")) BASE = BASE_PRIME; else if (!strcmp(val, "x2p2y2")) BASE = BASE_X2P2Y2; else return 0;
        } else if (!strcmp(key, "parity")) {
            if (fscanf(f, "%63s", val) != 1) return 0;
            PARITY = !strcmp(val, "odd") ? 1 : !strcmp(val, "even") ? 2 : 0;
        } else if (!strcmp(key, "pmod")) {
            int n;
            if (fscanf(f, "%d %d", &PMOD, &n) != 2 || PMOD < 2 || PMOD > 8192) return 0;
            memset(PALLOW, 0, sizeof PALLOW);
            for (int i = 0; i < n; i++) { int r; if (fscanf(f, "%d", &r) != 1 || r < 0 || r >= PMOD) return 0; PALLOW[r] = 1; }
        } else if (!strcmp(key, "extra")) {
            if (fscanf(f, "%d", &NEXTRA) != 1 || NEXTRA < 0 || NEXTRA > 16) return 0;
            for (int i = 0; i < NEXTRA; i++) { unsigned long long v; if (fscanf(f, "%llu", &v) != 1) return 0; EXTRA[i] = v; }
        } else if (!strcmp(key, "A") || !strcmp(key, "B")) {
            int n;
            if (fscanf(f, "%d", &n) != 1 || n < 1 || n > MAXT) return 0;
            u64 *T = key[0] == 'A' ? A : B;
            for (int i = 0; i < n; i++) {
                unsigned long long v;
                if (fscanf(f, "%llu", &v) != 1) return 0;
                T[i] = v;
                if (v < 1 || (i && T[i] <= T[i - 1])) return 0;
            }
            if (key[0] == 'A') NA = n; else NB = n;
            ok++;
        } else return 0;
    }
    fclose(f);
    return ok == 2;
}

int main(int argc, char **argv) {
    if (argc < 5) { fprintf(stderr, "usage: %s spec Nlo Nhi outprefix [threads] [wlogA] [wlogB] [chunklog]\n", argv[0]); return 1; }
    if (!read_spec(argv[1])) { fputs("bad spec file\n", stderr); return 1; }
    NLO = (u64)strtold(argv[2], NULL);
    NHI = (u64)strtold(argv[3], NULL);
    const char *prefix = argv[4];
    int nth = argc > 5 ? atoi(argv[5]) : omp_get_max_threads();
    int wlogA = argc > 6 ? atoi(argv[6]) : 28;
    int wlogB = argc > 7 ? atoi(argv[7]) : 28;
    int chlog = argc > 8 ? atoi(argv[8]) : 26;
    if (NHI > (1ULL << 62) || NLO >= NHI || wlogA < 4 || wlogA > 30 || wlogB < 4 || wlogB > 30 || chlog < 15 || chlog > 28) {
        fputs("bad args\n", stderr); return 1;
    }
    omp_set_num_threads(nth);
    CHW = 1L << chlog;

    u64 maxA = 0, maxB = 0;
    for (int i = 0; i < NA; i++) if (A[i] <= (1ULL << wlogA)) { AW[NAW++] = A[i]; maxA = A[i]; }
    for (int i = 0; i < NB; i++) if (B[i] <= (1ULL << wlogB)) { BW[NBW++] = B[i]; maxB = B[i]; }
    HA = (((long)(maxA >> 6) + 2 + BLKW - 1) / BLKW) * BLKW;
    HB = (((long)(maxB >> 6) + 2 + BLKW - 1) / BLKW) * BLKW;

    long c_begin = (long)((NLO / 64) / BLKW) * BLKW;
    long c_end = (long)(((NHI + 63) / 64 + BLKW - 1) / BLKW) * BLKW;
    long nchunks = (c_end - c_begin + CHW - 1) / CHW;
    make_base_primes((u64)sqrtl((long double)c_end * 64.0L) + 2, BASE == BASE_X2P2Y2);
    if (BASE == BASE_X2P2Y2) make_small_primes((u64)cbrtl((long double)c_end * 64.0L) + 100);
    if (BASE == BASE_PRIME && PMOD) {
        long g = PMOD, h = 64;
        while (h) { long t = g % h; g = h; h = t; }
        CMW = PMOD / g;
        CMASK = calloc((size_t)CMW, 8);
        for (long w = 0; w < CMW; w++)
            for (int b = 0; b < 64; b++)
                if (PALLOW[(64 * w + b) % PMOD]) CMASK[w] |= 1ULL << b;
    }
    fprintf(stderr, "%s: base=%s pmod=%d extra=%d parity=%d |A|=%d |B|=%d window %d+%d  N in [%llu, %llu)  chunks=%ld threads=%d\n",
            NAME, BASE == BASE_PRIME ? "prime" : "x2p2y2", PMOD, NEXTRA, PARITY, NA, NB, NAW, NBW,
            (unsigned long long)NLO, (unsigned long long)NHI, nchunks, nth);

    P = aligned_alloc(64, (size_t)(HA + HB + CHW) * 8);
    Q = aligned_alloc(64, (size_t)(HB + CHW) * 8);
    if (!P || !Q) { fputs("out of memory\n", stderr); return 1; }

    char fn_state[512], fn_out[512];
    snprintf(fn_state, sizeof fn_state, "%s.state", prefix);
    snprintf(fn_out, sizeof fn_out, "%s.found", prefix);
    long c_next = 0;
    u64 tot_surv = 0, tot_found = 0;
    double el0 = 0;
    FILE *fs = fopen(fn_state, "r");
    int resumed = 0;
    if (fs) {
        long a; unsigned long long b, c, lo, hi; double e; char nm[128]; int wa, wb, cl;
        if (fscanf(fs, "%127s %llu %llu %d %d %d %ld %llu %llu %lf", nm, &lo, &hi, &wa, &wb, &cl, &a, &b, &c, &e) == 10 &&
            !strcmp(nm, NAME) && lo == NLO && hi == NHI && wa == wlogA && wb == wlogB && cl == chlog) {
            c_next = a; tot_surv = b; tot_found = c; el0 = e; resumed = 1;
            fprintf(stderr, "resuming at chunk %ld\n", a);
        }
        fclose(fs);
    }
    if (!resumed) { FILE *fo = fopen(fn_out, "w"); if (fo) fclose(fo); }

    tls_t **TL = calloc((size_t)nth, sizeof(tls_t *));
    for (int i = 0; i < nth; i++) {
        TL[i] = calloc(1, sizeof(tls_t));
        TL[i]->nxt = malloc((NBP + 1) * sizeof(u64));
        TL[i]->tm = malloc((NBP + 1) * sizeof(u32));
    }
    u64 pmask = PARITY == 1 ? 0xAAAAAAAAAAAAAAAAULL : PARITY == 2 ? 0x5555555555555555ULL : ~0ULL;
    double t0 = now();
    for (long ci = c_next; ci < nchunks; ci++) {
        long c0 = c_begin + ci * CHW;
        long c1 = c0 + CHW < c_end ? c0 + CHW : c_end;
        long pg0 = c0 - HA - HB, qg0 = c0 - HB;
        long ps0 = pg0 < 0 ? 0 : pg0, qs0 = qg0 < 0 ? 0 : qg0;
        if (ps0 > pg0) memset(P, 0, (size_t)(ps0 - pg0) * 8);
        if (qs0 > qg0) memset(Q, 0, (size_t)(qs0 - qg0) * 8);
        size_t nbp = 0;
        u64 nmax = (u64)c1 * 64;
        while (nbp < NBP && (u64)BP[nbp] * BP[nbp] <= nmax) nbp++;

        #pragma omp parallel
        {
            int id = omp_get_thread_num(), nt = omp_get_num_threads();
            long nb = (c1 - ps0) / BLKW;
            long b0 = nb * id / nt, b1 = nb * (id + 1) / nt;
            if (b1 > b0) sieve_slab(ps0 + b0 * BLKW, ps0 + b1 * BLKW, pg0, nbp, TL[id]->nxt, TL[id]->tm);
        }
        long nqb = (c1 - qs0) / BLKW;
        #pragma omp parallel for schedule(dynamic, 4)
        for (long bi = 0; bi < nqb; bi++) {
            long bw = qs0 + bi * BLKW;
            u64 *Qb = Q + (bw - qg0);
            memset(Qb, 0, (size_t)BLKW * 8);
            for (int i = 0; i < NAW; i++) or_shift(Qb, bw, BLKW, P, pg0, AW[i]);
        }
        long ncb = (c1 - c0) / BLKW;
        #pragma omp parallel for schedule(dynamic, 4)
        for (long bi = 0; bi < ncb; bi++) {
            tls_t *S = TL[omp_get_thread_num()];
            long bw = c0 + bi * BLKW;
            memset(S->C, 0, (size_t)BLKW * 8);
            for (int i = 0; i < NBW; i++) or_shift(S->C, bw, BLKW, Q, qg0, BW[i]);
            for (long i = 0; i < BLKW; i++) {
                u64 u = ~S->C[i] & pmask;
                while (u) {
                    int b = __builtin_ctzll(u);
                    u &= u - 1;
                    u64 N = (u64)(bw + i) * 64 + (u64)b;
                    if (N < NLO || N >= NHI) continue;
                    S->survivors++;
                    if (!has_rep(N)) {
                        if (S->nres == S->cap) { S->cap = S->cap ? 2 * S->cap : 64; S->res = realloc(S->res, S->cap * sizeof(u64)); }
                        S->res[S->nres++] = N;
                    }
                }
            }
        }
        FILE *fo = fopen(fn_out, "a");
        for (int i = 0; i < nth; i++) {
            tls_t *S = TL[i];
            for (size_t j = 0; j < S->nres; j++) {
                fprintf(fo, "%llu\n", (unsigned long long)S->res[j]);
                if (S->res[j] > 100000) fprintf(stderr, "  NO REPRESENTATION: %s  n = %llu\n", NAME, (unsigned long long)S->res[j]);
            }
            tot_found += S->nres; S->nres = 0;
            tot_surv += S->survivors; S->survivors = 0;
        }
        fclose(fo);
        double el = el0 + now() - t0;
        fs = fopen(fn_state, "w");
        fprintf(fs, "%s %llu %llu %d %d %d %ld %llu %llu %.1f\n", NAME, (unsigned long long)NLO, (unsigned long long)NHI,
                wlogA, wlogB, chlog, ci + 1, (unsigned long long)tot_surv, (unsigned long long)tot_found, el);
        fclose(fs);
        if ((ci & 7) == 7 || ci + 1 == nchunks)
            fprintf(stderr, "n < %.4Le done: survivors=%llu, no-rep found=%llu, %.0fs\n", (long double)c1 * 64.0L,
                    (unsigned long long)tot_surv, (unsigned long long)tot_found, el);
    }
    fprintf(stderr, "DONE %s [%llu,%llu): %llu numbers with no representation (see %s)\n", NAME,
            (unsigned long long)NLO, (unsigned long long)NHI, (unsigned long long)tot_found, fn_out);
    return 0;
}
