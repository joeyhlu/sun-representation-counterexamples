// direct_scan.c
//
// Independent exhaustive check, by plain table look-up, of representation conjectures of the form
//
//        n = u + v + t,      u in U,  v in V,  t in T,
//
// where U and V are value sets of quadratic polynomials and T is a short list of offsets.
// No factorisation criterion and no bit-parallel trick is used: the program marks every sum
// u + v <= N in a bit table S (double loop), then for each n looks up S[n - t] for every t.
//
// Forms (all conjectures by Zhi-Wei Sun, OEIS):
//   1  A302982         n = x^2 + 5y^2 + 2^z + 3*2^w                    (x, y, z, w >= 0)
//   2  A302984         n = x^2 + 2y^2 + 2^z + 5*2^w
//   3  A308644 Conj.2  n = (3^a*4^b)^2 + c(3c+1)/2 + d(7d+3)/2         (a, b >= 0; c, d integers)
//   4  A308644 Conj.1  n = (3^a*5^b)^2 + c(3c+1)/2 + d(7d+3)/2
//   5  A308656 Conj.2  n = (2^a*9^b)^2 + c(4c+1)   + d(3d+1)/2
//   6  A302983         n = x^2 + 2y^2 + 2^z + 3*2^w                    (control: no exception expected)
//   7  A308644 Conj.1  n = (3^a*5^b)^2 + c(3c+1)/2 + d(7d+1)/2         (control)
//
// Build:  gcc -O3 -march=native -fopenmp -o direct_scan direct_scan.c -lm
// Usage:  ./direct_scan FORM N        prints every n in [1, N] with no representation
//
#define _GNU_SOURCE
#include <omp.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef uint64_t u64;
typedef int64_t i64;

static u64 N;
static u64 *S;                                   // bit table of sums u + v
static u64 *U, *V, *T;
static size_t nU, nV, nT;

static void push(u64 **arr, size_t *n, size_t *cap, u64 v) {
    if (*n == *cap) { *cap = *cap ? 2 * *cap : 1024; *arr = realloc(*arr, *cap * sizeof(u64)); }
    (*arr)[(*n)++] = v;
}
static int cmp(const void *a, const void *b) { u64 x = *(const u64 *)a, y = *(const u64 *)b; return x < y ? -1 : x > y; }
static size_t uniq(u64 *a, size_t n) {
    qsort(a, n, sizeof(u64), cmp);
    size_t m = 0;
    for (size_t i = 0; i < n; i++) if (m == 0 || a[i] != a[m - 1]) a[m++] = a[i];
    return m;
}
// values of (A*c*c + B*c)/D for all integers c (positive and negative) up to N
static void quad_values(u64 **arr, size_t *n, i64 A, i64 B, i64 D) {
    size_t cap = 0;
    *n = 0; *arr = NULL;
    for (int sign = 1; sign >= -1; sign -= 2) {
        for (i64 c = sign == 1 ? 0 : -1;; c += sign) {
            __int128 v = ((__int128)A * c * c + (__int128)B * c) / D;
            if (v < 0) { fprintf(stderr, "negative value\n"); exit(1); }
            if (v > (__int128)N) break;
            push(arr, n, &cap, (u64)v);
        }
    }
    *n = uniq(*arr, *n);
}
// offsets t = m1 * p^i + m2 * q^j  (two-power sums), or (p^i * q^j)^2 (smooth squares)
static void offsets_two_powers(u64 p, u64 m1, u64 q, u64 m2) {
    size_t cap = 0; nT = 0; T = NULL;
    for (__int128 a = 1; m1 * a <= (__int128)N; a *= p)
        for (__int128 b = 1; m1 * a + m2 * b <= (__int128)N; b *= q)
            push(&T, &nT, &cap, (u64)(m1 * a + m2 * b));
    nT = uniq(T, nT);
}
static void offsets_smooth_squares(u64 p, u64 q) {
    size_t cap = 0; nT = 0; T = NULL;
    for (__int128 a = 1; a * a <= (__int128)N; a *= p)
        for (__int128 b = a; b * b <= (__int128)N; b *= q)
            push(&T, &nT, &cap, (u64)(b * b));
    nT = uniq(T, nT);
}

int main(int argc, char **argv) {
    if (argc < 3) { fprintf(stderr, "usage: %s FORM N\n", argv[0]); return 1; }
    int form = atoi(argv[1]);
    N = (u64)strtold(argv[2], NULL);
    switch (form) {
        case 1: quad_values(&U, &nU, 1, 0, 1); quad_values(&V, &nV, 5, 0, 1); offsets_two_powers(2, 1, 2, 3); break;
        case 2: quad_values(&U, &nU, 1, 0, 1); quad_values(&V, &nV, 2, 0, 1); offsets_two_powers(2, 1, 2, 5); break;
        case 3: quad_values(&U, &nU, 3, 1, 2); quad_values(&V, &nV, 7, 3, 2); offsets_smooth_squares(3, 4); break;
        case 4: quad_values(&U, &nU, 3, 1, 2); quad_values(&V, &nV, 7, 3, 2); offsets_smooth_squares(3, 5); break;
        case 5: quad_values(&U, &nU, 4, 1, 1); quad_values(&V, &nV, 3, 1, 2); offsets_smooth_squares(2, 9); break;
        case 6: quad_values(&U, &nU, 1, 0, 1); quad_values(&V, &nV, 2, 0, 1); offsets_two_powers(2, 1, 2, 3); break;
        case 7: quad_values(&U, &nU, 3, 1, 2); quad_values(&V, &nV, 7, 1, 2); offsets_smooth_squares(3, 5); break;
        default: fprintf(stderr, "unknown form\n"); return 1;
    }
    size_t words = (size_t)(N / 64 + 1);
    S = calloc(words, 8);
    if (!S) { fprintf(stderr, "out of memory\n"); return 1; }
    fprintf(stderr, "form %d, N = %llu: |U| = %zu, |V| = %zu, |T| = %zu\n", form, (unsigned long long)N, nU, nV, nT);
    // mark all sums u + v <= N  (each thread owns a range of the table to avoid write races)
    #pragma omp parallel
    {
        int id = omp_get_thread_num(), nt = omp_get_num_threads();
        u64 lo = (u64)((__int128)(N + 1) * id / nt), hi = (u64)((__int128)(N + 1) * (id + 1) / nt);   // sums in [lo, hi)
        lo = (lo / 64) * 64; if (id + 1 < nt) hi = (hi / 64) * 64;
        for (size_t i = 0; i < nU; i++) {
            u64 u = U[i];
            if (u >= hi) break;
            // v in [lo - u, hi - u)
            size_t a = 0, b = nV;
            u64 need = lo > u ? lo - u : 0;
            while (a < b) { size_t m = (a + b) / 2; if (V[m] < need) a = m + 1; else b = m; }
            for (size_t j = a; j < nV; j++) {
                u64 s = u + V[j];
                if (s >= hi || s > N) break;
                S[s >> 6] |= 1ULL << (s & 63);
            }
        }
    }
    // scan
    u64 total = 0;
    #pragma omp parallel for schedule(dynamic, 1) reduction(+ : total)
    for (u64 blk = 0; blk <= N / (1 << 20); blk++) {
        u64 n0 = blk << 20, n1 = n0 + (1 << 20);
        if (n0 == 0) n0 = 1;
        if (n1 > N + 1) n1 = N + 1;
        for (u64 n = n0; n < n1; n++) {
            int rep = 0;
            for (size_t k = 0; k < nT && T[k] <= n; k++) {
                u64 m = n - T[k];
                if (S[m >> 6] >> (m & 63) & 1) { rep = 1; break; }
            }
            if (!rep) {
                total++;
                #pragma omp critical
                printf("%llu\n", (unsigned long long)n);
            }
        }
    }
    fprintf(stderr, "form %d: %llu numbers in [1, %llu] with no representation\n", form, (unsigned long long)total, (unsigned long long)N);
    return 0;
}
