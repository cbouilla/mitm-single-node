#ifndef MITM_TYPES
#define MITM_TYPES

// this is included by nearly everything

#include <inttypes.h>

typedef uint8_t  u8 ;
typedef uint16_t u16;
typedef uint32_t u32;
typedef uint64_t u64;

#ifdef __AVX512F__
#include <immintrin.h>

typedef u32 v32 __attribute__ ((vector_size (64), aligned(64)));
typedef u64 v64 __attribute__ ((vector_size (64), aligned(64)));

static inline v32 v32bcast(u32 x) { return (v32) _mm512_set1_epi32(x); }
static inline v64 v64bcast(u64 x) { return (v64) _mm512_set1_epi64(x); }

static inline v32 v32load(const void *addr) { return (v32) _mm512_load_si512((__m512i *) addr);}
static inline void v32store(void *addr, v32 x) { _mm512_store_si512((__m512i *) addr, (__m512i) x);}
static inline v64 v64load(const void *addr) { return (v64) _mm512_load_si512((__m512i *) addr);}
static inline void v64store(void *addr, v64 x) { _mm512_store_si512((__m512i *) addr, (__m512i) x);}

static inline v32 v32zero() { return (v32) _mm512_setzero_si512(); }

// [a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p], [q,r,s,t,u,v,w,x,y,z,aa,bb,cc,dd,ee,ff] ---> [a,c,...,cc,ee], [b, d, ..., dd, ff]
static inline void v32desinterleave(v64 x, v64 y, v32 *fst, v32 *snd)
{ 
    const __m512i idx_fst = (__m512i) (v32) {0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30};
    const __m512i idx_snd = (__m512i) (v32) {1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31};
    *fst = (v32) _mm512_permutex2var_epi32((__m512i) x, idx_fst, (__m512i) y);
    *snd = (v32) _mm512_permutex2var_epi32((__m512i) x, idx_snd, (__m512i) y);
}

// [a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p], [q,r,s,t,u,v,w,x,y,z,aa,bb,cc,dd,ee,ff] ---> [a, q, b, r, ...], [..., o, ee, p, ff]
static inline void v32interleave(v32 lo, v32 hi, v64 mask, v64 *fst, v64 *snd) 
{ 
    const __m512i idx_fst = (__m512i) (v32) {0, 16, 1, 17, 2, 18, 3, 19, 4, 20, 5, 21, 6, 22, 7, 23};
    const __m512i idx_snd = (__m512i) (v32) {8, 24, 9, 25, 10, 26, 11, 27, 12, 28, 13, 29, 14, 30, 15, 31};
    __m512i x = _mm512_permutex2var_epi32((__m512i) lo, idx_fst, (__m512i) hi);
    __m512i y = _mm512_permutex2var_epi32((__m512i) lo, idx_snd, (__m512i) hi);
    *fst = (v64) x & mask;
    *snd = (v64) y & mask;
}

#else

#ifdef __AVX2__
#include <immintrin.h>

typedef u32 v32 __attribute__ ((vector_size (32), aligned(32)));
typedef u64 v64 __attribute__ ((vector_size (32), aligned(32)));
static inline v32 v32bcast(u32 x) { return (v32) _mm256_set1_epi32(x); }
static inline v64 v64bcast(u64 x) { return (v64) _mm256_set1_epi64x(x); }

static inline v32 v32load(const void *addr) { return (v32) _mm256_load_si256((__m256i *) addr);}
static inline void v32store(void *addr, v32 x) { _mm256_store_si256((__m256i *) addr, (__m256i) x);}
static inline v64 v64load(const void *addr) { return (v64) _mm256_load_si256((__m256i *) addr);}
static inline void v64store(void *addr, v64 x) { _mm256_store_si256((__m256i *) addr, (__m256i) x);}


// [a,b,c,d,e,f,g,h], [i,j,k,l,m,n,o,p] ---> [a, c, e, g, i, k, m, o], [b, d, f, h, j, l, n, p]
static inline void v32desinterleave(v64 x, v64 y, v32 *fst, v32 *snd)
{ 
    const __m256i idx = _mm256_set_epi32(7,5,3,1,6,4,2,0);
    __m256i u = _mm256_permutevar8x32_epi32((__m256i) x, idx);
    __m256i v = _mm256_permutevar8x32_epi32((__m256i) y, idx); 
    *fst = (v32) _mm256_permute2x128_si256(u, v, 0x20);
    *snd = (v32) _mm256_permute2x128_si256(u, v, 0x31); 
}

// [a,b,c,d,e,f,g,h], [i,j,k,l,m,n,o,p] ---> [a, i, b, j, c, k, d, l], [e, m, f, n, g, o, h, p]
static inline void v32interleave(v32 lo, v32 hi, v64 mask, v64 *fst, v64 *snd) 
{ 
    __m256i u = _mm256_unpacklo_epi32((__m256i) lo, (__m256i) hi);  // [a, i, b, j, e, m, f, n]
    __m256i v = _mm256_unpackhi_epi32((__m256i) lo, (__m256i) hi);  // [c, k, d, k, g, o, h, p]
    __m256i x = _mm256_permute2x128_si256(u, v, 0x20);
    __m256i y = _mm256_permute2x128_si256(u, v, 0x31);
    *fst = (v64) x & mask;
    *snd = (v64) y & mask;
}

static inline v32 v32zero() { return (v32) _mm256_setzero_si256(); }
#endif
#endif

#endif