#ifndef MITM_TOOLS
#define MITM_TOOLS

// this is included by nearly everything

#include <cinttypes>
#include <chrono>
#include <fstream>
#include <cstdio>
#include <string>
#include <vector>
#include <optional>

using std::vector;
using std::pair;
using std::tuple;
using std::optional;
using std::nullopt;

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
    static constexpr __m512i idx_fst = (__m512i) (v32) {0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30};
    static constexpr __m512i idx_snd = (__m512i) (v32) {1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31};
    *fst = (v32) _mm512_permutex2var_epi32((__m512i) x, idx_fst, (__m512i) y);
    *snd = (v32) _mm512_permutex2var_epi32((__m512i) x, idx_snd, (__m512i) y);
}

// [a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p], [q,r,s,t,u,v,w,x,y,z,aa,bb,cc,dd,ee,ff] ---> [a, q, b, r, ...], [..., o, ee, p, ff]
static inline void v32interleave(v32 lo, v32 hi, v64 mask, v64 *fst, v64 *snd) 
{ 
    static constexpr __m512i idx_fst = (__m512i) (v32) {0, 16, 1, 17, 2, 18, 3, 19, 4, 20, 5, 21, 6, 22, 7, 23};
    static constexpr __m512i idx_snd = (__m512i) (v32) {8, 24, 9, 25, 10, 26, 11, 27, 12, 28, 13, 29, 14, 30, 15, 31};
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
    static constexpr __m256i idx = (__m256i) (v32) {7,5,3,1,6,4,2,0};
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

#ifdef __ARM_NEON
typedef u32 v32 __attribute__ ((vector_size (16), aligned(16)));
typedef u64 v64 __attribute__ ((vector_size (16), aligned(16)));
#endif

#ifdef __ALTIVEC__
typedef u32 v32 __attribute__ ((vector_size (16), aligned(16)));
typedef u64 v64 __attribute__ ((vector_size (16), aligned(16)));
#endif


namespace mitm {

u64 make_mask(int n)
{
    return (n >= 64) ? 0xffffffffffffffffull : (1ull << n) - 1;
}

double wtime() /* with inline it doesn't violate one definition rule */
{

  auto clock = std::chrono::high_resolution_clock::now();
  auto nanoseconds = std::chrono::duration_cast<std::chrono::nanoseconds>(clock.time_since_epoch()).count();
  double seconds = nanoseconds / (static_cast<double>(1000000000.0));
  return seconds;
}

// murmur64 hash functions, tailorized for 64-bit ints / Cf. Daniel Lemire
u64 murmur64(u64 h)
{
    h ^= h >> 33;
    h *= 0xff51afd7ed558ccdull;
    h ^= h >> 33;
    h *= 0xc4ceb9fe1a85ec53ull;
    h ^= h >> 33;
    return h;
}

/* deterministic RNG based on TRIVIUM */
class PRNG {
private:
    u64 s11, s12, s21, s22, s31, s32 = 1;   /* internal state */

    void setseed()
    {
        s11 = this->seed;
        s12 = 0;
        s21 = this->seq;
        s22 = 0;
        s31 = 0;
        s32 = 0x700000000000;
        for (int i = 0; i < 18; i++)  /* blank rounds */
            rand();
    }

public:
    const u64 seed, seq;

    static u64 read_urandom()
    {
        union {
            u64 value;
            char cs[sizeof(u64)];
        } u;
        std::ifstream rfin("/dev/urandom");
        rfin.read(u.cs, sizeof(u.cs));
        rfin.close();
        return u.value;
    }

    u64 rand()
    {
        u64 s66 = (s12 << 62) ^ (s11 >> 2);
        u64 s93 = (s12 << 35) ^ (s11 >> 29);
        u64 s162 = (s22 << 59) ^ (s21 >> 5);
        u64 s177 = (s22 << 44) ^ (s21 >> 20);
        u64 s243 = (s32 << 62) ^ (s31 >> 2);
        u64 s288 = (s32 << 17) ^ (s31 >> 47);
        u64 s91 = (s12 << 37) ^ (s11 >> 27);
        u64 s92 = (s12 << 36) ^ (s11 >> 28);
        u64 s171 = (s22 << 50) ^ (s21 >> 14);
        u64 s175 = (s22 << 46) ^ (s21 >> 18);
        u64 s176 = (s22 << 45) ^ (s21 >> 19);
        u64 s264 = (s32 << 41) ^ (s31 >> 23);
        u64 s286 = (s32 << 19) ^ (s31 >> 45);
        u64 s287 = (s32 << 18) ^ (s31 >> 46);
        u64 s69 = (s12 << 59) ^ (s11 >> 5);
        u64 t1 = s66 ^ s93; /* update */
        u64 t2 = s162 ^ s177;
        u64 t3 = s243 ^ s288;
        u64 z = t1 ^ t2 ^ t3;
        t1 ^= (s91 & s92) ^ s171;
        t2 ^= (s175 & s176) ^ s264;
        t3 ^= (s286 & s287) ^ s69;
        s12 = s11;    /* rotate */
        s11 = t3;
        s22 = s21;
        s21 = t1;
        s32 = s31;
        s31 = t2;
        return z;
    }

    PRNG(u64 seed, u64 seq) : seed(seed), seq(seq)  { setseed(); }
    PRNG(u64 seed) : seed(seed), seq(0) { setseed(); }
    PRNG() : seed(read_urandom()), seq(0) { setseed(); }
};

/* represent n in 4 bytes */
void human_format(u64 n, char *target)
{
    if (n < 1000) {
        sprintf(target, "%" PRId64, n);
        return;
    }
    if (n < 1000000) {
        sprintf(target, "%.1fK", n / 1e3);
        return;
    }
    if (n < 1000000000) {
        sprintf(target, "%.1fM", n / 1e6);
        return;
    }
    if (n < 1000000000000ll) {
        sprintf(target, "%.1fG", n / 1e9);
        return;
    }
    if (n < 1000000000000000ll) {
        sprintf(target, "%.1fT", n / 1e12);
        return;
    }
}

/* represent n in 4 bytes */
u64 human_parse(const std::string &_h)
{
    std::string h(_h);
    int n = h.length();
    if (h[n - 1] == 'T') {
        h.pop_back();
        return 1000000000000ll * (u64) std::stoi(h);
    }
    if (h[n - 1] == 'G') {
        h.pop_back();
        return 1000000000ll * (u64) std::stoi(h);
    }
    if (h[n - 1] == 'M') {
        h.pop_back();
        return 1000000ll * (u64) std::stoi(h);
    }
    if (h[n - 1] == 'K') {
        h.pop_back();
        return 1000ll * (u64) std::stoi(h);
    }
    return std::stoull(h);
}
}
#endif
