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

namespace mitm {


static /* since it's defined in a header file */
double wtime() /* with inline it doesn't violate one definition rule */
{

  auto clock = std::chrono::high_resolution_clock::now();
  auto nanoseconds = std::chrono::duration_cast<std::chrono::nanoseconds>(clock.time_since_epoch()).count();
  double seconds = nanoseconds / (static_cast<double>(1000000000.0));
  return seconds;
}

/* deterministic RNG based on TRIVIUM */
class PRNG {
private:
    u64 s11, s12, s21, s22, s31, s32 = 1;   /* internal state */

	u64 read_urandom()
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
