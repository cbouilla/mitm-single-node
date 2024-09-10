#include <cassert>
#include <getopt.h>
#include <err.h>

#include "mitm.hpp"
#include "sequential/pcs_engine.hpp"

/* We would like to call C function defined in `sha256.c` */
extern "C"{
        void sha256_process(u32 state[8], const u8 data[], u32 length);
}


int n = 20;         // default problem size (easy)
u64 seed = 0x1337;  // default fixed seed


////////////////////////////////////////////////////////////////////////////////
class SHA2ClawProblem : public mitm::AbstractClawProblem
{
public:
  int n, m;
  u64 mask;

  /* cheating */
  u64 g_shift, golden_x, golden_y;

  const u32 sha256_IV[8] = {
        0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
        0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19
  };

  u64 f(u64 x) const
  {
    assert((x & mask) == x);
    u32 data[8];
    for (int i = 0; i < 8; i++)
        data[i] = sha256_IV[i];
    u64 msg[8];
    for (int i = 0; i < 8; i++)
        msg[i] = 0;
    msg[0] = x;
    sha256_process(data, (const u8*) msg, 64);
    return (data[0] ^ ((u64) data[1] << 32)) & mask;
  }

  u64 g(u64 x) const
  {
    assert((x & mask) == x);
    u32 data[8];
    for (int i = 0; i < 8; i++)
        data[i] = sha256_IV[i];
    u64 msg[8];
    for (int i = 0; i < 8; i++)
        msg[i] = 0xffffffff;
    msg[0] = x;
    sha256_process(data, (const u8*) msg, 64);
    return (data[0] ^ ((u64) data[1] << 32) ^ g_shift) & mask;
  }

  SHA2ClawProblem(int n, mitm::PRNG &prng) : n(n), m(n)
  {
    mask = (1ull << n) - 1;
    golden_x = prng.rand() & mask;
    golden_y = prng.rand() & mask;
    u64 y = f(golden_x);
    g_shift = 0;
    g_shift = g(golden_y) ^ y;
    
    assert(f(golden_x) == g(golden_y));
  }

  bool is_good_pair(u64 x, u64 y) const 
  {
    return (x == golden_x) && (y == golden_y);
  }
};


mitm::Parameters process_command_line_options(int argc, char **argv)
{
    struct option longopts[5] = {
        {"ram", required_argument, NULL, 'r'},
        {"difficulty", required_argument, NULL, 'd'},
        {"n", required_argument, NULL, 'n'},
        {"seed", required_argument, NULL, 's'},
        {NULL, 0, NULL, 0}
    };

    mitm::Parameters params;

    for (;;) {
        int ch = getopt_long(argc, argv, "", longopts, NULL);
        switch (ch) {
        case -1:
            return params;
        case 'r':
            params.nbytes_memory = mitm::human_parse(optarg);
            break;
        case 'd':
            params.theta = std::stof(optarg);
            break;
        case 'n':
            n = std::stoi(optarg);
            break;
        case 's':
            seed = std::stoull(optarg, 0);
            break;
        default:
            errx(1, "Unknown option %s\n", optarg);
        }
    }
}


int main(int argc, char* argv[])
{
    mitm::Parameters params = process_command_line_options(argc, argv);
    mitm::PRNG prng(seed);
    printf("sha2-claw demo! seed=%016" PRIx64 ", n=%d\n", seed, n); 

    SHA2ClawProblem pb(n, prng);
    auto claw = mitm::claw_search<mitm::ScalarSequentialEngine>(pb, params, prng);
    if (claw) {
        auto [x0, x1] = *claw;
        printf("f(%" PRIx64 ") = g(%" PRIx64 ")\n", x0, x1);
    } else {
        printf("Golden collision not found\n");
    }
        
    return EXIT_SUCCESS;
}
