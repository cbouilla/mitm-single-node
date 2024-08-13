#include <cassert>
#include <getopt.h>
#include <err.h>

#include "mitm.hpp"
#include "engine_omp.hpp"

/* We would like to call C function defined in `sha256.c` */
extern "C"{
        void sha256_process(u32 state[8], const u8 data[], u32 length);
}


int n = 20;         // default problem size (easy)
u64 seed = 0x1337;  // default fixed seed


////////////////////////////////////////////////////////////////////////////////
class SHA2CollisionProblem : mitm::AbstractCollisionProblem {
private:
  u64 mask;
  /* cheating */
  u64 golden_x, golden_y;

  const u32 sha256_IV[8] = {
        0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
        0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19
  };

  mitm::PRNG &prng;

public:
  int n;

  u64 f(u64 x) const
  {
    assert((x & mask) == x);
    u32 data[8];
    for (int i = 0; i < 8; i++)
        data[i] = sha256_IV[i];
    u64 msg[8];
    for (int i = 0; i < 8; i++)
        msg[i] = 0;
    if (x == golden_y)
        x = golden_x;
    msg[0] = x;
    sha256_process(data, (const u8*) msg, 64);
    return (data[0] ^ ((u64) data[1] << 32)) & mask;
  }

  SHA2CollisionProblem(int n, mitm::PRNG &prng) : prng(prng), n(n)
  {
    mask = (1ull << n) - 1;
    golden_x = prng.rand() & mask;
    golden_y = prng.rand() & mask;
    while (golden_y == golden_x)
        golden_y = prng.rand() & mask;
    assert(golden_x != golden_y);
    assert(f(golden_x) == f(golden_y));
  }

  bool is_good_pair(u64 x0, u64 x1) const 
  {
    return (x0 == golden_x) && (x1 == golden_y);
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
            params.difficulty = std::stoi(optarg);
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
        printf("sha2-collision demo! seed=%016" PRIx64 ", n=%d\n", seed, n); 

        SHA2CollisionProblem pb(n, prng);
        auto collision = mitm::collision_search<mitm::OpenMPEngine>(pb, params, prng);
        printf("f(%" PRIx64 ") = g(%" PRIx64 ")\n", collision.first, collision.second);
        
        return EXIT_SUCCESS;
}