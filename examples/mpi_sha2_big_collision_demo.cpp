#include <cassert>
#include <getopt.h>
#include <err.h>

#include <mpi.h>

#include "mitm.hpp"
#include "mpi/pcs_engine.hpp"

/* We would like to call C function defined in `sha256.c` */
extern "C"{
    void sha256_process(u32 state[8], const u8 data[], u32 length);
}


int n = 20;         // default problem size (easy)
u64 seed = 0x1337;  // default fixed seed


////////////////////////////////////////////////////////////////////////////////
class SHA2CollisionProblem : public mitm::AbstractCollisionProblem {
private:
  u64 mask;
 
  const u32 sha256_IV[8] = {
        0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
        0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19
  };

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
    msg[0] = x;
    sha256_process(data, (const u8*) msg, 64);
    return (data[0] ^ ((u64) data[1] << 32)) & mask;
  }

  SHA2CollisionProblem(int n) :  mask((1ull << n) - 1), n(n) {}
};


void process_command_line_options(int argc, char **argv, mitm::MpiParameters &params)
{
    struct option longopts[5] = {
        {"ram", required_argument, NULL, 'r'},
        {"n", required_argument, NULL, 'n'},
        {"recv-per-node", required_argument, NULL, 'e'},
        {NULL, 0, NULL, 0}
    };

    for (;;) {
        int ch = getopt_long(argc, argv, "", longopts, NULL);
        switch (ch) {
        case -1:
            return;
        case 'r':
            params.nbytes_memory = mitm::human_parse(optarg);
            break;
        case 'n':
            n = std::stoi(optarg);
            break;
        case 'e':
            params.recv_per_node = std::stoi(optarg);
            break;
        default:
            errx(1, "Unknown option %s\n", optarg);
        }
    }
}


int main(int argc, char* argv[])
{
    MPI_Init(NULL, NULL);

    mitm::MpiParameters params;
    process_command_line_options(argc, argv, params);
    params.setup(MPI_COMM_WORLD);

    if (params.role == mitm::CONTROLLER)
        printf("sha2-collision demo! n=%d\n", n); 
   
    SHA2CollisionProblem pb(n);
    mitm::PRNG prng(seed);
    auto collision = mitm::collision_search<mitm::MpiEngine>(pb, params, prng);
    if (params.role == mitm::CONTROLLER)
        printf("f(%" PRIx64 ") = g(%" PRIx64 ")\n", collision.first, collision.second);
    
    MPI_Finalize();
    return EXIT_SUCCESS;
}