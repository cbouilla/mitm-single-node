#include <cassert>
#include <getopt.h>
#include <err.h>

#include <mpi.h>

#include "mpi_common.hpp"
#include "mpi_naive_alltoall.hpp"
#include "mpi_naive_isend.hpp"

/* We would like to call C function defined in `sha256.c` */
extern "C"{
    void Speck64128KeySchedule(const u32 K[], u32 rk[]);
    void Speck64128Encrypt(const u32 Pt[],u32 Ct[], const u32 rk[]);
    void Speck64128Decrypt(u32 Pt[],const u32 Ct[],const u32 rk[]);
}

int n = 20;         // default problem size (easy)
u64 seed = 0x1337;  // default fixed seed
bool expensive;

////////////////////////////////////////////////////////////////////////////////
class DoubleSpeck64_Problem : mitm::AbstractClawProblem
{
public:
    int n;
    u64 mask;
    mitm::PRNG &prng;

    u32 P[2][2] = {{0, 0}, {0xffffffff, 0xffffffff}};         /* two plaintext-ciphertext pairs */
    u32 C[2][2];

    // speck32-64 encryption of P[0], using k
    u64 f(u64 k) const
    {
        assert((k & mask) == k);
        u32 K[4] = {(u32) (k & 0xffffffff), (u32) ((k >> 32)), 0, 0};
        u32 rk[22];
        Speck64128KeySchedule(K, rk);
        u32 Ct[2];
        Speck64128Encrypt(P[0], Ct, rk);
        return ((u64) Ct[0] ^ ((u64) Ct[1] << 32)) & mask;
    }

    // speck32-64 decryption of C[0], using k
    u64 g(u64 k) const
    {
        assert((k & mask) == k);
        u32 K[4] = {(u32) (k & 0xffffffff), (u32) ((k >> 32)), 0, 0};
        u32 rk[22];
        Speck64128KeySchedule(K, rk);
        u32 Pt[2];
        Speck64128Decrypt(Pt, C[0], rk);
        return ((u64) Pt[0] ^ ((u64) Pt[1] << 32)) & mask;
    }

  DoubleSpeck64_Problem(int n, mitm::PRNG &prng) : n(n), prng(prng)
  {
    assert(n <= 64);
    mask = (1ull << n) - 1;
    u64 khi = prng.rand() & mask;
    u64 klo = prng.rand() & mask;
    // printf("Secret keys = %016" PRIx64 " %016" PRIx64 "\n", khi, klo);
    u32 Ka[4] = {(u32) (khi & 0xffffffff), (u32) ((khi >> 32)), 0, 0};
    u32 Kb[4] = {(u32) (klo & 0xffffffff), (u32) ((klo >> 32)), 0, 0};
    u32 rka[27];
    u32 rkb[27];
    Speck64128KeySchedule(Ka, rka);
    Speck64128KeySchedule(Kb, rkb);
    u32 mid[2][2];
    Speck64128Encrypt(P[0], mid[0], rka);
    Speck64128Encrypt(mid[0], C[0], rkb);
    Speck64128Encrypt(P[1], mid[1], rka);
    Speck64128Encrypt(mid[1], C[1], rkb);

    assert(f(khi) == g(klo));
    assert(is_good_pair(khi, klo));
  }

  bool is_good_pair(u64 khi, u64 klo) const
  {
    u32 Ka[4] = {(u32) (khi & 0xffffffff), (u32) ((khi >> 32)), 0, 0};
    u32 Kb[4] = {(u32) (klo & 0xffffffff), (u32) ((klo >> 32)), 0, 0};
    u32 rka[27];
    u32 rkb[27];
    Speck64128KeySchedule(Ka, rka);
    Speck64128KeySchedule(Kb, rkb);
    u32 mid[2];
    u32 Ct[2];
    Speck64128Encrypt(P[1], mid, rka);
    Speck64128Encrypt(mid, Ct, rkb);
    return (Ct[0] == C[1][0]) && (Ct[1] == C[1][1]);
  }
};


void process_command_line_options(int argc, char **argv, mitm::MpiParameters &params)
{
    struct option longopts[5] = {
        {"n", required_argument, NULL, 'n'},
        {"seed", required_argument, NULL, 's'},
        {"recv-per-node", required_argument, NULL, 'e'},
        {"expensive", no_argument, NULL, 'p'},
        {NULL, 0, NULL, 0}
    };

    for (;;) {
        int ch = getopt_long(argc, argv, "", longopts, NULL);
        switch (ch) {
        case -1:
            return;
        case 'n':
            n = std::stoi(optarg);
            break;
        case 's':
            seed = std::stoull(optarg, 0);
            break;
        case 'e':
            params.recv_per_node = std::stoi(optarg);
            break;
        case 'p':
            expensive = 1;
            break;
        default:
            errx(1, "Unknown option %s\n", optarg);
        }
    }
}


int main(int argc, char* argv[])
{
    MPI_Init(NULL, NULL);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    mitm::MpiParameters params;
    process_command_line_options(argc, argv, params);
    mitm::PRNG prng(seed);
    if (rank == 0) {
        printf("************************************************************************\n");
        printf("double-speck64 demo! seed=%016" PRIx64 ", n=%d\n", prng.seed, n); 
    }

    params.setup(MPI_COMM_WORLD, 0);  // no controller process
    DoubleSpeck64_Problem Pb(n, prng);
    
#if 0
    if (params.verbose) {
        printf("==============================================================\n");
        printf("All-to-all version\n");
        if (expensive) printf("expensive f/g\n");
        printf("==============================================================\n");
    }
    std::vector<std::pair<u64, u64>> claws_alltoall, claws_isend;
    if (expensive)
        claws_alltoall = mitm::naive_mpi_claw_search_alltoall<true>(Pb, params);
    else
        claws_alltoall = mitm::naive_mpi_claw_search_alltoall<false>(Pb, params);
    if (params.verbose)
        for (auto it = claws_alltoall.begin(); it != claws_alltoall.end(); it++) {
            auto [x0, x1] = *it;
            assert(Pb.f(x0) == Pb.g(x1));
            printf("f(%" PRIx64 ") = g(%" PRIx64 ")\n", x0, x1);
        }
#endif
    if (params.verbose) {
        printf("==============================================================\n");
        printf("Isend version.\n");
        if (expensive) printf("expensive f/g\n");
        printf("==============================================================\n");
    }
    std::vector<std::pair<u64, u64>> claws_isend;
    if (expensive)
        claws_isend = mitm::naive_mpi_claw_search_isend<true>(Pb, params);
    else
        claws_isend = mitm::naive_mpi_claw_search_isend<false>(Pb, params);
    if (params.verbose)
        for (auto it = claws_isend.begin(); it != claws_isend.end(); it++) {
            auto [x0, x1] = *it;
            assert(Pb.f(x0) == Pb.g(x1));
            printf("f(%" PRIx64 ") = g(%" PRIx64 ")\n", x0, x1);
        }
    assert(claws_isend.size() == 1);
    MPI_Finalize();    
    return EXIT_SUCCESS;
}
