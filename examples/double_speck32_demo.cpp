#include <cassert>
#include <getopt.h>
#include <err.h>

#include "mitm.hpp"

/* We would like to call C function defined in `sha256.c` */
extern "C"{
    void Speck3264KeySchedule(const u16 K[], u16 rk[]);
    void Speck3264Encrypt(const u16 Pt[],u16 Ct[], const u16 rk[]);
    void Speck3264Decrypt(u16 Pt[],const u16 Ct[],const u16 rk[]);
}


int n = 20;         // default problem size (easy)
u64 seed = 0x1337;  // default fixed seed

////////////////////////////////////////////////////////////////////////////////
class DoubleSpeck32_Problem : mitm::AbstractClawProblem
{
public:
    int n;
    u64 mask;
    mitm::PRNG &prng;

    u16 P[2][2] = {{0, 0}, {0xffff, 0xffff}};         /* two plaintext-ciphertext pairs */
    u16 C[2][2];
    
    // speck32-64 encryption of P[0], using k
    u64 f(u64 k) const
    {
        assert((k & mask) == k);
        u16 K[4] = {(u16) (k & 0xffff), (u16) ((k >> 16) & 0xffff), 0, 0};
        u16 rk[22];
        Speck3264KeySchedule(K, rk);
        u16 Ct[2];
        Speck3264Encrypt(P[0], Ct, rk);
        return ((u64) Ct[0] ^ ((u64) Ct[1] << 16)) & mask;
    }

    // speck32-64 decryption of C[0], using k
    u64 g(u64 k) const
    {
        assert((k & mask) == k);
        u16 K[4] = {(u16) (k & 0xffff), (u16) ((k >> 16) & 0xffff), 0, 0};
        u16 rk[22];
        Speck3264KeySchedule(K, rk);
        u16 Pt[2];
        Speck3264Decrypt(Pt, C[0], rk);
        return ((u64) Pt[0] ^ ((u64) Pt[1] << 16)) & mask;
    }

  DoubleSpeck32_Problem(int n, mitm::PRNG &prng) : n(n), prng(prng)
  {
    assert(n <= 32);
    mask = (1ull << n) - 1;
    u64 khi = prng.rand() & mask;
    u64 klo = prng.rand() & mask;
    printf("Secret keys = %08" PRIx64 " %08" PRIx64 "\n", khi, klo);
    u16 Ka[4] = {(u16) (khi & 0xffff), (u16) ((khi >> 16) & 0xffff), 0, 0};
    u16 Kb[4] = {(u16) (klo & 0xffff), (u16) ((klo >> 16) & 0xffff), 0, 0};
    u16 rka[22];
    u16 rkb[22];
    Speck3264KeySchedule(Ka, rka);
    Speck3264KeySchedule(Kb, rkb);
    u16 mid[2][2];
    Speck3264Encrypt(P[0], mid[0], rka);
    Speck3264Encrypt(mid[0], C[0], rkb);
    Speck3264Encrypt(P[1], mid[1], rka);
    Speck3264Encrypt(mid[1], C[1], rkb);
  
    assert(f(khi) == g(klo));
    assert(is_good_pair(khi, klo));
  }

  bool is_good_pair(u64 khi, u64 klo) const 
  {
    u16 Ka[4] = {(u16) (khi & 0xffff), (u16) ((khi >> 16) & 0xffff), 0, 0};
    u16 Kb[4] = {(u16) (klo & 0xffff), (u16) ((klo >> 16) & 0xffff), 0, 0};
    u16 rka[22];
    u16 rkb[22];
    Speck3264KeySchedule(Ka, rka);
    Speck3264KeySchedule(Kb, rkb);
    u16 mid[2];
    u16 Ct[2];
    Speck3264Encrypt(P[1], mid, rka);
    Speck3264Encrypt(mid, Ct, rkb);
    return (Ct[0] == C[1][0]) && (Ct[1] == C[1][1]);
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
        printf("double-speck32 demo! seed=%016" PRIx64 ", n=%d\n", prng.seed, n); 

        DoubleSpeck32_Problem Pb(n, prng);            
        auto claw = mitm::claw_search<mitm::OpenMPEngine>(Pb, params, prng);
        printf("f(%" PRIx64 ") = g(%" PRIx64 ")\n", claw.first, claw.second);
        
        return EXIT_SUCCESS;
}
