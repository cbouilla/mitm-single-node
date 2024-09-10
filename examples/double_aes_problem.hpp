#ifndef MITM_SPECK64
#define MITM_SPECK64

#include "problem.hpp"


/* We would like to call C function defined in `sha256.c` */
extern "C"{
    void aes128_key_expansion(const u64 *user_key, u64 *key);
    void aes128_encrypt(const u64 *key, const u64 *plaintext, u64 *ciphertext) ;
    void aes128_decrypt(const u64 *key, const u64 *ciphertext, u64 *plaintext);
}


namespace mitm {

////////////////////////////////////////////////////////////////////////////////
class DoubleAES_Problem : public mitm::AbstractClawProblem
{
public:
    int n, m;
    u64 mask;

    u64 P[2][2] = {{0, 0}, {0xffffffffffffffffull, 0xffffffffffffffffull}};         /* two plaintext-ciphertext pairs */
    u64 C[2][2];

    // encryption of P[0], using k
    u64 f(u64 k) const
    {
        assert((k & mask) == k);
        u64 K[2] = {k, 0};
        u64 rk[22];
        aes128_key_expansion(K, rk);
        u64 Ct[2];
        aes128_encrypt(rk, P[0], Ct);
        return Ct[0] & mask;
    }

    // decryption of C[0], using k
    u64 g(u64 k) const
    {
        assert((k & mask) == k);
        u64 K[2] = {k, 0};
        u64 rk[22];
        aes128_key_expansion(K, rk);
        u64 Pt[2];
        aes128_decrypt(rk, C[0], Pt);
        return Pt[0] & mask;
    }


    bool is_good_pair(u64 khi, u64 klo) const
    {
        u64 Ka[4] = {khi, 0};
        u64 Kb[4] = {klo, 0};
        u64 rka[22];
        u64 rkb[22];
        aes128_key_expansion(Ka, rka);
        aes128_key_expansion(Kb, rkb);
        u64 mid[2];
        u64 Ct[2];
        aes128_encrypt(rka, P[1], mid);
        aes128_encrypt(rkb, mid, Ct);
        return (Ct[0] == C[1][0]) && (Ct[1] == C[1][1]);
    }

  DoubleAES_Problem(int n, mitm::PRNG &prng) : n(n), m(n)
  {
    assert(n <= 64);
    mask = (1ull << n) - 1;
    u64 khi = prng.rand() & mask;
    u64 klo = prng.rand() & mask;
    // printf("Secret keys = %016" PRIx64 " %016" PRIx64 "\n", khi, klo);
    u64 Ka[2] = {khi, 0};
    u64 Kb[2] = {klo, 0};
    u64 rka[22];
    u64 rkb[22];
    aes128_key_expansion(Ka, rka);
    aes128_key_expansion(Kb, rkb);
    u64 mid[2][2];
    aes128_encrypt(rka, P[0], mid[0]);
    aes128_encrypt(rka, P[1], mid[1]);
    aes128_encrypt(rkb, mid[0], C[0]);
    aes128_encrypt(rkb, mid[1], C[1]);

    assert(f(khi) == g(klo));
    assert(is_good_pair(khi, klo));
  }
};

}
#endif