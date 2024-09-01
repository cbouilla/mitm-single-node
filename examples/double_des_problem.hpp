#ifndef MITM_DESPROBLEM
#define MITM_DESPROBLEM

#include "problem.hpp"

/* We use the DES implementation from OpenSSL --- but we ship our own vectorized one... */
extern "C"{
    typedef struct DES_ks {
        union {
            u8 cblock[8];
            /*
             * make sure things are correct size on machines with 8 byte longs
             */
            u32 deslong[2];
        } ks[16];
    } DES_key_schedule;

    void DES_set_key_unchecked(const u8 *key, DES_key_schedule *schedule);
    void DES_encrypt1(u32 *data, DES_key_schedule *ks, int enc);
}

namespace mitm {

////////////////////////////////////////////////////////////////////////////////
class DoubleDES_Problem : public mitm::AbstractClawProblem
{
public:
    const int n;
    const int m = 64;
    const u64 in_mask;
    static constexpr int vlen = sizeof(v32) * 8;

    u32 P[2][2] = {{0, 0}, {0xffffffff, 0xffffffff}};         /* two plaintext-ciphertext pairs */
    u32 C[2][2];

    void des(u64 k, const u32 *in, u32 *out, int enc) const
    {
        u8 K[8] = {(u8) ((k <<  1) & 0xfe), (u8) ((k >>  6) & 0xfe), (u8) ((k >> 13) & 0xfe), (u8) ((k >> 20) & 0xfe),
                   (u8) ((k >> 27) & 0xfe), (u8) ((k >> 34) & 0xfe), (u8) ((k >> 41) & 0xfe), (u8) ((k >> 48) & 0xfe)}; 
        DES_key_schedule ks;
        DES_set_key_unchecked(K, &ks);
        out[0] = in[0];
        out[1] = in[1];
        DES_encrypt1(out, &ks, enc);
    }

    // encryption of P[0], using k
    u64 f(u64 k) const
    {
        assert((k & in_mask) == k);
        u32 cc[2];
        des(k, P[0], cc, 1);
        return cc[0] ^ (((u64) cc[1]) << 32);
    }

    // decryption of C[0], using k
    u64 g(u64 k) const
    {
        assert((k & in_mask) == k);
        u32 pp[2];
        des(k, C[0], pp, 0);
        return pp[0] ^ (((u64) pp[1]) << 32);
    }

/*
    void vfg(const u64 k[], u64 rf[], u64 rg[]) const
    {
        v64 key_ortho[KEY_SIZE];
        v64 plain_ortho[vlen];
        v64 cipher_ortho[vlen];
  
        for (int i = 0; i < vlen; i++)
            plain_std[i] = __builtin_bswap64(plain_std[i]);

        orthogonalize(plain_std, plain_ortho);
        des__(plain_ortho, key_ortho, cipher_ortho);
        unorthogonalize(cipher_ortho,plain_std);
    
        for (int i = 0; i < vlen; i++)
            plain_std[i] = __builtin_bswap64(plain_std[i]);

    }
*/
    bool is_good_pair(u64 k0, u64 k1) const
    {
        u32 mid[2], out[2];
        des(k0, P[1], mid, 1);
        des(k1, mid, out, 1);
        return (out[0] == C[1][0]) && (out[1] == C[1][1]);
    }

  DoubleDES_Problem(int n, mitm::PRNG &prng) : n(n), in_mask(make_mask(n))
  {
    assert(n <= 56);
    u64 k0 = prng.rand() & in_mask;
    u64 k1 = prng.rand() & in_mask;
    printf("Secret keys = %016" PRIx64 " %016" PRIx64 "\n", k0, k1);
    u32 mid[2];
    des(k0, P[0], mid, 1);
    des(k1, mid, C[0], 1);
    des(k0, P[1], mid, 1);
    des(k1, mid, C[1], 1);
    assert(f(k0) == g(k1));
    assert(is_good_pair(k0, k1));
  }
};

}
#endif