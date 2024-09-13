#ifndef MITM_DESPROBLEM
#define MITM_DESPROBLEM

#include "problem.hpp"

extern "C"{
    /* We use the (sequential) DES implementation from OpenSSL */
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
/* but we ship our own bitsliced version */
void des_both(const u64 *enc_in_ortho, const u64 *dec_in_ortho, const u64 *keys, const u64 *enc, u64 *outputs);

namespace mitm {

////////////////////////////////////////////////////////////////////////////////
class DoubleDES_Problem : public mitm::AbstractClawProblem
{
public:
    const int n;
    const int m = 64;
    u64 in_mask;
    static constexpr int vlen = sizeof(v32) * 8;

    u64 P[2] = {0, 0xffffffffffffffffull};         /* two plaintext-ciphertext pairs */
    u64 C[2];

    v64 vP0[64], vC0[64];

    // 0 <= k < 2**56;
    static u64 des56(u64 k, u64 in, int enc)
    {
        u8 K[8] = {
            (u8) ((k >> 48) & 0xfe), (u8) ((k >> 41) & 0xfe), (u8) ((k >> 34) & 0xfe), (u8) ((k >> 27) & 0xfe), 
            (u8) ((k >> 20) & 0xfe), (u8) ((k >> 13) & 0xfe), (u8) ((k >>  6) & 0xfe), (u8) ((k <<  1) & 0xfe)
        };
        DES_key_schedule ks;
        DES_set_key_unchecked(K, &ks);
        u32 out[2];
        in = __builtin_bswap64(in);
        out[0] = (u32) (in & 0xffffffff);
        out[1] = (u32) (in >> 32);
        DES_encrypt1(out, &ks, enc);
        return __builtin_bswap64(out[0] ^ (((u64) out[1]) << 32));
    }

    // 0 <= k < 2**64 (with parity); big-endian order
    static u64 des64(u64 k, u64 in, int enc)
    {
        u8 K[8] = {
            (u8) ((k >> 56) & 0xff), (u8) ((k >> 48) & 0xff), (u8) ((k >> 40) & 0xff), (u8) ((k >> 32) & 0xff), 
            (u8) ((k >> 24) & 0xff), (u8) ((k >> 16) & 0xff), (u8) ((k >>  8) & 0xff), (u8) ((k      ) & 0xff)
        }; 
        DES_key_schedule ks;
        DES_set_key_unchecked(K, &ks);
        u32 out[2];
        in = __builtin_bswap64(in);
        out[0] = (u32) (in & 0xffffffff);
        out[1] = (u32) (in >> 32);
        DES_encrypt1(out, &ks, enc);
        return __builtin_bswap64(out[0] ^ (((u64) out[1]) << 32));
    }


    // encryption of P[0], using k
    u64 f(u64 k) const
    {
        assert((k & in_mask) == k);
        return des56(k, P[0], 1);
    }

    // decryption of C[0], using k
    u64 g(u64 k) const
    {
        assert((k & in_mask) == k);
        return des56(k, C[0], 0);
    }

    void vfg(const u64 k[], const bool choice[], u64 out[]) const
    {
        u64 enc[vlen / 64] __attribute__ ((aligned(sizeof(u64) * vlen)));
        for (int i = 0; i < vlen / 64; i++) {
            enc[i] = 0;
            for (int j = 0; j < 64; j++)
                enc[i] ^= ((u64) choice[i*64 + j]) << j;
        }

        des_both((u64 *) vP0, (u64 *) vC0, k, enc, out);
    }

    bool is_good_pair(u64 k0, u64 k1) const
    {
        u64 mid = des56(k0, P[1], 1);
        u64 out = des56(k1, mid, 1);
        return (out == C[1]);
    }

    static constexpr u64 test_inputs[64] = {
        0x95F8A5E5DD31D900ull, 0xDD7F121CA5015619ull, 0x2E8653104F3834EAull, 0x4BD388FF6CD81D4Full, 
        0x20B9E767B2FB1456ull, 0x55579380D77138EFull, 0x6CC5DEFAAF04512Full, 0x0D9F279BA5D87260ull, 
        0xD9031B0271BD5A0Aull, 0x424250B37C3DD951ull, 0xB8061B7ECD9A21E5ull, 0xF15D0F286B65BD28ull, 
        0xADD0CC8D6E5DEBA1ull, 0xE6D5F82752AD63D1ull, 0xECBFE3BD3F591A5Eull, 0xF356834379D165CDull, 
        0x2B9F982F20037FA9ull, 0x889DE068A16F0BE6ull, 0xE19E275D846A1298ull, 0x329A8ED523D71AECull, 
        0xE7FCE22557D23C97ull, 0x12A9F5817FF2D65Dull, 0xA484C3AD38DC9C19ull, 0xFBE00A8A1EF8AD72ull, 
        0x750D079407521363ull, 0x64FEED9C724C2FAFull, 0xF02B263B328E2B60ull, 0x9D64555A9A10B852ull, 
        0xD106FF0BED5255D7ull, 0xE1652C6B138C64A5ull, 0xE428581186EC8F46ull, 0xAEB5F5EDE22D1A36ull, 
        0xE943D7568AEC0C5Cull, 0xDF98C8276F54B04Bull, 0xB160E4680F6C696Full, 0xFA0752B07D9C4AB8ull, 
        0xCA3A2B036DBC8502ull, 0x5E0905517BB59BCFull, 0x814EEB3B91D90726ull, 0x4D49DB1532919C9Full, 
        0x25EB5FC3F8CF0621ull, 0xAB6A20C0620D1C6Full, 0x79E90DBC98F92CCAull, 0x866ECEDD8072BB0Eull, 
        0x8B54536F2F3E64A8ull, 0xEA51D3975595B86Bull, 0xCAFFC6AC4542DE31ull, 0x8DD45A2DDF90796Cull, 
        0x1029D55E880EC2D0ull, 0x5D86CB23639DBEA9ull, 0x1D1CA853AE7C0C5Full, 0xCE332329248F3228ull, 
        0x8405D1ABE24FB942ull, 0xE643D78090CA4207ull, 0x48221B9937748A23ull, 0xDD7C0BBD61FAFD54ull, 
        0x2FBC291A570DB5C4ull, 0xE07C30D7E4E26E12ull, 0x0953E2258E8E90A1ull, 0x5B711BC4CEEBF2EEull, 
        0xCC083F1E6D9E85F6ull, 0xD2FD8867D50D2DFEull, 0x06E7EA22CE92708Full, 0x166B40B44ABA4BD6ull
    }; 
    
    static constexpr u64 test_outputs[64] = {
        0x95F8A5E5DD31D900ull, 0xDD7F121CA5015619ull, 0x2E8653104F3834EAull, 0x4BD388FF6CD81D4Full,
        0x20B9E767B2FB1456ull, 0x55579380D77138EFull, 0x6CC5DEFAAF04512Full, 0x0D9F279BA5D87260ull,
        0xD9031B0271BD5A0Aull, 0x424250B37C3DD951ull, 0xB8061B7ECD9A21E5ull, 0xF15D0F286B65BD28ull,
        0xADD0CC8D6E5DEBA1ull, 0xE6D5F82752AD63D1ull, 0xECBFE3BD3F591A5Eull, 0xF356834379D165CDull,
        0x2B9F982F20037FA9ull, 0x889DE068A16F0BE6ull, 0xE19E275D846A1298ull, 0x329A8ED523D71AECull,
        0xE7FCE22557D23C97ull, 0x12A9F5817FF2D65Dull, 0xA484C3AD38DC9C19ull, 0xFBE00A8A1EF8AD72ull,
        0x750D079407521363ull, 0x64FEED9C724C2FAFull, 0xF02B263B328E2B60ull, 0x9D64555A9A10B852ull,
        0xD106FF0BED5255D7ull, 0xE1652C6B138C64A5ull, 0xE428581186EC8F46ull, 0xAEB5F5EDE22D1A36ull,
        0xE943D7568AEC0C5Cull, 0xDF98C8276F54B04Bull, 0xB160E4680F6C696Full, 0xFA0752B07D9C4AB8ull,
        0xCA3A2B036DBC8502ull, 0x5E0905517BB59BCFull, 0x814EEB3B91D90726ull, 0x4D49DB1532919C9Full,
        0x25EB5FC3F8CF0621ull, 0xAB6A20C0620D1C6Full, 0x79E90DBC98F92CCAull, 0x866ECEDD8072BB0Eull,
        0x8B54536F2F3E64A8ull, 0xEA51D3975595B86Bull, 0xCAFFC6AC4542DE31ull, 0x8DD45A2DDF90796Cull,
        0x1029D55E880EC2D0ull, 0x5D86CB23639DBEA9ull, 0x1D1CA853AE7C0C5Full, 0xCE332329248F3228ull,
        0x8405D1ABE24FB942ull, 0xE643D78090CA4207ull, 0x48221B9937748A23ull, 0xDD7C0BBD61FAFD54ull,
        0x2FBC291A570DB5C4ull, 0xE07C30D7E4E26E12ull, 0x0953E2258E8E90A1ull, 0x5B711BC4CEEBF2EEull,
        0xCC083F1E6D9E85F6ull, 0xD2FD8867D50D2DFEull, 0x06E7EA22CE92708Full, 0x166B40B44ABA4BD6ull
    };

    static constexpr u64 test_vectors[19][3] = {
        {0x7CA110454A1A6E57ull, 0x01A1D6D039776742ull, 0x690F5B0D9A26939Bull }, 
        {0x0131D9619DC1376Eull, 0x5CD54CA83DEF57DAull, 0x7A389D10354BD271ull }, 
        {0x07A1133E4A0B2686ull, 0x0248D43806F67172ull, 0x868EBB51CAB4599Aull }, 
        {0x3849674C2602319Eull, 0x51454B582DDF440Aull, 0x7178876E01F19B2Aull }, 
        {0x04B915BA43FEB5B6ull, 0x42FD443059577FA2ull, 0xAF37FB421F8C4095ull }, 
        {0x0113B970FD34F2CEull, 0x059B5E0851CF143Aull, 0x86A560F10EC6D85Bull }, 
        {0x0170F175468FB5E6ull, 0x0756D8E0774761D2ull, 0x0CD3DA020021DC09ull }, 
        {0x43297FAD38E373FEull, 0x762514B829BF486Aull, 0xEA676B2CB7DB2B7Aull }, 
        {0x07A7137045DA2A16ull, 0x3BDD119049372802ull, 0xDFD64A815CAF1A0Full }, 
        {0x04689104C2FD3B2Full, 0x26955F6835AF609Aull, 0x5C513C9C4886C088ull }, 
        {0x37D06BB516CB7546ull, 0x164D5E404F275232ull, 0x0A2AEEAE3FF4AB77ull }, 
        {0x1F08260D1AC2465Eull, 0x6B056E18759F5CCAull, 0xEF1BF03E5DFA575Aull }, 
        {0x584023641ABA6176ull, 0x004BD6EF09176062ull, 0x88BF0DB6D70DEE56ull }, 
        {0x025816164629B007ull, 0x480D39006EE762F2ull, 0xA1F9915541020B56ull }, 
        {0x49793EBC79B3258Full, 0x437540C8698F3CFAull, 0x6FBF1CAFCFFD0556ull }, 
        {0x4FB05E1515AB73A7ull, 0x072D43A077075292ull, 0x2F22E49BAB7CA1ACull }, 
        {0x49E95D6D4CA229BFull, 0x02FE55778117F12Aull, 0x5A6B612CC26CCE4Aull }, 
        {0x018310DC409B26D6ull, 0x1D9D5C5018F728C2ull, 0x5F4C038ED12B2E41ull }, 
        {0x1C587F1C13924FEFull, 0x305532286D6F295Aull, 0x63FAC0D034D9F793ull }
    };

    static u64 reverseBits(u64 x) {
        constexpr u64 M4_HI = 0xf0f0f0f0f0f0f0f0;
        constexpr u64 M4_LO = 0x0f0f0f0f0f0f0f0f;
        constexpr u64 M5_HI = 0xcccccccccccccccc;
        constexpr u64 M5_LO = 0x3333333333333333;
        constexpr u64 M6_HI = 0xaaaaaaaaaaaaaaaa;
        constexpr u64 M6_LO = 0x5555555555555555;

        x = (x & M4_HI) >> 4 | (x & M4_LO) << 4;  // Inverse les moitiés de 4 bits
        x = (x & M5_HI) >> 2 | (x & M5_LO) << 2;  // Inverse les paires de 2 bits
        x = (x & M6_HI) >> 1 | (x & M6_LO) << 1;  // Inverse les bits voisins
        return x;
    }

    static void openssl_des_test()
    {
        u8 K[8] = {0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01}; 
        DES_key_schedule ks;
        DES_set_key_unchecked(K, &ks);
        u64 in = __builtin_bswap64(0x95F8A5E5DD31D900ull);
        u32 out[2];
        out[0] = (u32) (in & 0xffffffff);
        out[1] = (u32) (in >> 32);
        DES_encrypt1(out, &ks, 1);
        u64 res = out[0] ^ (((u64) out[1]) << 32);
        assert(res == __builtin_bswap64(0x8000000000000000ull)); 

        u64 k = __builtin_bswap64(0x04B915BA43FEB5B6ull);
        K[0] = (u8) (k & 0xff);
        K[1] = (u8) ((k >> 8) & 0xff);
        K[2] = (u8) ((k >> 16) & 0xff);
        K[3] = (u8) ((k >> 24) & 0xff);
        K[4] = (u8) ((k >> 32) & 0xff);
        K[5] = (u8) ((k >> 40) & 0xff);
        K[6] = (u8) ((k >> 48) & 0xff);
        K[7] = (u8) ((k >> 56) & 0xff);

        // 0x02 0x5c 0x0a 0x5d 0x21 0x7f 0x5a 0x5b

        assert(K[0] == 0x04);
        assert(K[7] == 0xB6);
        DES_set_key_unchecked(K, &ks);

        in = __builtin_bswap64(0x42FD443059577FA2ull);
        out[0] = (u32) (in & 0xffffffff);
        out[1] = (u32) (in >> 32);
        DES_encrypt1(out, &ks, 1);
        res = out[0] ^ (((u64) out[1]) << 32);
        assert(res == __builtin_bswap64(0xAF37FB421F8C4095ull)); 

        assert(des64(0x04B915BA43FEB5B6ull, 0x42FD443059577FA2ull, 1) == 0xAF37FB421F8C4095ull);
       
        in = 0xdeadbeefbaadcafeull;
        assert(des64(0, des64(0, in, 1), 0) == in);

        // test des64 vs test vectors
        for (int j = 0; j < 64; j++) {
            u64 check = des64(0, test_inputs[j], 1);
            // printf("openssl: %016" PRIx64 "\n", check);
            assert(check == 1ull << (63-j));
        }
        
        for (int j = 0; j < 64; j++) {
            u64 check = des64(0, 1ull << (63-j), 1);
            // printf("openssl: %016" PRIx64 "\n", check);
            assert(check == test_outputs[j]);
        }

        for (int j = 0; j < 64; j++) {
            u64 check = des64(0, test_outputs[j], 0);
            // printf("openssl: %016" PRIx64 "\n", check);
            assert(check == 1ull << (63-j));
        }

        for (int j = 0; j < 19; j++) {
            auto [k, in, out] = test_vectors[j];
            u64 check = des64(k, in, 1);
            // printf("openssl: %016" PRIx64 " vs expected: %016" PRIx64 "\n", check, out);
            assert(check == out);
        }

        // test decryption
        for (int j = 0; j < 19; j++) {
            auto [k, in, out] = test_vectors[j];
            u64 check = des64(k, out, 0);
            // printf("openssl: %016" PRIx64 " vs expected: %016" PRIx64 "\n", check, in);
            assert(check == in);
        }

        // consistency between des56 and des64
        PRNG vprng;
        in = vprng.rand();
        u64 k64 = vprng.rand();
        u64 k56 = 0;
        for (int i = 0; i < 8; i++)
            k56 ^= (((k64 >> (1 + 8*i)) & 0x7f) << (7*i));
        assert(des64(k64, in, 1) == des56(k56, in, 1)); 
    }

    static void usuba_des_test()
    {
        PRNG vprng;
        u64 K[vlen] __attribute__ ((aligned(sizeof(u64) * vlen)));
        u64 out[vlen] __attribute__ ((aligned(sizeof(u64) * vlen)));

        for (int i = 0; i < vlen; i++)
            K[i] = vprng.rand() & 0x00ffffffffffffffull;
        
        u64 choice[vlen / 64] __attribute__ ((aligned(sizeof(u64) * vlen)));
        for (int i = 0; i < vlen / 64; i++)
            choice[i] = vprng.rand();

        v64 vin1[64], vin2[64];
        for (int i = 0; i < 64; i++) {
            vin1[i] = v64zero();
            vin2[i] = v64zero();
        }

        des_both((const u64*) vin1, (const u64*) vin2, K, choice, out);

        for (int i = 0; i < vlen; i++) {
            bool bit = (choice[i / 64] >> (i % 64)) & 1;
            u64 check = des56(K[i], 0, bit);
            // printf("in==0, usuba: %016" PRIx64 " vs openssl: %016" PRIx64 "\n", enc[i], check_enc);
            assert(out[i] == check);
        }

        u64 in1 = vprng.rand();
        u64 in2 = vprng.rand();
        for (int i = 0; i < 64; i++) {
            vin1[63-i] = (in1 >> i) & 1 ? v64bcast(-1) : v64bcast(0);
            vin2[63-i] = (in2 >> i) & 1 ? v64bcast(-1) : v64bcast(0);
        }
        des_both((const u64*) vin1, (const u64*) vin2, K, choice, out);
        for (int i = 0; i < vlen; i++) {
            bool bit = (choice[i / 64] >> (i % 64)) & 1;
            u64 check = bit ? des56(K[i], in1, 1) : des56(K[i], in2, 0);
            // printf("in==$$ usuba: %016" PRIx64 " vs openssl: %016" PRIx64 "\n", out[i], check);
            assert(out[i] == check);
        }
    }


    DoubleDES_Problem(int n, mitm::PRNG &prng) : n(n)
    {
        openssl_des_test();
        usuba_des_test();
        printf("DES: tests pass\n");

        in_mask = make_mask(n);
        assert(n <= 56);
        u64 k0 = prng.rand() & in_mask;
        u64 k1 = prng.rand() & in_mask;
        // printf("Secret keys = %016" PRIx64 " %016" PRIx64 "\n", k0, k1);
        u64 mid = des56(k0, P[0], 1);
        C[0] = des56(k1, mid, 1);
        mid = des56(k0, P[1], 1);
        C[1] = des56(k1, mid, 1);
        assert(f(k0) == g(k1));
        assert(is_good_pair(k0, k1));

        for (int i = 0; i < 64; i++) {
            vP0[63-i] = (P[0] >> i) & 1 ? v64bcast(-1) : v64bcast(0);
            vC0[63-i] = (C[0] >> i) & 1 ? v64bcast(-1) : v64bcast(0);
        }
    }
};
}
#endif