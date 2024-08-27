#ifndef MITM_SPECK64
#define MITM_SPECK64

#include "problem.hpp"

namespace mitm {

// derived from the SIMON and SPECK Implementation Guide

#define ROTL32(x,r) (((x)<<(r)) | (x>>(32-(r))))
#define ROTR32(x,r) (((x)>>(r)) | ((x)<<(32-(r))))

#define ER32(x,y,k) (x=ROTR32(x,8), x+=y, x^=k, y=ROTL32(y,3), y^=x)
#define DR32(x,y,k) (y^=x, y=ROTR32(y,3), x^=k, x-=y, x=ROTL32(x,8))

void Speck64128KeySchedule(const u32 K[],u32 rk[])
{
    u32 i,D=K[3],C=K[2],B=K[1],A=K[0];
    for(i=0;i<27;){
        rk[i]=A; ER32(B,A,i++);
        rk[i]=A; ER32(C,A,i++);
        rk[i]=A; ER32(D,A,i++);
    }
}

void Speck64128Encrypt(const u32 Pt[], u32 Ct[], const u32 rk[])
{
    u32 i;
    Ct[0]=Pt[0]; Ct[1]=Pt[1];
    for(i=0;i<27;)
        ER32(Ct[1],Ct[0],rk[i++]);
}

void Speck64128Decrypt(u32 Pt[], const u32 Ct[], u32 const rk[])
{
    int i;
    Pt[0]=Ct[0]; Pt[1]=Ct[1];
    for(i=26;i>=0;)
        DR32(Pt[1],Pt[0],rk[i--]);
}

// vectorized implementation using GCC vector extensions

void vSpeck64128KeySchedule(const v32 K[], v32 rk[])
{
    v32 D = K[3], C = K[2], B = K[1], A = K[0];
    for(int i = 0; i < 27; ) {
        rk[i] = A; ER32(B, A, i++);
        rk[i] = A; ER32(C, A, i++);
        rk[i] = A; ER32(D, A, i++);
    }
}

void vSpeck64128Encrypt(const v32 Pt[], v32 Ct[], const v32 rk[])
{
    Ct[0] = Pt[0]; Ct[1] = Pt[1];
    for(int i = 0; i < 27;)
        ER32(Ct[1], Ct[0], rk[i++]);
}

void vSpeck64128Decrypt(v32 Pt[], const v32 Ct[], const v32 rk[])
{
    Pt[0] = Ct[0]; Pt[1] = Ct[1];
    for(int i = 26; i >= 0;)
        DR32(Pt[1], Pt[0], rk[i--]);
}

// v32 vBcast(u32 x) {
//     constexpr int vlen = sizeof(v32) / sizeof(u32);
//     u32 vx[vlen]  __attribute__ ((aligned(sizeof(v32))));
//     for (int i = 0; i < vlen; i++)
//         vx[i] = x;
//     return *(v32 *) vx;
// }

////////////////////////////////////////////////////////////////////////////////
class DoubleSpeck64_Problem : public mitm::AbstractClawProblem
{
public:
    int n, m;
    u64 in_mask, out_mask;
    static constexpr int vlen = sizeof(v32) / sizeof(u32);

    u32 P[2][2] = {{0, 0}, {0xffffffff, 0xffffffff}};         /* two plaintext-ciphertext pairs */
    u32 C[2][2];

    // speck32-64 encryption of P[0], using k
    u64 f(u64 k) const
    {
        assert((k & in_mask) == k);
        u32 K[4] = {(u32) (k & 0xffffffff), (u32) ((k >> 32)), 0, 0};
        u32 rk[27];
        Speck64128KeySchedule(K, rk);
        u32 Ct[2];
        Speck64128Encrypt(P[0], Ct, rk);
        return ((u64) Ct[0] ^ ((u64) Ct[1] << 32)) & out_mask;
    }

    // speck32-64 decryption of C[0], using k
    u64 g(u64 k) const
    {
        assert((k & in_mask) == k);
        u32 K[4] = {(u32) (k & 0xffffffff), (u32) ((k >> 32)), 0, 0};
        u32 rk[27];
        Speck64128KeySchedule(K, rk);
        u32 Pt[2];
        Speck64128Decrypt(Pt, C[0], rk);
        return ((u64) Pt[0] ^ ((u64) Pt[1] << 32)) & out_mask;
    }

    void vfg(const u64 k[], u64 rf[], u64 rg[]) const
    {
        v32 zero = v32zero();
        v64 klo = v64load(k);
        v64 khi = v64load(&k[vlen / 2]);
        v32 K[4];
        v32desinterleave(klo, khi, &K[0], &K[1]);
        K[2] = zero;
        K[3] = zero;
        v32 rk[27];
        vSpeck64128KeySchedule(K, rk);
        
        v32 vP[2] = {zero, zero};
        v32 vMf[2];
        vSpeck64128Encrypt(vP, vMf, rk);
        v64 vmask = v64bcast(out_mask);
        v32interleave(vMf[0], vMf[1], vmask, (v64 *) &rf[0], (v64 *) &rf[vlen / 2]);
        
        v32 vC[2] = {v32bcast(C[0][0]), v32bcast(C[0][1])};
        v32 vMg[2];
        vSpeck64128Decrypt(vMg, vC, rk);
        v32interleave(vMg[0], vMg[1], vmask, (v64 *) &rg[0], (v64 *) &rg[vlen / 2]);
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

    DoubleSpeck64_Problem(int n, mitm::PRNG &prng) : n(n), m(n)
    {
        assert(n <= 64);
        in_mask = make_mask(n);
        out_mask = make_mask(m);
        u64 khi = prng.rand() & in_mask;
        u64 klo = prng.rand() & in_mask;
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

};

}
#endif