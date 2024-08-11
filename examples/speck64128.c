#include <inttypes.h>
typedef uint32_t u32;

// derived from the SIMON and SPECK Implementation Guide

#define ROTL32(x,r) (((x)<<(r)) | (x>>(32-(r))))
#define ROTR32(x,r) (((x)>>(r)) | ((x)<<(32-(r))))

#define ER32(x,y,k) (x=ROTR32(x,8), x+=y, x^=k, y=ROTL32(y,3), y^=x)
#define DR32(x,y,k) (y^=x, y=ROTR32(y,3), x^=k, x-=y, x=ROTL32(x,8))

void Speck64128KeySchedule(u32 K[],u32 rk[])
{
    u32 i,D=K[3],C=K[2],B=K[1],A=K[0];
    for(i=0;i<27;){
        rk[i]=A; ER32(B,A,i++);
        rk[i]=A; ER32(C,A,i++);
        rk[i]=A; ER32(D,A,i++);
    }
}

void Speck64128Encrypt(u32 Pt[],u32 Ct[],u32 rk[])
{
    u32 i;
    Ct[0]=Pt[0]; Ct[1]=Pt[1];
    for(i=0;i<27;)
        ER32(Ct[1],Ct[0],rk[i++]);
}

void Speck64128Decrypt(u32 Pt[],u32 Ct[],u32 rk[])
{
    int i;
    Pt[0]=Ct[0]; Pt[1]=Ct[1];
    for(i=26;i>=0;)
        DR32(Pt[1],Pt[0],rk[i--]);
}

#if 0
#include <stdio.h>
int main(void)
{
  u32 plaintext[2] = {0x7475432d, 0x3b726574};
  u32 key[4] = {0x03020100, 0x0b0a0908, 0x13121110, 0x1b1a1918,};
  u32 ciphertext[2];
  u32 rk[27];
  Speck64128KeySchedule(key, rk);
  Speck64128Encrypt(plaintext, ciphertext, rk);
  printf("Plaintext:  0x%08x 0x%08x\n", plaintext[0], plaintext[1]);
  printf("Key:        0x%08x 0x%08x 0x%08x 0x%08x\n", key[0], key[1], key[2], key[3]);
  printf("Ciphertext: 0x%08x 0x%08x\n", ciphertext[0], ciphertext[1]);
  printf("Expected:   0x454e028b 0x8c6fa548\n");
  
  plaintext[0] = plaintext[1] = 0;
  Speck64128Decrypt(plaintext, ciphertext, rk);
  printf("Plaintext:  0x%08x 0x%08x\n", plaintext[0], plaintext[1]);
  
  return 0;
}
#endif