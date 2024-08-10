#include <inttypes.h>
typedef uint16_t u16;

// derived from the SIMON and SPECK Implementation Guide

#define ROTL16(x,r) (((x)<<(r)) | (x>>(16-(r))))
#define ROTR16(x,r) (((x)>>(r)) | ((x)<<(16-(r))))

#define ER16(x,y,k) (x=ROTR16(x,7), x+=y, x^=k, y=ROTL16(y,2), y^=x)
#define DR16(x,y,k) (y^=x, y=ROTR16(y,2), x^=k, x-=y, x=ROTL16(x,7))

void Speck3264KeySchedule(u16 K[], u16 rk[])
{
	u16 i, D=K[3], C=K[2], B=K[1], A=K[0];
	for(i=0;i<22;){
		rk[i]=A; ER16(B,A,i++);
		rk[i]=A; ER16(C,A,i++);
		rk[i]=A; ER16(D,A,i++);
	}
}

void Speck3264Encrypt(u16 Pt[],u16 Ct[],u16 rk[])
{
	int i;
	Ct[0]=Pt[0]; Ct[1]=Pt[1];
	for(i=0;i<22;) 
		ER16(Ct[1],Ct[0],rk[i++]);
}

void Speck3264Decrypt(u16 Pt[],u16 Ct[],u16 rk[])
{
	int i;
	Pt[0]=Ct[0]; 
	Pt[1]=Ct[1];
	for(i=21;i>=0;) 
		DR16(Pt[1],Pt[0],rk[i--]);
}

#if 0
int main(void)
{
  uint16_t plaintext[2] = {0x694c, 0x6574};
  uint16_t key[4] = {0x0100, 0x0908, 0x1110, 0x1918};
  uint16_t ciphertext[2];
  u16 rk[22];
  Speck3264KeySchedule(key, rk);
  Speck3264Encrypt(plaintext, ciphertext, rk);
  printf("Plaintext:  0x%04x 0x%04x\n", plaintext[0], plaintext[1]);
  printf("Key:        0x%04x 0x%04x 0x%04x 0x%04x\n", key[0], key[1], key[2], key[3]);
  printf("Ciphertext: 0x%04x 0x%04x\n", ciphertext[0], ciphertext[1]);
  printf("Expected:   0x42f2 0xa868\n");
  
  plaintext[0] = plaintext[1] = 0;
  Speck3264Decrypt(plaintext, ciphertext, rk);
  printf("Plaintext:  0x%04x 0x%04x\n", plaintext[0], plaintext[1]);
  
  return 0;
}
#endif