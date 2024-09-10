/* source: https://github.com/mrdcvlsc/AES */
#include <emmintrin.h>
#include <immintrin.h>
#include <xmmintrin.h>

static const size_t Nr = 10;


static inline __m128i AES_128_ASSIST(__m128i tmp1, __m128i tmp2)
{
	__m128i tmp3;
	tmp2 = _mm_shuffle_epi32(tmp2, 0xff);
	tmp3 = _mm_slli_si128(tmp1, 0x4);
	tmp1 = _mm_xor_si128(tmp1, tmp3);
	tmp3 = _mm_slli_si128(tmp3, 0x4);
	tmp1 = _mm_xor_si128(tmp1, tmp3);
	tmp3 = _mm_slli_si128(tmp3, 0x4);
	tmp1 = _mm_xor_si128(tmp1, tmp3);
	tmp1 = _mm_xor_si128(tmp1, tmp2);
	return tmp1;
}

void aes128_key_expansion(const unsigned char *user_key, unsigned char *key)
{
	__m128i tmp1, tmp2;
	__m128i *key_sched = (__m128i *) key;
	tmp1 = _mm_loadu_si128((__m128i *) user_key);
	key_sched[0] = tmp1;
	tmp2 = _mm_aeskeygenassist_si128(tmp1, 0x1);
	tmp1 = AES_128_ASSIST(tmp1, tmp2);
	key_sched[1] = tmp1;
	tmp2 = _mm_aeskeygenassist_si128(tmp1, 0x2);
	tmp1 = AES_128_ASSIST(tmp1, tmp2);
	key_sched[2] = tmp1;
	tmp2 = _mm_aeskeygenassist_si128(tmp1, 0x4);
	tmp1 = AES_128_ASSIST(tmp1, tmp2);
	key_sched[3] = tmp1;
	tmp2 = _mm_aeskeygenassist_si128(tmp1, 0x8);
	tmp1 = AES_128_ASSIST(tmp1, tmp2);
	key_sched[4] = tmp1;
	tmp2 = _mm_aeskeygenassist_si128(tmp1, 0x10);
	tmp1 = AES_128_ASSIST(tmp1, tmp2);
	key_sched[5] = tmp1;
	tmp2 = _mm_aeskeygenassist_si128(tmp1, 0x20);
	tmp1 = AES_128_ASSIST(tmp1, tmp2);
	key_sched[6] = tmp1;
	tmp2 = _mm_aeskeygenassist_si128(tmp1, 0x40);
	tmp1 = AES_128_ASSIST(tmp1, tmp2);
	key_sched[7] = tmp1;
	tmp2 = _mm_aeskeygenassist_si128(tmp1, 0x80);
	tmp1 = AES_128_ASSIST(tmp1, tmp2);
	key_sched[8] = tmp1;
	tmp2 = _mm_aeskeygenassist_si128(tmp1, 0x1b);
	tmp1 = AES_128_ASSIST(tmp1, tmp2);
	key_sched[9] = tmp1;
	tmp2 = _mm_aeskeygenassist_si128(tmp1, 0x36);
	tmp1 = AES_128_ASSIST(tmp1, tmp2);
	key_sched[10] = tmp1;
}


// Performs AES encryption to a 16 byte block of memory.
void aes128_encrypt(const unsigned char *key, const unsigned char *plaintext, unsigned char *ciphertext) 
{
	unsigned char round_keys[176];
	aes128_key_expansion(key, round_keys);
	
	  // load the current block & current round key into the registers
	  __m128i *xmm_round_keys = (__m128i *) round_keys;
	  __m128i state = _mm_loadu_si128((__m128i *) &plaintext[0]);

	  // original key
	  state = _mm_xor_si128(state, xmm_round_keys[0]);

	  // perform usual rounds
	  for (size_t i = 1; i < Nr - 1; i += 2) {
		state = _mm_aesenc_si128(state, xmm_round_keys[i]);
		state = _mm_aesenc_si128(state, xmm_round_keys[i + 1]);
	  }

	  // last round
	  state = _mm_aesenc_si128(state, xmm_round_keys[Nr - 1]);
	  state = _mm_aesenclast_si128(state, xmm_round_keys[Nr]);

	  // store from register to array
	  _mm_storeu_si128((__m128i *) (ciphertext), state);
}

// AES decryption to a 16 byte block of memory.
void aes128_decrypt(const unsigned char *key, const unsigned char *ciphertext, unsigned char *plaintext)
{
	unsigned char round_keys[176];
	aes128_key_expansion(key, round_keys);

	  // load the current block & current round key into the registers
	  __m128i *xmm_round_keys = (__m128i *) round_keys;
	  __m128i state = _mm_loadu_si128((__m128i *) &ciphertext[0]);

	  // first round
	  state = _mm_xor_si128(state, xmm_round_keys[Nr]);

	  // usual rounds
	  for (size_t i = Nr - 1; i > 1; i -= 2) {
		state = _mm_aesdec_si128(state, _mm_aesimc_si128(xmm_round_keys[i]));
		state = _mm_aesdec_si128(state, _mm_aesimc_si128(xmm_round_keys[i - 1]));
	  }

	  // last round
	  state = _mm_aesdec_si128(state, _mm_aesimc_si128(xmm_round_keys[1]));
	  state = _mm_aesdeclast_si128(state, xmm_round_keys[0]);

	  // store from register to array
	  _mm_storeu_si128((__m128i *) plaintext, state);
}
