/// Primitive adaptation of speck blockcipher using 16-bit words.
/// speck32* : {0, 1}^32 x {0, 1}^32 -> {0, 1}^32, the official documentation
/// specifies different length for k, we are using shorter key for simplicity
/// since security is not a concern.

#include <stdint.h>
#include <iostream>
#include <array>


/* rot to right to the right by r amount */
inline uint16_t rot_r(uint16_t inp, uint16_t r) { return (inp>>r) | (inp<<(16-r)); }
/* rot to right to the left by r amount */
inline uint16_t rot_l(uint16_t inp, uint16_t r) { return (inp<<r) | (inp>>(16-r)); }

/* official parameters mentioned in speck paper */
constexpr int alpha = 7;
constexpr int beta = 2;
constexpr int nrounds = 22;


void round(std::array<uint16_t, 2>& text, uint16_t key) {
    /* Round function modifies the input one has to be careful */
    text[0]  = rot_r(text[0], alpha);
    text[0] += text[1] ;
    text[1] ^= key;
    text[1]  = rot_l(text[1], beta);
    text[1] ^= text[0];
}

void round_inv(std::array<uint16_t, 2>& text, uint16_t key) {
    /* Round function modifies the input one has to be careful */
    text[1] ^= text[0];
    text[1]  = rot_r(text[1], beta);
    text[1] ^= key;
    text[0] -= text[1];
    text[0]  = rot_l(text[0], alpha);
}


inline void next_round_key(std::array<uint16_t, 2> &keys, uint16_t i)
{
    keys[0]  = rot_r(keys[0], alpha);
    keys[0] += keys[1];
    keys[0] ^= i;
    keys[1]  = rot_l(keys[1], beta);
    keys[1] ^= keys[0];
}

void encrypt(std::array<uint16_t, 2>& text,
             std::array<uint16_t, 2>& cipher_text,
             std::array<uint16_t, 2>& key)
{
    cipher_text[0] = text[0];
    cipher_text[1] = text[1];
    std::array<uint16_t, 2> key_tmp;

    key_tmp[0] = key[0];
    key_tmp[1] = key[1];

    /* save the expanded key here */
    std::array<uint16_t, nrounds> keys{};
    keys[0] = key_tmp[1];

    /* key expansion */
    for (int i = 1; i < nrounds; ++i) {
        next_round_key(key_tmp, i);
        keys[i] = key_tmp[1];
    }

    /* actual encryption */
    for (int i = 0; i < nrounds; ++i) {
        round(cipher_text, keys[i]);
    }
}


void decrypt(std::array<uint16_t, 2>& cipher_text,
             std::array<uint16_t, 2>& plain_text,
             std::array<uint16_t, 2>& key)
{
    plain_text[0] = cipher_text[0];
    plain_text[1] = cipher_text[1];
    std::array<uint16_t, 2> key_tmp;

    key_tmp[0] = key[0];
    key_tmp[1] = key[1];

    std::array<uint16_t, nrounds> keys{};
    keys[0] = key_tmp[1];

    /* key expansion */
    for (int i = 1; i < nrounds; ++i) {
        next_round_key(key_tmp, i);
        keys[i] = key_tmp[1];
    }

    /* decryption */
    for (int i = 0; i < nrounds; ++i) {
        round_inv(plain_text, keys[nrounds - (i+1) ]);
    }
}


// -------------- END of SPECK32/32 specification ---------------
