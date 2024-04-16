#ifndef MITM_BITS_LIB
#define MITM_BITS_LIB
/* Thank you ChatGPT4 */


#include <cstdint>
#include <cstring>

inline
void bits_memcpy(void *dst,  void const* const src, size_t nbits) {
  size_t nbytes = nbits / 8;
  size_t rem_bits = nbits % 8;

  uint8_t *src_u8 = (uint8_t*) src;
  uint8_t *dst_u8 = (uint8_t*) dst;
  

  // Copy full bytes
  std::memcpy(dst_u8, src_u8, nbytes);

  // Handle the remaining bits, if any
  if (rem_bits > 0) {
    uint8_t mask = (1 << rem_bits) - 1;
    dst_u8[nbytes] &= ~mask;
    dst_u8[nbytes] |= src_u8[nbytes] & mask;
  }
}


inline
int bits_memcmp(void const* const dst, void const* const src, size_t nbits) {
  size_t nbytes = nbits / 8;
  size_t rem_bits = nbits % 8;

  uint8_t *src_u8 = (uint8_t*) src;
  uint8_t *dst_u8 = (uint8_t*) dst;


  // Compare full bytes
  int result = std::memcmp(dst_u8, src_u8, nbytes);
  if (result != 0 || rem_bits == 0) {
    return result;
  }

  // Compare the remaining bits
  uint8_t mask = (1 << rem_bits) - 1;
  uint8_t dst_last_bits = dst_u8[nbytes] & mask;
  uint8_t src_last_bits = src_u8[nbytes] & mask;

  if (dst_last_bits < src_last_bits) {
    return -1;
  } else if (dst_last_bits > src_last_bits) {
    return 1;
  }

  return 0;
}


#endif
