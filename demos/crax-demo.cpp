#ifdef __cplusplus
extern "C" {
#endif
#define ROT(x, n) (((x) >> (n)) | ((x) << (32-(n))))

#define ALZETTE(x, y, c)                    \
  (x) += ROT((y), 31), (y) ^= ROT((x), 24), \
  (x) ^= (c),                               \
  (x) += ROT((y), 17), (y) ^= ROT((x), 17), \
  (x) ^= (c),                               \
  (x) += (y),          (y) ^= ROT((x), 31), \
  (x) ^= (c),                               \
  (x) += ROT((y), 24), (y) ^= ROT((x), 16), \
  (x) ^= (c)

#define ALZETTE_INV(x, y, c)                \
  (x) ^= (c),                               \
  (y) ^= ROT((x), 16), (x) -= ROT((y), 24), \
  (x) ^= (c),                               \
  (y) ^= ROT((x), 31), (x) -= (y),          \
  (x) ^= (c),                               \
  (y) ^= ROT((x), 17), (x) -= ROT((y), 17), \
  (x) ^= (c),                               \
  (y) ^= ROT((x), 24), (x) -= ROT((y), 31)

#define N_STEPS 10

static const uint32_t RCON[5] = {                             \
  0xB7E15162, 0xBF715880, 0x38B4DA56, 0x324E7738, 0xBB1185EB };

void craxs10_enc_ref(uint32_t *xword, uint32_t *yword, const uint32_t *key)
{
    int step;

    for (step = 0; step < NSTEPS; step++) {
        xword[0] ^= step;
        if ((step % 2) == 0) {
            xword[0] ^= key[0];
            yword[0] ^= key[1];
        } else {
            xword[0] ^= key[2];
            yword[0] ^= key[3];
        }
        ALZETTE(xword[0], yword[0], RCON[step%5]);
    }
    xword[0] ^= key[0];
    yword[0] ^= key[1];
}

void craxs10_dec_ref(uint32_t *xword, uint32_t *yword, const uint32_t *key)
{
    int step;

    xword[0] ^= key[0];
    yword[0] ^= key[1];
    for (step = NSTEPS-1; step >= 0; step--) {
        ALZETTE_INV(xword[0], yword[0], RCON[step%5]);
        if ((step % 2) == 0) {
            xword[0] ^= key[0];
            yword[0] ^= key[1];
        } else {
            xword[0] ^= key[2];
            yword[0] ^= key[3];
        }
        xword[0] ^= step;
    }
}
#ifdef __cplusplus
}
#endif