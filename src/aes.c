#include "aes.h"
//compile using gcc and following arguments: -g;-O0;-Wall;-msse2;-msse;-march=native;-maes

/*** MACROS ***/

#define ENC_BLOCK(m, k) \
    m = _mm_xor_si128(m, k[0x0]); \
    m = _mm_aesenc_si128(m, k[0x1]); \
    m = _mm_aesenc_si128(m, k[0x2]); \
    m = _mm_aesenc_si128(m, k[0x3]); \
    m = _mm_aesenc_si128(m, k[0x4]); \
    m = _mm_aesenc_si128(m, k[0x5]); \
    m = _mm_aesenc_si128(m, k[0x6]); \
    m = _mm_aesenc_si128(m, k[0x7]); \
    m = _mm_aesenc_si128(m, k[0x8]); \
    m = _mm_aesenc_si128(m, k[0x9]); \
    m = _mm_aesenclast_si128(m, k[0xa]);


#define DEC_BLOCK(m, k)\
    m = _mm_xor_si128(m, k[0xa + 0x0]); \
    m = _mm_aesdec_si128(m, k[0xa + 0x1]); \
    m = _mm_aesdec_si128(m, k[0xa + 0x2]); \
    m = _mm_aesdec_si128(m, k[0xa + 0x3]); \
    m = _mm_aesdec_si128(m, k[0xa + 0x4]); \
    m = _mm_aesdec_si128(m, k[0xa + 0x5]); \
    m = _mm_aesdec_si128(m, k[0xa + 0x6]); \
    m = _mm_aesdec_si128(m, k[0xa + 0x7]); \
    m = _mm_aesdec_si128(m, k[0xa + 0x8]); \
    m = _mm_aesdec_si128(m, k[0xa + 0x9]); \
    m = _mm_aesdeclast_si128(m, k[0x0]);

#define aes128_key_exp(k, rcon) aes128_key_expansion(k, _mm_aeskeygenassist_si128(k, rcon))


/*** PRIVATE ***/

static __m128i aes128_key_expansion(__m128i key, __m128i keygened)
{
    keygened = _mm_shuffle_epi32(keygened, _MM_SHUFFLE(0x3, 0x3, 0x3, 0x3));
    key = _mm_xor_si128(key, _mm_slli_si128(key, 0x4));
    key = _mm_xor_si128(key, _mm_slli_si128(key, 0x4));
    key = _mm_xor_si128(key, _mm_slli_si128(key, 0x4));
    return _mm_xor_si128(key, keygened);
}

static void aes128_loadkey_enc(uint8_t* enc_key, __m128i* key_schedule)
{
    key_schedule[0x00] = _mm_loadu_si128((const __m128i*) enc_key);
	key_schedule[0x01] = aes128_key_exp(key_schedule[0x0], 0x01);
	key_schedule[0x02] = aes128_key_exp(key_schedule[0x1], 0x02);
	key_schedule[0x03] = aes128_key_exp(key_schedule[0x2], 0x04);
	key_schedule[0x04] = aes128_key_exp(key_schedule[0x3], 0x08);
	key_schedule[0x05] = aes128_key_exp(key_schedule[0x4], 0x10);
	key_schedule[0x06] = aes128_key_exp(key_schedule[0x5], 0x20);
	key_schedule[0x07] = aes128_key_exp(key_schedule[0x6], 0x40);
	key_schedule[0x08] = aes128_key_exp(key_schedule[0x7], 0x80);
	key_schedule[0x09] = aes128_key_exp(key_schedule[0x8], 0x1b);
	key_schedule[0x0a] = aes128_key_exp(key_schedule[0x9], 0x36);
}

static void aes128_load_key(uint8_t* enc_key, __m128i* key_schedule)
{
    aes128_loadkey_enc(enc_key, key_schedule);
	key_schedule[0x0b] = _mm_aesimc_si128(key_schedule[0x9]);
	key_schedule[0x0c] = _mm_aesimc_si128(key_schedule[0x8]);
	key_schedule[0x0d] = _mm_aesimc_si128(key_schedule[0x7]);
	key_schedule[0x0e] = _mm_aesimc_si128(key_schedule[0x6]);
	key_schedule[0x0f] = _mm_aesimc_si128(key_schedule[0x5]);
	key_schedule[0x10] = _mm_aesimc_si128(key_schedule[0x4]);
	key_schedule[0x11] = _mm_aesimc_si128(key_schedule[0x3]);
	key_schedule[0x12] = _mm_aesimc_si128(key_schedule[0x2]);
	key_schedule[0x13] = _mm_aesimc_si128(key_schedule[0x1]);
}

static void aes128_enc(__m128i* key_schedule, uint8_t* plainText, uint8_t* cipherText)
{
    __m128i m = _mm_loadu_si128((__m128i*) plainText);
    ENC_BLOCK(m,key_schedule);
    _mm_storeu_si128((__m128i*) cipherText, m);
}


static void aes128_dec(__m128i* key_schedule, uint8_t* cipherText, uint8_t* plainText)
{
    __m128i m = _mm_loadu_si128((__m128i*) cipherText);
    DEC_BLOCK(m,key_schedule);
    _mm_storeu_si128((__m128i*) plainText, m);
}

/*** PUBLIC ***/

uint8_t* encrypt(uint8_t* key, uint8_t* plainText)
{
    __m128i key_schedule[0x14];
    aes128_loadkey_enc(key, key_schedule);
    uint8_t* cipherText;
    cipherText = (uint8_t*) malloc(sizeof(plainText));
    if(cipherText == NULL)
    {
        printf("ERROR !!!\n");
        exit(0x1);
    }
    aes128_enc(key_schedule, plainText, cipherText);
    return cipherText;
}

uint8_t* decrypt(uint8_t* key, uint8_t* cipherText)
{
    __m128i key_schedule[0x14];
    aes128_load_key(key, key_schedule);
    uint8_t* plainText;
    plainText = (uint8_t*) malloc(sizeof(cipherText));
    if(plainText == NULL)
    {
        printf("ERROR !!!\n");
        exit(0x1);
    }
    aes128_dec(key_schedule, cipherText, plainText);
    return plainText;
}

/* Return 0x0 if OK, 0x1 if encryption failed, 0x2 if decryption failed, 0x3 if both failed */
int self_test(void)
{
    uint8_t plain[0x10] = {0x32, 0x43, 0xf6, 0xa8, 0x88, 0x5a, 0x30, 0x8d, 0x31, 0x31, 0x98, 0xa2, 0xe0, 0x37, 0x07, 0x34};
    uint8_t enc_key[0x10] = {0x2b, 0x7e, 0x15, 0x16, 0x28, 0xae, 0xd2, 0xa6, 0xab, 0xf7, 0x15, 0x88, 0x09, 0xcf, 0x4f, 0x3c};
    uint8_t cipher[0x10] = {0x39, 0x25, 0x84, 0x1d, 0x02, 0xdc, 0x09, 0xfb, 0xdc, 0x11, 0x85, 0x97, 0x19, 0x6a, 0x0b, 0x32};
    uint8_t computed_cipher[0x10];
    uint8_t computed_plain[0x10];
    int out = 0x0;
    __m128i key_schedule[0x14];
    aes128_load_key(enc_key, key_schedule);
    aes128_enc(key_schedule, plain, computed_cipher);
    aes128_dec(key_schedule, cipher, computed_plain);
    if(memcmp(cipher, computed_cipher, sizeof(cipher)))
    {
        out = 0x1;
    }
    if(memcmp(plain, computed_plain, sizeof(plain)))
    {
        out |= 0x2;
    }
    return out;
}

