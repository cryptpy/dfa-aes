/* 
 *  AES lib using AES-NI instructions.
 */
 
#ifndef _AES_H
#define _AES_H

#include <stdio.h>
#include <stdint.h>     //for int8_t
#include <string.h>     //for memcmp
#include <wmmintrin.h>  //for intrinsics for AES-NI

uint8_t* encrypt(uint8_t *key, uint8_t* plainTex);
uint8_t* decrypt(uint8_t *key, uint8_t* cipherText);
int self_test(void);

#endif
