/**
 *  Licensed by "The MIT License". See file LICENSE.
 */

#ifndef DFA_H
#define DFA_H

#include <algorithm>
#include <array>
#include <climits>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <omp.h>
#include <set>
#include <sstream>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <string.h>
#include <vector>

#include "constant.hpp"

using namespace std;

/* 4-byte array */
using KeyTuple = array<uint8_t, 0x4>;

/* 16-byte state */
using State = array<uint8_t, 0x10>;

/* Data structure for differentials */
using DiffStat = array<multimap<uint8_t,uint8_t>, 0x10>;

/* Data structure for vector of key candidate tuples */
using VKeyTuple = vector<KeyTuple>;

/* Galois multiplication tables */
static map<uint8_t, const uint8_t*> gmt =
{
    {0x01, gm_01},
    {0x09, gm_09},
    {0x0b, gm_0b},
    {0x0d, gm_0d},
    {0x0e, gm_0e},
    {0x8d, gm_8d},
    {0xf6, gm_f6}
};

/* Related bytes (column-wise) */
static const uint8_t rb[0x4][0x4] =
{
    {0x0, 0x7, 0xa, 0xd},
    {0x1, 0x4, 0xb, 0xe},
    {0x2, 0x5, 0x8, 0xf},
    {0x3, 0x6, 0x9, 0xc}
};

/* Maps a fault location 'l' to the correct set of fault deltas for the standard filter. Note: 'l' is enumerated column-wise */
static const size_t map_fault[0x10] = {0x0, 0x1, 0x2, 0x3, 0x3, 0x0, 0x1, 0x2, 0x2, 0x3, 0x0, 0x1, 0x1, 0x2, 0x3, 0x0};

/* Pointer to inverses of fault deltas in GF(256) for the standard filter (depend on the fault location) */
static const uint8_t* ideltas1[0x4][0x10] = 
{
    {gm_8d, gm_01, gm_8d, gm_01, gm_01, gm_f6, gm_01, gm_f6, gm_01, gm_8d, gm_01, gm_8d, gm_f6, gm_01, gm_f6, gm_01},
    {gm_01, gm_f6, gm_01, gm_f6, gm_01, gm_8d, gm_01, gm_8d, gm_f6, gm_01, gm_f6, gm_01, gm_8d, gm_01, gm_8d, gm_01},
    {gm_01, gm_8d, gm_01, gm_8d, gm_f6, gm_01, gm_f6, gm_01, gm_8d, gm_01, gm_8d, gm_01, gm_01, gm_f6, gm_01, gm_f6},
    {gm_f6, gm_01, gm_f6, gm_01, gm_8d, gm_01, gm_8d, gm_01, gm_01, gm_f6, gm_01, gm_f6, gm_01, gm_8d, gm_01, gm_8d}
};

/* Pointer to inverses of fault deltas in GF(256) for the improved filter (depend on the fault location) */
static const uint8_t* ideltas2[0x4][0x4] =
{
    {gm_8d, gm_01, gm_01, gm_f6},
    {gm_f6, gm_8d, gm_01, gm_01},
    {gm_01, gm_f6, gm_8d, gm_01},
    {gm_01, gm_01, gm_f6, gm_8d}
};

/* indices for c,d and the 0xa-th round key k (improved fault equations) */
static const size_t indices_x[0x4][0x10] =
{
    {0x0, 0xd, 0xa, 0x7, 0xc, 0x9, 0x6, 0x3, 0x8, 0x5, 0x2, 0xf, 0x4, 0x1, 0xe, 0xb},
    {0xc, 0x9, 0x6, 0x3, 0x8, 0x5, 0x2, 0xf, 0x4, 0x1, 0xe, 0xb, 0x0, 0xd, 0xa, 0x7},
    {0x8, 0x5, 0x2, 0xf, 0x4, 0x1, 0xe, 0xb, 0x0, 0xd, 0xa, 0x7, 0xc, 0x9, 0x6, 0x3},
    {0x4, 0x1, 0xe, 0xb, 0x0, 0xd, 0xa, 0x7, 0xc, 0x9, 0x6, 0x3, 0x8, 0x5, 0x2, 0xf}
};

/* Indices for the 0x9-th round key h (improved fault equations)*/
static const size_t indices_y[0x4][16] =
{
    {0x0, 0x1, 0x2, 0x3, 0xc, 0xd, 0xe, 0xf, 0x8, 0x9, 0xa, 0xb, 0x4, 0x5, 0x6, 0x7},
    {0xc, 0xd, 0xe, 0xf, 0x8, 0x9, 0xa, 0xb, 0x4, 0x5, 0x6, 0x7, 0x0, 0x1, 0x2, 0x3},
    {0x8, 0x9, 0xa, 0xb, 0x4, 0x5, 0x6, 0x7, 0x0, 0x1, 0x2, 0x3, 0xc, 0xd, 0xe, 0xf},
    {0x4, 0x5, 0x6, 0x7, 0x0, 0x1, 0x2, 0x3, 0xc, 0xd, 0xe, 0xf, 0x8, 0x9, 0xa, 0xb}
};

vector<State> analyse(State &c, State &d, const size_t l, const size_t cores);

DiffStat differentials(State &c, State &d, const size_t l);

DiffStat standard_filter(DiffStat x);

vector<VKeyTuple> combine(DiffStat x);

vector<vector<VKeyTuple>> preproc(vector<VKeyTuple> cmb, const size_t cores);

vector<State> improved_filter(State &c, State &d, vector<VKeyTuple> &v, const size_t l);

vector<State> postproc(vector<vector<State>> &v);

State reconstruct(State &k);

uint32_t ks_core(uint32_t t, size_t r);

vector<uint8_t> getKeys(multimap<uint8_t, uint8_t> m);

void printState(State x);

vector<pair<pair<State, State>, State>> readfile(const string file, int bf);

void writefile(State plaintext, State ciphertext, vector<State> keys, const string file);

static uint8_t EQ(const uint8_t c, const uint8_t d, const uint8_t k, const uint8_t* gm)
{
    return gm[isbox[c ^ k] ^ isbox[d ^ k]];
}

template<typename T>
static constexpr size_t bits(T = T{})
{
    return sizeof(T) * CHAR_BIT;
}

template<typename T>
static constexpr T ROTL(const T x, const unsigned c)
{
    return (x << c) | (x >> (bits(x) - c));
}

template<typename T>
static constexpr T ROTR(const T x, const unsigned c)
{
    return (x >> c) | (x << (bits(x) - c));
}

#endif
