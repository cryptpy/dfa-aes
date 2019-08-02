/**
 *  Licensed by "The MIT License". See file LICENSE.
 */

#include "aes.h"
#include "dfa.hpp"

/* Start of differential fault analysis */
vector<State> analyse(State &c, State &d, const size_t l, const size_t cores)
{
    printf("Applying standard filter.");
    vector<VKeyTuple> cmb = combine(standard_filter(differentials(c, d, l)));
    printf("Done.\n");
    size_t n = cmb[0x0].size() * cmb[0x1].size() * cmb[0x2].size() * cmb[0x3].size();
    printf("Size of keyspace: %lu = 2^%f \n", n, log2(n));

    printf("Applying improved filter.");
    fflush(stdout);

    /* Pre-processing */
    vector<vector<VKeyTuple>> sliced_cmb = preproc(cmb, cores);

    /* Set openmp parameters and feed sliced key space in parallel to the improved filter */
    vector<vector<State>> r;
    omp_set_num_threads(cores);

#pragma omp parallel for ordered
    for(size_t i = 0x0; i < sliced_cmb.size(); ++i)
    {
        vector<State> v = improved_filter(c, d, sliced_cmb[i], l);
#pragma omp ordered
        r.push_back(v);
    }
    printf("Done.\n");

    /* Post-processing */
    vector<State> v = postproc(r);
    printf("Size of keyspace: %lu = 2^%f \n", v.size(), log2(v.size()));
    return v;
}

DiffStat differentials(State &c, State &d, const size_t l) 
{
    /* Choose inverse deltas depending on the fault location 'l' */
    const uint8_t** gm = ideltas1[map_fault[l]];

    /* Init differential matrix */
    DiffStat x;
    for(size_t i = 0x0; i < 0x10; ++i)
    {
        multimap<uint8_t, uint8_t> v;
        x[i] = v;
    }

    for(size_t i = 0x0; i < 0x100; ++i)
    {
        uint8_t k = (uint8_t) i;

        x[0x0].insert(pair<uint8_t, uint8_t>(EQ(c[0x0], d[0x0], k, gm[0x0]), k));
        x[0x1].insert(pair<uint8_t, uint8_t>(EQ(c[0x1], d[0x1], k, gm[0x1]), k));
        x[0x2].insert(pair<uint8_t, uint8_t>(EQ(c[0x2], d[0x2], k, gm[0x2]), k));
        x[0x3].insert(pair<uint8_t, uint8_t>(EQ(c[0x3], d[0x3], k, gm[0x3]), k));

        x[0x4].insert(pair<uint8_t, uint8_t>(EQ(c[0x4], d[0x4], k, gm[0x4]), k));
        x[0x5].insert(pair<uint8_t, uint8_t>(EQ(c[0x5], d[0x5], k, gm[0x5]), k));
        x[0x6].insert(pair<uint8_t, uint8_t>(EQ(c[0x6], d[0x6], k, gm[0x6]), k));
        x[0x7].insert(pair<uint8_t, uint8_t>(EQ(c[0x7], d[0x7], k, gm[0x7]), k));

        x[0x8].insert(pair<uint8_t, uint8_t>(EQ(c[0x8], d[0x8], k, gm[0x8]), k));
        x[0x9].insert(pair<uint8_t, uint8_t>(EQ(c[0x9], d[0x9], k, gm[0x9]), k));
        x[0xa].insert(pair<uint8_t, uint8_t>(EQ(c[0xa], d[0xa], k, gm[0xa]), k));
        x[0xb].insert(pair<uint8_t, uint8_t>(EQ(c[0xb], d[0xb], k, gm[0xb]), k));

        x[0xc].insert(pair<uint8_t, uint8_t>(EQ(c[0xc], d[0xc], k, gm[0xc]), k));
        x[0xd].insert(pair<uint8_t, uint8_t>(EQ(c[0xd], d[0xd], k, gm[0xd]), k));
        x[0xe].insert(pair<uint8_t, uint8_t>(EQ(c[0xe], d[0xe], k, gm[0xe]), k));
        x[0xf].insert(pair<uint8_t, uint8_t>(EQ(c[0xf], d[0xf], k, gm[0xf]), k));
    }
    return x;
}

DiffStat standard_filter(DiffStat x)
{
    /* Iterate over columns */
    for(size_t i = 0x0; i < 0x4; ++i)
    {
        /* Init set of valid indices */
        set<uint8_t> s;
        for(size_t j = 0x0; j < 0x100; ++j)
        {
            s.insert((uint8_t) j); 
        }

        /* Compute valid indices for a given column i */
        for(size_t j = 0x0; j < 0x4; ++j)
        {
            vector<uint8_t> keys = getKeys(x[rb[i][j]]);
            set<uint8_t> t(begin(keys), end(keys));
            set<uint8_t> r;
            set_intersection(s.begin(), s.end(), t.begin(), t.end(), inserter(r, r.begin()));
            s = r;
        }

        /* Erase invalid entries */
        for(size_t j = 0x0; j < 0x4; ++j)
        {
            vector<uint8_t> keys = getKeys(x[rb[i][j]]);
            set<uint8_t> t(begin(keys), end(keys));
            set<uint8_t> r;
            set_difference(t.begin(), t.end(), s.begin(), s.end(), inserter(r, r.begin()));
            vector<uint8_t> toDel(begin(r), end(r));
            for(size_t k = 0x0; k < toDel.size(); ++k) 
            {
                x[rb[i][j]].erase(toDel[k]);
            }
        }
    }
    return x;
}

/* Computes the cartesian product of all remaining key candidates of related elements */
vector<VKeyTuple> combine(DiffStat m)
{
    vector<VKeyTuple> result;

    /* Iterate over columns */
    for(size_t i = 0x0; i < 0x4; ++i)
    {
        VKeyTuple v;

        auto w = m[rb[i][0x0]];
        auto x = m[rb[i][0x1]];
        auto y = m[rb[i][0x2]];
        auto z = m[rb[i][0x3]];

        vector<uint8_t> wk = getKeys(w);

        for(size_t j = 0x0; j < wk.size(); ++j)
        {
            /* Extract elements with the same key and combine them */
            auto retw = w.equal_range(wk[j]);
            auto retx = x.equal_range(wk[j]);
            auto rety = y.equal_range(wk[j]);
            auto retz = z.equal_range(wk[j]);

            for(auto iw = retw.first; iw != retw.second; ++iw)
            {
                for(auto ix = retx.first; ix != retx.second; ++ix)
                {
                    for(auto iy = rety.first; iy != rety.second; ++iy)
                    {
                        for(auto iz = retz.first; iz != retz.second; ++iz)
                        {
                            KeyTuple t;
                            t[0x0] = iw->second;
                            t[0x1] = ix->second;
                            t[0x2] = iy->second;
                            t[0x3] = iz->second;
                            v.push_back(t);
                        }
                    }
                }
            }
        }
        result.push_back(v);
    }
    return result;
}

/* Prepare data for application of improved filter on multiple cores */
vector<vector<VKeyTuple>> preproc(vector<VKeyTuple> cmb, size_t cores)
{
    vector<VKeyTuple> slices;
    for(size_t i = 0x0; i < cores; ++i)
    {
        VKeyTuple v; 
        slices.push_back(v);
    }

    size_t n = cmb[0x0].size() / cores;
    size_t m = cmb[0x0].size() % cores;

    /* Distribute elements of the first column */
    for(size_t i = 0x0; i < n; ++i)
    {
        for(size_t j = 0x0; j < cores; ++j)
        {
            slices[j].push_back(cmb[0x0][cores * i + j]);
        }
    }
    
    /* Distribute the remaining elements, too */
    for(size_t j = 0x0; j < m; ++j)
    {
        slices[j].push_back(cmb[0x0][cores * n + j]);
    }

    vector<vector<VKeyTuple>> slices_cmb;
    for (size_t i = 0x0; i < slices.size(); ++i)
    {
        vector<VKeyTuple> v;
        v.push_back(slices[i]);
        v.push_back(cmb[0x1]);
        v.push_back(cmb[0x2]);
        v.push_back(cmb[0x3]);
        slices_cmb.push_back(v);
    }
    return slices_cmb;
}

vector<State> improved_filter(State &c, State &d, vector<VKeyTuple> &v, const size_t l)
{
    /* Configure fault equations depending on the fault location 'l' */
    const uint8_t** gm = ideltas2[l % 0x4];      // inverse deltas
    const size_t* x = indices_x[map_fault[l]];   // indices for c, d and k
    const size_t* y = indices_y[map_fault[l]];   // indices for h

    vector<State> candidates;

    for (size_t i0 = 0x0; i0 < v[0x0].size(); ++i0)
    {
        for (size_t i1 = 0x0; i1 < v[0x1].size(); ++i1)
        {
            for (size_t i2 = 0x0; i2 < v[0x2].size(); ++i2)
            {
                for (size_t i3 = 0x0; i3 < v[0x3].size(); ++i3)
                {
                    /* Index order of the tuples in 'key': (0x0, 0x7, Oxa, Oxd), (0x1, 0x4, Oxb, Oxe), (0x2, 0x5, 0x8, 0xf), (0x3, 0x6, 0x9, 0xc) */
                    
                    /* 10-th round key */
                    State k = 
                    {
                        v[0x0][i0][0x0], v[0x1][i1][0x0], v[0x2][i2][0x0], v[0x3][i3][0x0],   /* 0x0 0x1 0x2 0x3 */
                        v[0x1][i1][0x1], v[0x2][i2][0x1], v[0x3][i3][0x1], v[0x0][i0][0x1],   /* 0x4 0x5 0x6 0x7 */
                        v[0x2][i2][0x2], v[0x3][i3][0x2], v[0x0][i0][0x2], v[0x1][i1][0x2],   /* 0x8 0x9 0xa 0xb */
                        v[0x3][i3][0x3], v[0x0][i0][0x3], v[0x1][i1][0x3], v[0x2][i2][0x3]    /* 0xc 0xd 0xe 0xf */
                    };

                    /* 9-th round key */
                    State h;
                    h[0x0] = k[0x0] ^ sbox[k[0x9] ^ k[0xd]] ^ 0x36;
                    h[0x1] = k[0x1] ^ sbox[k[0xa] ^ k[0xe]];
                    h[0x2] = k[0x2] ^ sbox[k[0xb] ^ k[0xf]];
                    h[0x3] = k[0x3] ^ sbox[k[0x8] ^ k[0xc]];
                    h[0x4] = k[0x0] ^ k[0x4];
                    h[0x5] = k[0x1] ^ k[0x5];
                    h[0x6] = k[0x2] ^ k[0x6];
                    h[0x7] = k[0x3] ^ k[0x7];
                    h[0x8] = k[0x4] ^ k[0x8];
                    h[0x9] = k[0x5] ^ k[0x9];
                    h[0xa] = k[0x6] ^ k[0xa];
                    h[0xb] = k[0x7] ^ k[0xb];
                    h[0xc] = k[0x8] ^ k[0xc];
                    h[0xd] = k[0x9] ^ k[0xd];
                    h[0xe] = k[0xa] ^ k[0xe];
                    h[0xf] = k[0xb] ^ k[0xf];

                    uint8_t f0 = gm[0x0][
                        isbox[
                            gm_0e[isbox[c[x[0x0]] ^ k[x[0x0]]] ^ h[y[0x0]]] ^ gm_0b[isbox[c[x[0x1]] ^ k[x[0x1]]] ^ h[y[0x1]]] ^
                            gm_0d[isbox[c[x[0x2]] ^ k[x[0x2]]] ^ h[y[0x2]]] ^ gm_09[isbox[c[x[0x3]] ^ k[x[0x3]]] ^ h[y[0x3]]]
                        ] ^ isbox[
                            gm_0e[isbox[d[x[0x0]] ^ k[x[0x0]]] ^ h[y[0x0]]] ^ gm_0b[isbox[d[x[0x1]] ^ k[x[0x1]]] ^ h[y[0x1]]] ^
                            gm_0d[isbox[d[x[0x2]] ^ k[x[0x2]]] ^ h[y[0x2]]] ^ gm_09[isbox[d[x[0x3]] ^ k[x[0x3]]] ^ h[y[0x3]]]
                        ]
                    ];

                    uint8_t f1 = gm[0x1][
                        isbox[
                            gm_09[isbox[c[x[0x4]] ^ k[x[0x4]]] ^ h[y[0x4]]] ^ gm_0e[isbox[c[x[0x5]] ^ k[x[0x5]]] ^ h[y[0x5]]] ^
                            gm_0b[isbox[c[x[0x6]] ^ k[x[0x6]]] ^ h[y[0x6]]] ^ gm_0d[isbox[c[x[0x7]] ^ k[x[0x7]]] ^ h[y[0x7]]]
                        ] ^ isbox[
                            gm_09[isbox[d[x[0x4]] ^ k[x[0x4]]] ^ h[y[0x4]]] ^ gm_0e[isbox[d[x[0x5]] ^ k[x[0x5]]] ^ h[y[0x5]]] ^
                            gm_0b[isbox[d[x[0x6]] ^ k[x[0x6]]] ^ h[y[0x6]]] ^ gm_0d[isbox[d[x[0x7]] ^ k[x[0x7]]] ^ h[y[0x7]]]
                        ]
                    ];

                    uint8_t f2 = gm[0x2][
                        isbox[
                            gm_0d[isbox[c[x[0x8]] ^ k[x[0x8]]] ^ h[y[0x8]]] ^ gm_09[isbox[c[x[0x9]] ^ k[x[0x9]]] ^ h[y[0x9]]] ^
                            gm_0e[isbox[c[x[0xa]] ^ k[x[0xa]]] ^ h[y[0xa]]] ^ gm_0b[isbox[c[x[0xb]] ^ k[x[0xb]]] ^ h[y[0xb]]]
                        ] ^ isbox[
                            gm_0d[isbox[d[x[0x8]] ^ k[x[0x8]]] ^ h[y[0x8]]] ^ gm_09[isbox[d[x[0x9]] ^ k[x[0x9]]] ^ h[y[0x9]]] ^
                            gm_0e[isbox[d[x[0xa]] ^ k[x[0xa]]] ^ h[y[0xa]]] ^ gm_0b[isbox[d[x[0xb]] ^ k[x[0xb]]] ^ h[y[0xb]]]
                        ]
                    ];

                    uint8_t f3 = gm[0x3][
                        isbox[
                            gm_0b[isbox[c[x[0xc]] ^ k[x[0xc]]] ^ h[y[0xc]]] ^ gm_0d[isbox[c[x[0xd]] ^ k[x[0xd]]] ^ h[y[0xd]]] ^
                            gm_09[isbox[c[x[0xe]] ^ k[x[0xe]]] ^ h[y[0xe]]] ^ gm_0e[isbox[c[x[0xf]] ^ k[x[0xf]]] ^ h[y[0xf]]]
                        ] ^ isbox[
                            gm_0b[isbox[d[x[0xc]] ^ k[x[0xc]]] ^ h[y[0xc]]] ^ gm_0d[isbox[d[x[0xd]] ^ k[x[0xd]]] ^ h[y[0xd]]] ^
                            gm_09[isbox[d[x[0xe]] ^ k[x[0xe]]] ^ h[y[0xe]]] ^ gm_0e[isbox[d[x[0xf]] ^ k[x[0xf]]] ^ h[y[0xf]]]
                        ]
                    ];

                    if((f0 == f1) && (f1 == f2) && (f2 == f3))
                    {
                        candidates.push_back(k);
                    }
                }
            }
        }
    }
    return candidates;
}

/* Post-processing of subkey candidates */
vector<State> postproc(vector<vector<State>> &v)
{
    vector<State> r;
    for(size_t i = 0x0; i < v.size(); ++i)
    {
        for(size_t j = 0x0; j < v[i].size(); ++j)
        {
            r.push_back(reconstruct(v[i][j]));
        }
    }
    return r;
}

/* Reconstructs the master key from the 10-th round subkey */
State reconstruct(State &k)
{
    array<uint32_t, 0x2c> sk;
    sk.fill(0x0);

    for(size_t i = 0x0; i < 0x4; ++i)
    {
        for(size_t j = 0x0; j < 0x4; ++j)
        {
            sk[i] ^= k[0x4 * i + j] << (0x18 - 0x8 * j);
        }
    }

    size_t n = 0xa;
    for(size_t i = 0x0; i < n; ++i)
    {
        sk[0x4 * (i + 0x1)] = sk[0x4 * i] ^ ks_core(sk[0x4 * i + 0x2] ^ sk[0x4 * i + 0x3], n - i);
        sk[0x4 * (i + 0x1) + 0x1] = sk[0x4 * i] ^ sk[0x4 * i + 0x1];
        sk[0x4 * (i + 0x1) + 0x2] = sk[0x4 * i + 0x1] ^ sk[0x4 * i + 0x2];
        sk[0x4 * (i + 0x1) + 0x3] = sk[0x4 * i + 0x2] ^ sk[0x4 * i + 0x3];
    }

    n = sk.size();
    State mk;
    for(size_t i = 0x0; i < 0x4; ++i)
    {
        for(size_t j = 0x0; j < 0x4; ++j)
        {
            mk[0x4 * i + j] =  (sk[n - 0x4 + i] >> (24 - 0x8 * j)) & 0xff;
        }
    }

    return mk;
}

/* Core of AES key schedule */
uint32_t ks_core(uint32_t t, size_t r)
{
    uint32_t b = ROTL<uint32_t>(t, 0x8);
    uint32_t c = rcon[r] << 0x18;
    for(size_t i = 0x0; i < 0x4; ++i)
    {
        c ^= sbox[(b >> (0x18 - 0x8 * i)) & 0xff] << (0x18 - 0x8 * i);
    }
    return c;
}

vector<uint8_t> getKeys(multimap<uint8_t,uint8_t> m)
{
    vector<uint8_t> keys;
    for(auto it = m.begin(), end = m.end(); it != end; it = m.upper_bound(it->first))
    {
        keys.push_back(it->first);
    }
    return keys;
}

void printState(State x)
{
    for(size_t i = 0x0; i < x.size(); ++i)
    {
        printf("%02x", x[i]);
    }
    printf(" %lu", x.size());
}

vector<pair<pair<State, State>, State>> readfile(const string file, int bf)
{
    vector<pair<pair<State, State>, State>> pairs;

    ifstream infile(file);

    if(infile.is_open())
    {
        string in;
        while(getline(infile, in))
        {
            pair<pair<State, State>, State> cts;
            istringstream iin(in);
            string x, y, z;
            getline(iin, x, ' ');
            getline(iin, y, ' ');
            getline(iin, z, ' ');
            State c,d, e;
            c.fill(0x0);
            d.fill(0x0);
            e.fill(0x0);
            for(size_t i = 0x0; i < x.length(); i += 0x2)
            {
                c[i / 0x2] = strtol(x.substr(i , 0x2).c_str(), 0x0, 0x10);
                d[i / 0x2] = strtol(y.substr(i , 0x2).c_str(), 0x0, 0x10);
                if(bf)
                {
                    e[i / 0x2] = strtol(z.substr(i , 0x2).c_str(), 0x0, 0x10);
                }
            }
            cts.first.first = c;
            cts.first.second = d;
            cts.second = e;
            pairs.push_back(cts);
        }
        infile.close();
    }
    return pairs;
}

void writefile(State plaintext, State ciphertext, vector<State> keys, const string file)
{
    FILE * outfile;
    outfile = fopen(file.c_str(), "a");

    if((plaintext.empty() && !ciphertext.empty()) || (!plaintext.empty() && ciphertext.empty()))
    {
        printf("ERROR !!!\n");
        exit(0x1);
    }

    if(!plaintext.empty())
    {
        for(size_t i = 0x0; i < plaintext.size(); ++i)
        {
            fprintf(outfile, "%02x", plaintext[i]);
        }
        fprintf(outfile, "\n");
        for(size_t i = 0x0; i < ciphertext.size(); ++i)
        {
            fprintf(outfile, "%02x", ciphertext[i]);
        }
        fprintf(outfile, "\n");
    }
    
    for(size_t i = 0x0; i < keys.size(); ++i)
    {
        for(size_t j = 0x0; j < keys[i].size(); ++j)
        {
            fprintf(outfile, "%02x", keys[i][j]);
        }
        fprintf(outfile, "\n");
    }
    fclose(outfile);
}

void help()
{
    printf("Usage: ./dfa c l b f\n\n" );
    printf("Parameters\n");
    printf("%2sc: Number of cores >= 1.\n", "");
    printf("%2sl: Byte number of the AES state affected by the fault.\n%5sMust be in {-1, 0,..., 15}, where -1 means unknown.\n", "", "");
    printf("%2sb: Indicate if a brute-force search is needed over remainding master keys.\n%5sMust be 'bf' or 'nobf'.\n", "", "");
    printf("%2sf: Input file with one or more pairs of correct and faulty ciphertexts; and corresponding plaintext if 'bf'.\n\n","");
}

void printerror()
{
    printf("ERROR !!!\n");
    exit(0x1);
}

void convert(char* buff, uint8_t* data)
{
    for(size_t i = 0x0; i < 0x10; ++i)
    {
        if(buff[0x2 * i] >= 0x30 && buff[0x2 * i] <= 0x39)
        {
            data[i] = buff[0x2 * i] - 0x30;
        }
        if(buff[0x2 * i] >= 0x41 && buff[0x2 * i] <= 0x5a)
        {
            data[i] = buff[0x2 * i] - 0x41 + 0xa;
        }
        if(buff[0x2 * i] >= 0x61 && buff[0x2 * i] <= 0x7a)
        {
            data[i] = buff[0x2 * i] - 0x61 + 10;
        }
        
        if(buff[0x2 * i + 0x1] >= 0x30 && buff[0x2 * i + 0x1] <= 0x39)
        {
            data[i] |= (buff[0x2 * i + 0x1] - 0x30) << 0x4;
        }
        if(buff[0x2 * i + 0x1] >= 0x41 && buff[0x2 * i + 0x1] <= 0x5a)
        {
            data[i] |= (buff[0x2 * i + 0x1] - 0x41 + 0xa) << 0x4;
        }
        if(buff[0x2 * i + 0x1] >= 0x61 && buff[0x2 * i + 0x1] <= 0x7a)
        {
            data[i] |= (buff[0x2 * i + 0x1] - 0x61 + 10) << 0x4;
        }
    }
}

void bruteforce(const string name)
{
    char buff[0x21];
    uint8_t plaintext[0x10];
    uint8_t expected[0x10];
    uint8_t key[0x10];
    uint8_t* ciphertext;
    printf("PLOP\n");
    FILE* file = fopen(name.c_str(), "r");
    if(file == NULL)
    {
        printerror();
    }
    
    /* READ PLAINTEXT */
    if(fread(buff, 0x21, 0x1, file) != 0x1)
    {
        printerror();
    }
    convert(buff, plaintext);
    
    /* READ EXPECTED CIPHERTEXT */
    if(fread(buff, 0x21, 0x1, file) != 0x1)
    {
        printerror();
    }
    convert(buff, expected);
    
    /* TESTING EACH KEY */
    while(fread(buff, 0x21, 0x1, file) == 0x1)
    {
        convert(buff, key);
        ciphertext = encrypt(key, plaintext);
        if(!memcmp(ciphertext, expected, (size_t) sizeof(expected)))
        {
            printf("THE ONE KEY FOUND !!!\n");
            for(size_t i = 0x0; i < (size_t) sizeof(key); ++i)
            {
                printf("%02x", key[i]);
            }
            exit(0x0);
        }
    }
}

int main(int argc, char **argv)
{
    if(argc != 0x5)
    {
        help();
        return -0x1;
    }

    const size_t c = atoi(argv[0x1]);   // number of cores
    const int l = atoi(argv[0x2]);      // fault location
    const char* b = argv[0x3];          // brute-force
    const string f = argv[0x4];         // input file

    if(c < 0x0 || l < -0x1 || l > 0xf || (strcmp(b, "bf") && strcmp(b, "nobf")) || f.empty())
    {
        help();
        return -0x1;
    }

    vector<pair<pair<State, State>, State>> pairs = readfile(f, !strcmp(b, "bf"));

    /* Set fault location range */
    size_t j = 0x0;
    size_t n = 0x0;
    if(l == -0x1)
    {
        n = 16;
    } else {
        j = l;
        n = l + 0x1;
    }

    for(size_t i = 0x0; i < pairs.size(); ++i)
    {
        printf("(%lu) Analysing ciphertext pair:\n\n", i);
        printState(pairs[i].first.first);
        printf(" ");
        printState(pairs[i].first.second);
        printf("\n\nNumber of core(s): %lu \n", c);

        /* Create new output file */
        stringstream ss;
        ss << "res/" << i << ".csv";
        FILE * outfile;
        const string name = ss.str();
        outfile = fopen(name.c_str(), "w");
        fclose(outfile);

        size_t count = 0x0;
        while(j < n){
            printf("----------------------------------------------------\n");
            printf("Fault location: %lu\n", j);
            vector<State> keys = analyse(pairs[i].first.first, pairs[i].first.second, j, c);
            count += keys.size();
            State plaintext, expected;
            if(strcmp(b, "bf"))
            {
                plaintext = pairs[i].second;
                expected = pairs[i].first.first;
            }
            writefile(plaintext, expected, keys, name);
            j++;
        }
        if(l == -0x1)   // Reset j
        {
            j = 0x0;
        }
        else
        {
            j = l;
        }
        if(!strcmp(b, "bf"))
        {
            bruteforce(name);
        }
        printf("\n\n%lu masterkeys written to %s\n\n", count, name.c_str());
    }
    return 0x0;
}
