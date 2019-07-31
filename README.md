#### C++ implementation of DFA of AES-128 using a single fault injection.
---

#### LICENSE
This program is free software; see [LICENSE][3] for more details.

#### REQUIREMENTS
* gcc-4.8.1 or newer.
* Multi-core support via [OpenMP](http://openmp.org/).
* Workstation with at least 32 cores (recommended).
* Python 3.0 or newer.
* Numpy 1.16 or newer.

#### ATTACK
This repository is the implementation of the DFA attack from Tunstall, Mukhopadhyay and Ali using a single fault injection.

To execute the attack, some conditions are necessary. Since it's a DFA attack, a pair of ciphertext are needed as following:
* a correct ciphertext from the target AES device,
* a faulty ciphertext where one byte is disturbed between the 7th round *MixColumns* and the 8th round *MixColumns*.

The attack reduces to about 256 possible master keys.

For more informations the original can be found [here](https://eprint.iacr.org/2009/575.pdf).

#### USAGE
**Note:** Some Computations might take quite some time, especially in the case where the fault location is not known and if too few cores are available.

**Building**
```
make
```

**Cleaning**
```
make clean
```

**Single pair**

A fault was injected in byte 0 during the 8th round of the AES.
```
./tma_dfa 32 0 tests/single.csv
```

After the computation is finished all remaining master keys are written to the file `res/0.csv` including the correct one.

**Multiple pairs**

A fault was injected somewhere during the 8th round of the AES.
```
./tma_dfa 32 -1 tests/multiple.csv
```

After the computation is finished all remaining master keys are written to `res/{0, 1, 2}.csv`, *i.e.* one per pairs in the input file.

#### REFERENCES
[Original README](https://github.com/Daeinar/dfa-aes/blob/master/README.md)
