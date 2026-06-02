# TTC: Toolbox for Theory Counting

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

TTC is a C++ tool for counting models of SMT-LIB formulas.  Given a formula
with projection variables prefixed by `proj`, TTC reports the number of models
projected onto those variables.  The tool can perform exact or approximate
counting and uses tree decompositions and component caching under the hood.

`ttc` also covers a set of experimental features like sampling and exact projected counting. See [extras.md]() for more details.

This complements our tool [csb](https://github.com/meelgroup/ttc), which is designed for counting over bitvectors.

## How to Build a Binary
To build on Linux, install the required dependencies:
```
sudo apt install build-essential cmake git libboost-program-options-dev liblpsolve55-dev libeigen3-dev libisl-dev libgmp-dev
```

Clone the repository and build the project:
```
git clone https://github.com/meelgroup/ttc
cd ttc
./configure.sh
cd build && make -j8
```

This compiles the tool against a modified version of `cvc5` from [this repository](https://github.com/meelgroup/ttc). If you want something else, please see `./configure.sh -h`, and [extras.md](). Detailed instructions on building on MacOS is also there.

## Usage: Volume Computation
Translate your problem into an SMT-LIB 2 file.  Mark the Boolean variables to be
projected with the prefix `proj` and pass the file to TTC:
```
$ ./ttc path/to/formula.smt2
...
s mc 42
...
```

## Usage: Projected Counting
By default if there are combination of theories, where a subset of theories is BV or Boolean, then `ttc` considers those BV/Boolean variables as projection variables and computes the projected count. If you want a subset of variables as projection variables then use `declare-projvar` keyword. If your formula is a pure LRA formula, and you want the counts to be projected on the Boolean variables, the use `-P` options.

### Example
```
$ ./ttc example/test.smt2
...
s mc 1
...
```



## Issues, questions, bugs, etc.
Please [open an issue](https://github.com/meelgroup/ttc/issues/new) for
support or bug reports. Please get in touch with Arijit Shaw or Kuldeep Meel in case you want a new feature / theory to be covered.

## How to Cite
If you use TTC in your research, please cite the following:


- **Efficient Volume Computation for SMT Formulas**

  Arijit Shaw, Uddalok Sarkar, and Kuldeep S. Meel

  *In Proceedings of Knowledge Representation and Reasoning (KR), Nov 2025*

- **Approximate SMT Counting Beyond Discrete Domains**

  Arijit Shaw, and Kuldeep S. Meel

  *In Proceedings of Design Automation Conference (DAC), Jun 2025*

