# TTC: Toolbox for Theory Counting

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

TTC is a C++ tool for counting models of SMT-LIB formulas.  Given a formula
with projection variables prefixed by `proj`, TTC reports the number of models
projected onto those variables.  The tool can perform exact or approximate
counting and uses tree decompositions and component caching under the hood.

## Documentation

A documentation on usage and building is available at https://arijitsh.github.io/docs-ttc/

## Different Usage

### Usage: Bitvector Counting
To counts the models of formulas over the theory of fixed-size bit-vectors (logics BV and UFBV), use this mode. The formula is eagerly bit-blasted to a model-preserving CNF in memory and counted with Arjun + ApproxMC in-process, so no intermediate file is written. This engine is auto-selected for bit-vector logics.


### Usage: Volume Computation
Computes the volume of the solution space of linear real arithmetic (QF_LRA) formulas — that is, the measure of the region of real assignments that satisfy the formula. The polytopes are estimated with the volesti random-walk sampler. This engine is auto-selected for QF_LRA inputs with no projection variables.

### Usage: Projected Counting
Hybrid formulas mix discrete and continuous variables. Engine 2 handles them with projection-based approximate counting (PACT): it counts the satisfying assignments to the discrete (bit-vector / Boolean) projection variables, while the continuous part is discharged by the SMT solver.

This engine is auto-selected for non-BV, non-LRA logics that have bit-vector or Boolean projection variables. 

### Usage: Skolem Function Counting

To count the interpretations (Skolem functions) of an uninterpreted function that satisfy a formula, use the SkolemFC mode. It is intended for QF_UFBV formulas that constrain the behaviour of an uninterpreted function.

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

