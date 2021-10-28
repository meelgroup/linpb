# LinPB: PB-XOR Solver

[![CP](https://img.shields.io/badge/CP-2021-blue.svg)](https://drops.dagstuhl.de/opus/volltexte/2021/15349/)
[![Dataset](https://img.shields.io/badge/paper-Dataset-yellow.svg)](https://doi.org/10.5281/zenodo.5526835)
[![Application](https://img.shields.io/badge/application-ApproxMCPB-orange.svg)](https://github.com/meelgroup/approxmcpb)

SAT Solver for linear pseudo boolean and xor constraints.

The branch is for library usage. Underlying PB Solver: [RoundingSat](https://gitlab.com/miao_research/roundingsat). We support XOR constraints based on Gauss-Jordan Elimination which is adapted from [CryptoMiniSat](https://github.com/msoos/cryptominisat#gauss-jordan-elimination).


## Input format:
   - Pseudo Boolean constraints: [OPB format](InputFormats.md)
   - XOR constraints: `* xor <variable>+ 0|1`
   
### Example input file

```
* #variable= 4 #constraint= 2
*
* xor x1 x2 1
* xor x3 x4 0
* 
+1 x1 +2 x2 >= 1;
+1 x1 +2 x3 -3 x4 = -1;
```

Note: We recommend to encode short xors into pb constraints.

## How to Build

In the root directory of LinPB:

    cd build
    cmake -DCMAKE_BUILD_TYPE=Release ..
    make

For a debug build:

    cd build_debug
    cmake -DCMAKE_BUILD_TYPE=Debug ..
    make

For more builds, similar build directories can be created.

### Dependencies

- C++17 (i.e., a reasonably recent compiler)
- Boost library: https://www.boost.org
- Optionally: SoPlex LP solver (see below)

### SoPlex

RoundingSat underlying LinPB supports an integration with the LP solver SoPlex to improve its search routine.
For this, first download SoPlex at [here](https://soplex.zib.de/download.php?fname=soplex-5.0.1.tgz) and place the downloaded file in the root directory of LinPB.
Next, follow the above build process, but configure with the cmake option `-Dsoplex=ON`:

    cd build
    cmake -DCMAKE_BUILD_TYPE=Release -Dsoplex=ON ..
    make

The location of the SoPlex package can be configured with the cmake option `-Dsoplex_pkg=<location>`.

## How to Use

Run the following command:
```
$ ./build/linpb [OPTIONS] /path/to/file
```
Find more options: 
```
$ ./build/linpb --help
```

## References

### Our work: LinPB

**[YM21]** Jiong Yang, Kuldeep S. Meel. Engineering an Efficient PB-XOR Solver. *CP 2021*

### RoundingSat

Original paper with a focus on cutting planes conflict analysis:  
**[EN18]** J. Elffers, J. Nordström. Divide and Conquer: Towards Faster Pseudo-Boolean Solving. *IJCAI 2018*

Integration with SoPlex:  
**[DGN20]** J. Devriendt, A. Gleixner, J. Nordström. Learn to Relax: Integrating 0-1 Integer Linear Programming with Pseudo-Boolean Conflict-Driven Search. *CPAIOR 2020 / Constraints journal*

Watched propagation:  
**[D20]** J. Devriendt. Watched Propagation for 0-1 Integer Linear Constraints. *CP 2020*

### Gauss-Jordan Elimination:

**[SN09]** Soos M., Nohl K., Castelluccia C. Extending SAT Solvers to Cryptographic Problems. *SAT 2009*

**[HJ12]** Han CS., Jiang JH.R. When Boolean Satisfiability Meets Gaussian Elimination in a Simplex Way. *CAV 2012*

**[SGM20]** Soos M., Gocht S., Meel K.S. Tinted, Detached, and Lazy CNF-XOR Solving and Its Applications to Counting and Sampling. *CAV 2020*


