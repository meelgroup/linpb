# Input formats

## OPB format

All variables are boolean (have values 0 or 1).

A *literal* is written as `x5` or `~x10`. We have `~x10 = 1 - x10`.

A *linear inequality* is written as `{<coef> <literal> }* <cmp> <degree>;`.
- `<coef>` is a coefficient (a positive or negative integer).
- `<literal>` is a literal.
- `<cmp>` is either `>=` or `=`.
- `<degree>` is a positive or negative integer.

Non-linear terms are not supported.

An OPB file consists of a header line `* #variable= <N> #constraint= <M>`, comment lines (lines starting with `*`)
and linear inequalities as defined above (using the variables `x1` to `xN`). The number of constraints is ignored by the solver; it reads until EOF.

For more information, see the [OPB description](http://www.cril.univ-artois.fr/PB16/format.pdf) for the SAT competition.
Note that RoundingSat only supports linear OPB constraints.

### Example OPB file

```
* #variable= 12 #constraint= 7
*
* Pigeonhole principle formula for 4 pigeons and 3 holes
* Generated with `cnfgen`
* (C) 2012-2016 Massimo Lauria <lauria.massimo@gmail.com>
* https://massimolauria.github.io/cnfgen
*
* COMMAND LINE: cnfgen -of opb php 4 3
*
+1 x1 +1 x2 +1 x3 >= 1;
+1 x4 +1 x5 +1 x6 >= 1;
+1 x7 +1 x8 +1 x9 >= 1;
+1 x10 +1 x11 +1 x12 >= 1;
-1 x1 -1 x4 -1 x7 -1 x10 >= -1;
-1 x2 -1 x5 -1 x8 -1 x11 >= -1;
-1 x3 -1 x6 -1 x9 -1 x12 >= -1;
```

### Encoding a clause

A clause `x1 v ~x2` can be encoded as `1 x1 1 ~x2 >= 1;`.

### Optimization problems

It is possible to add an objective function to the OPB file. Roundingsat will then optimization this function.
The objective function has to be specified before all constraints. The format is

```
min: 1 x1 1 -x2 2 ~x3 ;
```

## CNF format

All variables are boolean, and all constraints are clauses.

A *literal* is written as `5` for `x5` or `-10` for `~x10`.

A *clause* is written as `{<literal> }* 0`.

A DIMACS file consists of a header line `p cnf <N> <M>` where `<N>` is the number of variables and `<M>` is the number of clauses,
comment lines starting with `c`, and clauses as defined above. The variables are `x1` to `xN`, and the number of clauses is ignored and the solver reads until EOF.

For more information, see the [CNF description](https://www.satcompetition.org/2011/format-benchmarks2011.html) for the SAT competition.

### Example CNF file

```
c Pigeonhole principle formula for 4 pigeons and 3 holes
c Generated with `cnfgen`
c (C) 2012-2016 Massimo Lauria <lauria.massimo@gmail.com>
c https://massimolauria.github.io/cnfgen
c
c COMMAND LINE: cnfgen php 4 3
c
p cnf 12 22
1 2 3 0
4 5 6 0
7 8 9 0
10 11 12 0
-1 -4 0
-1 -7 0
-1 -10 0
-4 -7 0
-4 -10 0
-7 -10 0
-2 -5 0
-2 -8 0
-2 -11 0
-5 -8 0
-5 -11 0
-8 -11 0
-3 -6 0
-3 -9 0
-3 -12 0
-6 -9 0
-6 -12 0
-9 -12 0
```

### WCNF format

An extension of the CNF format to allow for soft (weighted) clauses. See the [WCNF description](https://maxsat-evaluations.github.io/2018/rules.html#input) for the MaxSAT competition.
