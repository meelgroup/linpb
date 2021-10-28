/***********************************************************************
Copyright (c) 2014-2020, Jan Elffers
Copyright (c) 2019-2020, Jo Devriendt
Copyright (c) 2020, Stephan Gocht

Parts of the code were copied or adapted from MiniSat.

MiniSAT -- Copyright (c) 2003-2006, Niklas Een, Niklas Sorensson
           Copyright (c) 2007-2010  Niklas Sorensson

Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be included
in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
***********************************************************************/

#pragma once

#include <boost/multiprecision/cpp_int.hpp>
#if WITHGMP
#include <boost/multiprecision/gmp.hpp>
#endif  // WITHGMP
#include <cassert>
#include <exception>
#include <iostream>
#include <limits>
#include <unordered_map>
#include <vector>

namespace rs {

#if WITHGMP
using int128 = boost::multiprecision::int128_t;  // NOTE: a bit slower than __int128, but plays nice with mpz_int
using bigint = boost::multiprecision::mpz_int;   // NOTE: requires GMP
#else
using int128 = __int128;
using bigint = boost::multiprecision::cpp_int;
#endif  // WITHGMP
using int256 = boost::multiprecision::int256_t;
using BigCoef = bigint;
using BigVal = bigint;

using ID = uint64_t;
const ID ID_Undef = std::numeric_limits<ID>::max();
const ID ID_Unsat = ID_Undef - 1;
const ID ID_Trivial = 1;  // represents constraint 0 >= 0

using Var = int;
using Lit = int;
inline Var toVar(Lit l) { return std::abs(l); }

const int resize_factor = 2;

const int INF = 1e9 + 1;  // 1e9 is the maximum number of variables in the system, anything beyond is infinity
const long long INFLPINT = 1e15 + 1;  // based on max long range captured by double

const int limit32 = 1e9;         // 2^29-2^30
const long long limit64 = 2e18;  // 2^60-2^61
const double limit96 = 8e27;     // 2^92-2^93, so 46 bits is less than half
const double limit128 = 32e36;   // 2^124-2^125, so 62 bits is less than half
const double limit256 = 1e76;    // 2^252-2^253, so 126 bits is less than half
const int conflLimit32 = 14;
const int conflLimit64 = 30;
const int conflLimit96 = 46;
const int conflLimit128 = 62;

using IntVecIt = std::vector<int>::iterator;

using ActValV = long double;
const ActValV actLimitV = (ActValV)1e300 * (ActValV)1e300 * (ActValV)1e300 * (ActValV)1e300 * (ActValV)1e300 *
                          (ActValV)1e300 * (ActValV)1e300 * (ActValV)1e300;  // ~1e2400 << 2^(2^13)
using ActValC = float;
const ActValC actLimitC = 1e30;  // ~1e30 << 2^(2^7)

/*
 * UNKNOWN: uninitialized value
 * FORMULA: original input formula constraints
 * LEARNED: learned from regular conflict analysis
 * FARKAS: LP solver infeasibility witness
 * LEARNEDFARKAS: constraint learned from conflict analysis on FARKAS
 * GOMORY: Gomory cut
 * UPPERBOUND: upper and lower bounds on the objective function
 *
 * max number of types is 16, as the type is stored with 4 bits in Constr
 */
enum class Origin {
  UNKNOWN,
  FORMULA,
  LEARNED,
  FARKAS,
  LEARNEDFARKAS,
  GOMORY,
  UPPERBOUND,
  GAUSS,
};

template <typename SMALL, typename LARGE>
struct ConstrExp;
using ConstrExp32 = ConstrExp<int, long long>;
using ConstrExp64 = ConstrExp<long long, int128>;
using ConstrExp96 = ConstrExp<int128, int128>;
using ConstrExp128 = ConstrExp<int128, int256>;
using ConstrExpArb = ConstrExp<bigint, bigint>;
struct ConstrExpSuper;

template <typename CE>
struct CePtr;
using Ce32 = CePtr<ConstrExp32>;
using Ce64 = CePtr<ConstrExp64>;
using Ce96 = CePtr<ConstrExp96>;
using Ce128 = CePtr<ConstrExp128>;
using CeArb = CePtr<ConstrExpArb>;
using CeSuper = CePtr<ConstrExpSuper>;
using CeNull = CePtr<ConstrExp32>;

template <typename CF, typename DG>
struct ConstrSimple;
using ConstrSimple32 = ConstrSimple<int, long long>;
using ConstrSimple64 = ConstrSimple<long long, int128>;
using ConstrSimple96 = ConstrSimple<int128, int128>;
using ConstrSimple128 = ConstrSimple<int128, int256>;
using ConstrSimpleArb = ConstrSimple<bigint, bigint>;
struct ConstrSimpleSuper;

struct Constr;
struct Clause;
struct Cardinality;

template <typename CF, typename DG>
struct Counting;
using Counting32 = Counting<int, long long>;
using Counting64 = Counting<long long, int128>;
using Counting96 = Counting<int128, int128>;

template <typename CF, typename DG>
struct Watched;
using Watched32 = Watched<int, long long>;
using Watched64 = Watched<long long, int128>;
using Watched96 = Watched<int128, int128>;

template <typename CF, typename DG>
struct CountingSafe;
using CountingSafe32 = CountingSafe<int, long long>;
using CountingSafe64 = CountingSafe<long long, int128>;
using CountingSafe96 = CountingSafe<int128, int128>;
using CountingSafeArb = CountingSafe<bigint, bigint>;

template <typename CF, typename DG>
struct WatchedSafe;
using WatchedSafe32 = WatchedSafe<int, long long>;
using WatchedSafe64 = WatchedSafe<long long, int128>;
using WatchedSafe96 = WatchedSafe<int128, int128>;
using WatchedSafeArb = WatchedSafe<bigint, bigint>;

template <typename CF>
struct Term {
  Term() : c(0), l(0) {}
  Term(const CF& x, Lit y) : c(x), l(y) {}
  CF c;
  Lit l;
};

template <typename CF>
std::ostream& operator<<(std::ostream& o, const Term<CF>& t) {
  return o << t.c << "x" << t.l;
}

inline class AsynchronousInterrupt : public std::exception {
 public:
  virtual const char* what() const throw() { return "Program interrupted by user."; }
} asynchInterrupt;

}  // namespace rs
