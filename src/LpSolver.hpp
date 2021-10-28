/***********************************************************************
Copyright (c) 2014-2020, Jan Elffers
Copyright (c) 2019-2020, Jo Devriendt
Copyright (c) 2020, Stephan Gocht

Parts of the code were copied or adapted from MiniSat.

MiniSat -- Copyright (c) 2003-2006, Niklas Een, Niklas Sorensson
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

#include "ConstrExp.hpp"
#include "ConstrSimple.hpp"
#include "aux.hpp"
#include "globals.hpp"
#include "typedefs.hpp"

#if WITHSOPLEX

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpragmas"
#pragma GCC diagnostic ignored "-Wunknown-warning-option"
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wstrict-overflow"
#pragma GCC diagnostic ignored "-Wdeprecated-copy"
#pragma GCC diagnostic ignored "-Wclass-memaccess"
#include "soplex.h"
#pragma GCC diagnostic pop

#endif  // WITHSOPLEX

namespace rs {

enum class LpStatus { INFEASIBLE, OPTIMAL, PIVOTLIMIT, UNDETERMINED };

struct RowData {
  ID id;
  bool removable;
  RowData(){};
  RowData(ID i, bool r) : id(i), removable(r){};
};

struct AdditionData {
  ConstrSimple64 cs;
  bool removable;
};

struct BoundData {
  ID id = ID_Trivial;
  ConstrSimple64 cs;
};

#if WITHSOPLEX

class LpSolver;
struct CandidateCut {
  ConstrSimple64 simpcons;
  CRef cr = CRef_Undef;
  double norm = 1;
  double ratSlack = 0;

  CandidateCut(){};
  CandidateCut(CeSuper in, const std::vector<double>& sol);
  CandidateCut(const Constr& in, CRef cr, const std::vector<double>& sol, ConstrExpPools& pools);
  double cosOfAngleTo(const CandidateCut& other) const;

 private:
  void initialize(const std::vector<double>& sol);
};
std::ostream& operator<<(std::ostream& o, const CandidateCut& cc);

class Solver;
class LpSolver {
  friend class Solver;
  friend struct CandidateCut;

  soplex::SoPlex lp;
  Solver& solver;

  double lpPivotMult = 1;
  constexpr static double INFTY = 1e100;
  constexpr static double maxMult = 1e15;  // sufficiently large to reduce rounding errors

  soplex::DVectorReal lpSol;
  std::vector<double> lpSolution;
  soplex::DVectorReal lpSlackSolution;
  soplex::DVectorReal lpMultipliers;
  soplex::DVectorReal upperBounds;
  soplex::DVectorReal lowerBounds;
  soplex::DSVectorReal lpRow;

  std::vector<RowData> row2data;
  std::vector<int> toRemove;  // rows
  std::unordered_map<ID, AdditionData> toAdd;
  BoundData boundsToAdd[2];  // [0] is upper bound, [1] lower bound

  std::vector<CandidateCut> candidateCuts;

 public:
  LpSolver(Solver& solver, const CeArb objective);
  void setNbVariables(int n);

  std::pair<LpStatus, CeSuper> checkFeasibility(bool inProcessing = false);  // TODO: don't use objective function here?
  void inProcess();

  void addConstraint(CeSuper c, bool removable, bool upperbound = false, bool lowerbound = false);
  void addConstraint(CRef cr, bool removable, bool upperbound = false, bool lowerbound = false);

 private:
  int getNbVariables() const;
  int getNbRows() const;

  void flushConstraints();

  void convertConstraint(const ConstrSimple64& c, soplex::DSVectorReal& row, double& rhs);
  void resetBasis();
  CeSuper createLinearCombinationFarkas(soplex::DVectorReal& mults);
  CandidateCut createLinearCombinationGomory(soplex::DVectorReal& mults);
  double getScaleFactor(soplex::DVectorReal& mults, bool removeNegatives);
  Ce64 rowToConstraint(int row);
  void constructGomoryCandidates();
  void constructLearnedCandidates();
  void addFilteredCuts();
  void pruneCuts();

  inline static double nonIntegrality(double a) { return aux::abs(std::round(a) - a); }
  inline static bool validVal(double a) { return std::round(a) == a && std::abs(a) < INFLPINT; }
  // NOTE: double type can only store ranges of integers up to ~9e15
};

#else

// TODO: check correspondence to above
class LpSolver {
 public:
  // LpSolver([[maybe_unused]] Solver& slvr, [[maybe_unused]] const CeArb objective){
  // See https://stackoverflow.com/questions/52263141/maybe-unused-and-constructors
  LpSolver(Solver& slvr, const CeArb objective) {
    _unused(slvr);
    _unused(objective);
  };
  void setNbVariables([[maybe_unused]] int n){};

  std::pair<LpStatus, CeSuper> checkFeasibility([[maybe_unused]] bool inProcessing = false) {
    assert(false);
    return {LpStatus::UNDETERMINED, CeNull()};
  }
  void inProcess() {}

  void addConstraint([[maybe_unused]] CeSuper c, [[maybe_unused]] bool removable, [[maybe_unused]] bool upperbound,
                     [[maybe_unused]] bool lowerbound) {}
  void addConstraint([[maybe_unused]] CRef cr, [[maybe_unused]] bool removable, [[maybe_unused]] bool upperbound,
                     [[maybe_unused]] bool lowerbound) {}
};

#endif  // WITHSOPLEX

}  // namespace rs
