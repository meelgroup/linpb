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

#include "Solver.hpp"
#include "typedefs.hpp"

namespace rs {

namespace run {

extern Solver solver;

template <typename SMALL, typename LARGE>
class Optimization {
  const CePtr<ConstrExp<SMALL, LARGE>> origObj;

  LARGE upper_bound;
  ID lastUpperBound = ID_Undef;
  ID lastUpperBoundUnprocessed = ID_Undef;

 public:
  Optimization(CePtr<ConstrExp<SMALL, LARGE>> obj) : origObj(obj) {
    assert(origObj->vars.size() > 0);
    // NOTE: -origObj->getDegree() keeps track of the offset of the reformulated objective (or after removing unit lits)
    upper_bound = origObj->absCoeffSum() - origObj->getRhs() + 1;
  };

  LARGE normalizedUpperBound() { return upper_bound + origObj->getDegree(); }

  void printObjBounds() {
    if (options.verbosity.get() == 0) return;
    std::cout << "c bounds ";
    if (solver.foundSolution()) {
      std::cout << bigint(upper_bound);  // TODO: remove bigint(...) hack
    } else {
      std::cout << "-";
    }
    std::cout << " @ " << stats.getTime() << "\n";
  }

  void handleNewSolution(const std::vector<Lit>& sol) {
    [[maybe_unused]] LARGE prev_val = upper_bound;
    upper_bound = -origObj->getRhs();
    for (Var v : origObj->vars) upper_bound += origObj->coefs[v] * (int)(sol[v] > 0);
    assert(upper_bound < prev_val);

    CePtr<ConstrExp<SMALL, LARGE>> aux = solver.cePools.take<SMALL, LARGE>();
    origObj->copyTo(aux);
    aux->invert();
    aux->addRhs(-upper_bound + 1);
    solver.dropExternal(lastUpperBound, true, true);
    std::pair<ID, ID> res = solver.addConstraint(aux, Origin::UPPERBOUND);
    lastUpperBoundUnprocessed = res.first;
    lastUpperBound = res.second;
    if (lastUpperBound == ID_Unsat) quit::exit_UNSAT(solver, upper_bound);
  }

  void optimize() {
    SolveState reply = SolveState::SAT;

    while (true) {
      if (reply != SolveState::INPROCESSED && reply != SolveState::RESTARTED) printObjBounds();

      SolveAnswer out = aux::timeCall<SolveAnswer>([&] { return solver.solve(); }, stats.SOLVETIME);
      reply = out.state;
      if (reply == SolveState::RESTARTED) continue;
      if (reply == SolveState::UNSAT) {
        printObjBounds();
        quit::exit_UNSAT(solver, upper_bound);
      }
      assert(solver.decisionLevel() == 0);
      if (reply == SolveState::SAT) {
        assert(solver.foundSolution());
        ++stats.NSOLS;
        handleNewSolution(out.solution);
      } else {
        assert(reply == SolveState::INPROCESSED);  // keep looping
      }
    }
  }
};

void decide();
void run(CeArb objective);

}  // namespace run

}  // namespace rs
