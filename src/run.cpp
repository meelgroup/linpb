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

#include "run.hpp"
#include "ConstrExp.hpp"
#include "Solver.hpp"

namespace rs {

DLL_PUBLIC Solver run::solver;

void run::decide() {
  while (true) {
    SolveState reply = aux::timeCall<SolveState>([&] { return solver.solve().state; }, stats.SOLVETIME);
    if (reply == SolveState::SAT) {
      quit::exit_SAT(solver);
      // return when used as a library
      return;
    }
    else if (reply == SolveState::UNSAT) {
      quit::exit_UNSAT(solver);
      // return when used as a library
      return;
    }
  }
}

DLL_PUBLIC void run::run(CeArb objective) {
  stats.RUNSTARTTIME = aux::cpuTime();
  if (options.verbosity.get() > 0)
    std::cout << "c #variables " << solver.getNbOrigVars() << " #constraints " << solver.getNbConstraints()
              << std::endl;
  try {
    if (objective->vars.size() > 0) {
      objective->stopLogging();
      objective->removeUnitsAndZeroes(solver.getLevel(), solver.getPos(), false);

      BigVal maxVal = objective->getCutoffVal();
      if (maxVal <= limit32) {  // TODO: try to internalize this check in ConstrExp
        Ce32 result = solver.cePools.take32();
        objective->copyTo(result);
        Optimization optim(result);
        optim.optimize();
      } else if (maxVal <= limit64) {
        Ce64 result = solver.cePools.take64();
        objective->copyTo(result);
        Optimization optim(result);
        optim.optimize();
      } else if (maxVal <= BigVal(limit96)) {
        Ce96 result = solver.cePools.take96();
        objective->copyTo(result);
        Optimization optim(result);
        optim.optimize();
      } else if (maxVal <= BigVal(limit128)) {
        Ce128 result = solver.cePools.take128();
        objective->copyTo(result);
        Optimization optim(result);
        optim.optimize();
      } else {
        CeArb result = solver.cePools.takeArb();
        objective->copyTo(result);
        Optimization<bigint, bigint> optim(result);
        optim.optimize();
      }
    } else {
      decide();
    }
  } catch (const AsynchronousInterrupt& ai) {
    std::cout << "c " << ai.what() << std::endl;
    quit::exit_INDETERMINATE(solver);
  }
}

}  // namespace rs
