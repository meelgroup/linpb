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

#include "quit.hpp"
#include <iostream>
#include "Solver.hpp"
#include "globals.hpp"

namespace rs {

void quit::printSol(const std::vector<Lit>& sol) {
  printf("v");
  for (Var v = 1; v < (Var)sol.size(); ++v) printf(sol[v] > 0 ? " x%d" : " -x%d", v);
  printf("\n");
}

void quit::printSolAsOpb(const std::vector<Lit>& sol) {
  for (Var v = 1; v < (Var)sol.size(); ++v) {
    if (sol[v] > 0)
      std::cout << "+1 x" << v << " >= 1 ;\n";
    else
      std::cout << "-1 x" << v << " >= 0 ;\n";
  }
}

void quit::exit_SAT(const Solver& solver) {
  assert(solver.foundSolution());
  if (solver.logger) solver.logger->flush();
  if (options.verbosity.get() > 0) stats.print();
  puts("s SATISFIABLE");	
  if (options.printSol) printSol(solver.lastSol);
  exit(10);
}

template <typename LARGE>
void quit::exit_UNSAT(const Solver& solver, const LARGE& bestObjVal) {
  if (solver.logger) solver.logger->flush();
  if (options.verbosity.get() > 0) stats.print();
  if (solver.foundSolution()) {
    std::cout << "o " << bestObjVal << std::endl;
    std::cout << "s OPTIMUM FOUND" << std::endl;
    if (options.printSol) printSol(solver.lastSol);
    exit(30);
  } else {
    puts("s UNSATISFIABLE");
    exit(20);
  }
}
template void quit::exit_UNSAT<int>(const Solver& solver, const int& bestObjVal);
template void quit::exit_UNSAT<long long>(const Solver& solver, const long long& bestObjVal);
template void quit::exit_UNSAT<int128>(const Solver& solver, const int128& bestObjVal);
template void quit::exit_UNSAT<int256>(const Solver& solver, const int256& bestObjVal);
template void quit::exit_UNSAT<bigint>(const Solver& solver, const bigint& bestObjVal);

void quit::exit_UNSAT(const Solver& solver) { quit::exit_UNSAT<int>(solver, 0); }

void quit::exit_INDETERMINATE(const Solver& solver) {
  if (solver.foundSolution()) exit_SAT(solver);
  if (solver.logger) solver.logger->flush();
  if (options.verbosity.get() > 0) stats.print();
  puts("s UNKNOWN");
  exit(0);
}

DLL_PUBLIC void quit::exit_ERROR(const std::initializer_list<std::string>& messages) {
  std::cout << "Error: ";
  for (const std::string& m : messages) std::cout << m;
  std::cout << std::endl;
  exit(1);
}

}  // namespace rs
