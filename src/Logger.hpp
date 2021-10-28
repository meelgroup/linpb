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

#include <fstream>
#include "Stats.hpp"
#include "typedefs.hpp"

namespace rs {

struct Logger {
  std::ofstream formula_out;
  std::ofstream proof_out;
  ID last_formID = 0;
  ID last_proofID = 0;
  std::vector<ID> unitIDs;

  Logger(const std::string& proof_log_name) {
    formula_out = std::ofstream(proof_log_name + ".formula");
    formula_out << "* #variable= 0 #constraint= 0\n";
    formula_out << " >= 0 ;\n";
    ++last_formID;
    proof_out = std::ofstream(proof_log_name + ".proof");
    proof_out << "pseudo-Boolean proof version 1.0\n";
    proof_out << "l 1\n";
    ++last_proofID;
  }

  void flush() {
    formula_out.flush();
    proof_out.flush();
  }

  void logComment([[maybe_unused]] const std::string& comment, [[maybe_unused]] const Stats& sts) {
#if !NDEBUG
    proof_out << "* " << sts.getDetTime() << " " << comment << "\n";
#endif
  }

  void logComment([[maybe_unused]] const std::string& comment) {
#if !NDEBUG
    proof_out << "* "
              << " " << comment << "\n";
#endif
  }
};

}  // namespace rs
