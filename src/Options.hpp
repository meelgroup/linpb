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

#include "aux.hpp"
#include "quit.hpp"

namespace rs {

class Option {
 public:
  const std::string name = "";
  const std::string description = "";

  Option(const std::string& n, const std::string& d) : name(n), description(d) {}

  virtual void printUsage(int colwidth) const = 0;
  virtual void parse(const std::string& v) = 0;
};

class VoidOption : public Option {
  bool val = false;

 public:
  VoidOption(const std::string& n, const std::string& d) : Option{n, d} {}

  explicit operator bool() const { return val; }

  void printUsage(int colwidth) const override {
    std::cout << " --" << name;
    for (int i = name.size(); i < colwidth + 3; ++i) std::cout << " ";
    std::cout << description << "\n";
  }

  void parse([[maybe_unused]] const std::string& v) override {
    assert(v == "");
    val = true;
  }
};

class BoolOption : public Option {
  bool val = false;

 public:
  BoolOption(const std::string& n, const std::string& d, bool v) : Option{n, d}, val(v) {}

  explicit operator bool() const { return val; }

  void printUsage(int colwidth) const override {
    std::cout << " --" << name << "=? ";
    for (int i = name.size(); i < colwidth; ++i) std::cout << " ";
    std::cout << description << " (0 or 1; default " << val << ")\n";
  }

  void parse(const std::string& v) override {
    try {
      val = std::stod(v);
    } catch (const std::invalid_argument& ia) {
      quit::exit_ERROR({"Invalid value for ", name, ": ", v, ".\nCheck usage with --help option."});
    }
    if (val != 0 && val != 1)
      quit::exit_ERROR({"Invalid value for ", name, ": ", v, ".\nCheck usage with --help option."});
  }
};

template <typename T>
class ValOption : public Option {
  T val;
  std::string checkDescription;
  std::function<bool(const T&)> check;

 public:
  ValOption(const std::string& n, const std::string& d, const T& v, const std::string& cd,
            const std::function<bool(const T&)>& c)
      : Option{n, d}, val(v), checkDescription(cd), check(c) {}

  T& get() { return val; }

  void printUsage(int colwidth) const override {
    std::cout << " --" << name << "=? ";
    for (int i = name.size(); i < colwidth; ++i) std::cout << " ";
    std::cout << description << " (" << checkDescription << "; default " << val << ")\n";
  }

  void parse(const std::string& v) override {
    try {
      val = aux::sto<T>(v);
    } catch (const std::invalid_argument& ia) {
      quit::exit_ERROR({"Invalid value for ", name, ": ", v, ".\nCheck usage with --help option."});
    }
    if (!check(val)) quit::exit_ERROR({"Invalid value for ", name, ": ", v, ".\nCheck usage with --help option."});
  }
};

class EnumOption : public Option {
  std::string val;
  std::vector<std::string> values;

 public:
  EnumOption(const std::string& n, const std::string& d, const std::string& v, const std::vector<std::string>& vs)
      : Option{n, d}, val(v), values(vs) {}

  void printUsage(int colwidth) const override {
    std::cout << " --" << name << "=? ";
    for (int i = name.size(); i < colwidth; ++i) std::cout << " ";
    std::cout << description << " (";
    for (int i = 0; i < (int)values.size(); ++i) {
      if (i > 0) std::cout << ", ";
      if (values[i] == val) std::cout << "default=";
      std::cout << values[i];
    }
    std::cout << ")\n";
  }

  void parse(const std::string& v) override {
    if (std::find(std::begin(values), std::end(values), v) == std::end(values))
      quit::exit_ERROR({"Invalid value for ", name, ": ", v, ".\nCheck usage with --help option."});
    val = v;
  }

  bool is(const std::string& v) const {
    assert(std::find(std::begin(values), std::end(values), v) != std::end(values));
    return val == v;
  }
};

struct Options {
  VoidOption help{"help", "Print this help message"};
  BoolOption printSol{"print-sol", "Print the solution if found", 0};
  ValOption<int> verbosity{"verbosity", "Verbosity of the output", 1, "0 =< int",
                           [](const int& x) -> bool { return x >= 0; }};
  ValOption<std::string> proofLog{"proof-log", "Filename for the proof logs, left unspecified disables proof logging",
                                  "", "/path/to/file", [](const std::string&) -> bool { return true; }};
  ValOption<double> lubyBase{"luby-base", "Base of the Luby restart sequence", 2, "1 =< float",
                             [](double x) -> bool { return 1 <= x; }};
  ValOption<int> lubyMult{"luby-mult", "Multiplier of the Luby restart sequence", 100, "1 =< int",
                          [](const int& x) -> bool { return x >= 1; }};
  ValOption<double> varDecay{"vsids-var", "VSIDS variable decay factor", 0.95, "0.5 =< float < 1",
                             [](const double& x) -> bool { return 0.5 <= x && x < 1; }};
  ValOption<double> clauseDecay{"vsids-clause", "VSIDS clause decay factor", 0.999, "0.5 =< float < 1",
                                [](const double& x) -> bool { return 0.5 <= x && x < 1; }};
  ValOption<int> dbCleanInc{"db-inc", "Database cleanup interval increment", 100, "1 =< int",
                            [](const int& x) -> bool { return 1 <= x; }};
  ValOption<double> propCounting{"prop-counting", "Counting propagation instead of watched propagation", 0.7,
                                 "0 (no counting) =< float =< 1 (always counting)",
                                 [](const double& x) -> bool { return 0 <= x && x <= 1; }};
  BoolOption propClause{"prop-clause", "Optimized two-watched propagation for clauses", 1};
  BoolOption propCard{"prop-card", "Optimized two-watched propagation for clauses", 1};
  BoolOption propIdx{"prop-idx", "Optimize index of watches during propagation", 1};
  BoolOption propSup{"prop-sup", "Avoid superfluous watch checks", 1};
  ValOption<double> lpPivotRatio{
      "lp", "Ratio of #pivots/#conflicts limiting LP calls (negative means infinite, 0 means no LP solving)", 1,
      "-1 =< float", [](const double& x) -> bool { return x >= -1; }};
  ValOption<int> lpPivotBudget{"lp-budget", "Base LP call pivot budget", 1000, "1 =< int",
                               [](const int& x) -> bool { return x >= 1; }};
  ValOption<double> lpIntolerance{"lp-intolerance", "Intolerance for floating point artifacts", 1e-6, "0 < float",
                                  [](const double& x) -> bool { return x > 1; }};
  BoolOption addGomoryCuts{"lp-cut-gomory", "Generate Gomory cuts", 1};
  BoolOption addLearnedCuts{"lp-cut-learned", "Use learned constraints as cuts", 1};
  ValOption<int> gomoryCutLimit{"lp-cut-gomlim",
                                "Max number of tableau rows considered for Gomory cuts in a single round", 100,
                                "1 =< int", [](const int& x) -> bool { return 1 <= x; }};
  ValOption<double> maxCutCos{
      "lp-cut-maxcos",
      "Upper bound on cosine of angle between cuts added in one round, higher means cuts can be more parallel", 0.1,
      "0 =< float =< 1", [](const double& x) -> bool { return 0 <= x && x <= 1; }};
  BoolOption slackdiv{"ca-slackdiv",
                      "Use slack+1 as divisor for reason constraints (instead of the asserting coefficient)", 0};
  BoolOption weakenFull{"ca-weaken-full",
                        "Weaken non-divisible non-falsified literals in reason constraints completely", 0};
  BoolOption weakenNonImplying{"ca-weaken-nonimplying",
                               "Weaken non-implying falsified literals from reason constraints", 0};
  BoolOption bumpOnlyFalse{"bump-onlyfalse",
                           "Bump activity of literals encountered during conflict analysis only when falsified", 0};
  BoolOption bumpCanceling{"bump-canceling",
                           "Bump activity of literals encountered during conflict analysis when canceling", 1};
  BoolOption bumpLits{
      "bump-lits",
      "Bump activity of literals encountered during conflict analysis, twice when occurring with opposing sign", 1};
  ValOption<int> bitsOverflow{"bits-overflow",
                              "Bit width of maximum coefficient during conflict analysis calculations (0 is unlimited, "
                              "unlimited or greater than 62 may use slower arbitrary precision implementations)",
                              62, "0 =< int", [](const int& x) -> bool { return x >= 0; }};
  ValOption<int> bitsReduced{"bits-reduced",
                             "Bit width of maximum coefficient after reduction when exceeding bits-overflow (0 is "
                             "unlimited, 1 reduces to cardinalities)",
                             29, "0 =< int", [](const int& x) -> bool { return x >= 0; }};
  ValOption<int> bitsLearned{
      "bits-learned",
      "Bit width of maximum coefficient for learned constraints (0 is unlimited, 1 reduces to cardinalities)", 29,
      "0 =< int", [](const int& x) -> bool { return x >= 0; }};
  ValOption<int> bitsInput{
      "bits-input",
      "Bit width of maximum coefficient for input constraints (0 is unlimited, 1 allows only cardinalities)", 0,
      "0 =< int", [](const int& x) -> bool { return x >= 0; }};

  ValOption<int> configGJE{"config-GJE",
        	      "Configure the tactic of Gaussian Jordan Elimination. "
		      "(0 disables GJE and uses PB encoding of XOR, 1 is lazy GJE, 2 is eager GJE, 3 is shared watches, 4 is mixed watches)",
		      4,
		      "0 =< int =< 4", [](const int& x) -> bool { return 0 <= x && x <= 4; }};

  const std::vector<Option*> options = {
      &help,           &printSol,       &verbosity,    &proofLog,      &lubyBase,      &lubyMult,
      &varDecay,       &clauseDecay,    &dbCleanInc,   &propCounting,  &propClause,    &propCard,
      &propIdx,        &propSup,        &lpPivotRatio, &lpPivotBudget, &lpIntolerance, &addGomoryCuts,
      &addLearnedCuts, &gomoryCutLimit, &maxCutCos,    &slackdiv,      &weakenFull,    &weakenNonImplying,
      &bumpOnlyFalse,  &bumpCanceling,  &bumpLits,     &bitsOverflow,  &bitsReduced,   &bitsLearned,
      &bitsInput,      &configGJE,
  };
  std::unordered_map<std::string, Option*> name2opt;

  Options() {
    for (Option* opt : options) name2opt[opt->name] = opt;
  }

  std::string formulaName;

  void parseCommandLine(int argc, char** argv) {
    std::unordered_map<std::string, std::string> opt_val;
    for (int i = 1; i < argc; i++) {
      std::string arg = argv[i];
      if (arg.substr(0, 2) != "--") {
        formulaName = arg;
      } else {
        size_t eqindex = std::min(arg.size(), arg.find('='));
        std::string inputopt = arg.substr(2, eqindex - 2);
        if (name2opt.count(inputopt) == 0) {
          quit::exit_ERROR({"Unknown option: ", inputopt, ".\nCheck usage with ", argv[0], " --help"});
        } else {
          name2opt[inputopt]->parse(arg.substr(std::min(arg.size(), eqindex + 1)));
        }
      }
    }

    if (help) {
      usage(argv[0]);
      exit(0);
    }
  }

  void usage(char* name) {
    printf("Usage: %s [OPTIONS] instance.(opb|cnf|wcnf)\n", name);
    printf("or     cat instance.(opb|cnf|wcnf) | %s [OPTIONS]\n", name);
    printf("\n");
    printf("Options:\n");
    for (Option* opt : options) opt->printUsage(14);
  }
};

}  // namespace rs
