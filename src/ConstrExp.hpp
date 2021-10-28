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

#include <memory>
#include <sstream>
#include "ConstrSimple.hpp"
#include "IntSet.hpp"
#include "Logger.hpp"
#include "SolverStructs.hpp"
#include "globals.hpp"
#include "typedefs.hpp"

namespace rs {

enum class AssertionStatus { NONASSERTING, ASSERTING, FALSIFIED };

// shared_ptr-like wrapper around ConstrExp, ensuring it gets released back to the pool when no longer needed.
template <typename CE>
struct CePtr {
  CE* ce;

  // default constructor
  CePtr() : ce(nullptr) {}
  // regular constructor
  CePtr(CE* c) : ce(c) {
    if (ce) ce->increaseUsage();
  }
  // copy constructor
  CePtr(const CePtr<CE>& other) : ce{other.ce} {
    if (ce) ce->increaseUsage();
  }
  // copy constructor allowing for polymorphism
  template <typename T, typename = std::enable_if_t<std::is_convertible_v<T&, CE&>>>
  CePtr(const CePtr<T>& other) : ce{other.ce} {
    if (ce) ce->increaseUsage();
  }
  // move constructor
  CePtr(CePtr<CE>&& other) : ce{other.ce} { other.ce = nullptr; }
  // move constructor allowing for polymorphism
  template <typename T, typename = std::enable_if_t<std::is_convertible_v<T&, CE&>>>
  CePtr(CePtr<T>&& other) : ce{other.ce} {
    other.ce = nullptr;
  }
  // destructor
  ~CePtr() {
    if (ce) ce->decreaseUsage();
  }
  // assignment operator
  CePtr<CE>& operator=(const CePtr<CE>& other) {
    if (this == &other) return *this;
    if (ce) ce->decreaseUsage();
    ce = other.ce;
    if (ce) ce->increaseUsage();
    return *this;
  }
  // move assignment operator
  CePtr<CE>& operator=(CePtr<CE>&& other) {
    if (this == &other) return *this;
    if (ce) ce->decreaseUsage();
    ce = other.ce;
    other.ce = nullptr;
    return *this;
  }

  CE& operator*() const { return *ce; }
  CE* operator->() const { return ce; }
  explicit operator bool() const { return ce; }
};

struct ConstraintAllocator;
class ConstrExpPools;
class Solver;

struct ConstrExpSuper {
  std::vector<Var> vars;
  Origin orig = Origin::UNKNOWN;

  virtual ~ConstrExpSuper() {}

  virtual void increaseUsage() = 0;
  virtual void decreaseUsage() = 0;

  virtual void copyTo(Ce32 ce) const = 0;
  virtual void copyTo(Ce64 ce) const = 0;
  virtual void copyTo(Ce96 ce) const = 0;
  virtual void copyTo(Ce128 ce) const = 0;
  virtual void copyTo(CeArb ce) const = 0;

  virtual CeSuper clone(ConstrExpPools& ce) const = 0;
  virtual CRef toConstr(ConstraintAllocator& ca, bool locked, ID id) const = 0;
  virtual std::unique_ptr<ConstrSimpleSuper> toSimple() const = 0;

  virtual void resize(size_t s) = 0;
  virtual void resetBuffer(ID proofID = ID_Trivial) = 0;
  virtual void initializeLogging(std::shared_ptr<Logger>& l) = 0;
  virtual void stopLogging() = 0;
  virtual bool isReset() const = 0;
  virtual void reset() = 0;

  virtual Lit getLit(Lit l) const = 0;
  virtual bool hasLit(Lit l) const = 0;

  virtual void weaken(Var v) = 0;
  virtual void weakenLast() = 0;

  virtual bool hasNegativeSlack(const IntVecIt& level) const = 0;
  virtual bool hasNegativeSlack(const IntSet& assumptions) const = 0;
  virtual bool isTautology() const = 0;
  virtual bool isInconsistency() const = 0;
  virtual bool isSatisfied(const IntVecIt& level) const = 0;

  virtual void removeUnitsAndZeroes(const IntVecIt& level, const std::vector<int>& pos, bool doSaturation = true) = 0;
  virtual bool hasNoUnits(const IntVecIt& level) const = 0;
  virtual void removeZeroes() = 0;
  virtual bool hasNoZeroes() const = 0;

  virtual void saturate(const std::vector<Var>& vs, bool check = true) = 0;
  virtual void saturate(bool check = true) = 0;
  virtual bool isSaturated() const = 0;
  virtual void saturateAndFixOverflow(const IntVecIt& level, bool fullWeakening, int bitOverflow, int bitReduce,
                                      Lit asserting) = 0;
  virtual void saturateAndFixOverflowRational(const std::vector<double>& lpSolution) = 0;
  virtual bool fitsInDouble() const = 0;
  virtual bool largestCoefFitsIn(int bits) const = 0;

  virtual void weakenDivideRound(const IntVecIt& level, Lit l, bool slackdiv, bool fullWeakening) = 0;

  virtual bool divideByGCD() = 0;
  virtual void postProcess(const IntVecIt& level, const std::vector<int>& pos, bool sortFirst, Stats& sts) = 0;
  virtual AssertionStatus isAssertingBefore(const IntVecIt& level, int lvl) const = 0;
  virtual int getAssertionLevel(const IntVecIt& level, const std::vector<int>& pos) const = 0;
  virtual void heuristicWeakening(const IntVecIt& level, const std::vector<int>& pos, Stats& sts) = 0;

  virtual bool simplifyToCardinality(bool equivalencePreserving) = 0;
  virtual void simplifyToClause() = 0;
  virtual bool isCardinality() const = 0;
  virtual int getCardinalityDegree() const = 0;
  virtual int getCardinalityDegreeWithZeroes() = 0;
  virtual void simplifyToMinLengthCardinality() = 0;
  virtual bool isClause() const = 0;
  virtual void sortInDecreasingCoefOrder(const std::function<bool(Var, Var)>& tiebreaker = [](Var, Var) {
    return false;
  }) = 0;
  virtual bool isSortedInDecreasingCoefOrder() const = 0;
  virtual void sort(const std::function<bool(Var, Var)>& comp) = 0;

  virtual ID logAsInput() = 0;
  virtual ID logProofLine() = 0;
  virtual ID logProofLineWithInfo([[maybe_unused]] std::string&& info, [[maybe_unused]] const Stats& sts) = 0;
  virtual void logUnit(const IntVecIt& level, const std::vector<int>& pos, Var v_unit, const Stats& sts) = 0;
  virtual void logInconsistency(const IntVecIt& level, const std::vector<int>& pos, const Stats& sts) = 0;

  virtual void toStreamAsOPB(std::ostream& o) const = 0;
  virtual void toStreamWithAssignment(std::ostream& o, const IntVecIt& level, const std::vector<int>& pos) const = 0;

  virtual void resolveWith(Ce32 c, Lit l, IntSet* actSet, const IntVecIt& Level, const std::vector<int>& Pos) = 0;
  virtual void resolveWith(Ce64 c, Lit l, IntSet* actSet, const IntVecIt& Level, const std::vector<int>& Pos) = 0;
  virtual void resolveWith(Ce96 c, Lit l, IntSet* actSet, const IntVecIt& Level, const std::vector<int>& Pos) = 0;
  virtual void resolveWith(Ce128 c, Lit l, IntSet* actSet, const IntVecIt& Level, const std::vector<int>& Pos) = 0;
  virtual void resolveWith(CeArb c, Lit l, IntSet* actSet, const IntVecIt& Level, const std::vector<int>& Pos) = 0;
  virtual void resolveWith(const Clause& c, Lit l, IntSet* actSet, const IntVecIt& Level,
                           const std::vector<int>& Pos) = 0;
  virtual void resolveWith(const Cardinality& c, Lit l, IntSet* actSet, const IntVecIt& Level,
                           const std::vector<int>& Pos) = 0;
};

template <typename CE>
class ConstrExpPool;

template <typename SMALL, typename LARGE>  // LARGE should be able to fit sums of SMALL
struct ConstrExp final : public ConstrExpSuper {
 private:
  ConstrExpPool<ConstrExp<SMALL, LARGE>>& pool;
  long long usageCount = 0;

 public:
  LARGE degree = 0;
  LARGE rhs = 0;
  std::vector<SMALL> coefs;
  std::vector<bool> used;
  std::stringstream proofBuffer;
  std::shared_ptr<Logger> plogger;

 private:
  void remove(Var v);
  bool increasesSlack(const IntVecIt& level, Var v) const;
  LARGE calcDegree() const;
  LARGE calcRhs() const;
  bool testConstraint() const;
  bool falsified(const IntVecIt& level, Var v) const;
  template <typename T>
  std::string proofMult(const T& mult) {
    std::stringstream ss;
    if (mult != 1) ss << mult << " * ";
    return ss.str();
  }
  void logIfUnit(Lit l, const SMALL& c, const IntVecIt& level, const std::vector<int>& pos);

 public:
  ConstrExp(ConstrExpPool<ConstrExp<SMALL, LARGE>>& cep);
  void increaseUsage() {
    ++usageCount;
    assert(usageCount > 0);
  }
  void decreaseUsage() {
    assert(usageCount > 0);
    if (--usageCount == 0) {
      pool.release(this);
    }
  }

  void copyTo(Ce32 ce) const { copyTo_(ce); }
  void copyTo(Ce64 ce) const { copyTo_(ce); }
  void copyTo(Ce96 ce) const { copyTo_(ce); }
  void copyTo(Ce128 ce) const { copyTo_(ce); }
  void copyTo(CeArb ce) const { copyTo_(ce); }

  CeSuper clone(ConstrExpPools& ce) const;
  CRef toConstr(ConstraintAllocator& ca, bool locked, ID id) const;
  std::unique_ptr<ConstrSimpleSuper> toSimple() const;

  void resize(size_t s);
  void resetBuffer(ID proofID = ID_Trivial);
  void initializeLogging(std::shared_ptr<Logger>& l);
  void stopLogging();

  bool isReset() const;
  void reset();

  LARGE getRhs() const;
  LARGE getDegree() const;
  SMALL getCoef(Lit l) const;
  SMALL getLargestCoef() const;
  SMALL getSmallestCoef() const;
  LARGE getCutoffVal() const;
  Lit getLit(Lit l) const;
  bool hasLit(Lit l) const;

  void addRhs(const LARGE& r);
  void addLhs(const SMALL& cf, Lit l);  // TODO: Term?
  void weaken(const SMALL& m, Var v);
  void weaken(Var v);
  void weakenLast();

  LARGE getSlack(const IntVecIt& level) const;
  bool hasNegativeSlack(const IntVecIt& level) const;
  LARGE getSlack(const IntSet& assumptions) const;
  bool hasNegativeSlack(const IntSet& assumptions) const;
  bool isTautology() const;
  bool isInconsistency() const;
  bool isSatisfied(const IntVecIt& level) const;

  // @post: preserves order of vars
  void removeUnitsAndZeroes(const IntVecIt& level, const std::vector<int>& pos, bool doSaturation = true);
  bool hasNoUnits(const IntVecIt& level) const;
  // @post: mutates order of vars
  void removeZeroes();
  bool hasNoZeroes() const;

  // @post: preserves order of vars
  void saturate(const std::vector<Var>& vs, bool check = true);
  void saturate(bool check = true);
  bool isSaturated() const;
  /*
   * Fixes overflow
   * @post: saturated
   * @post: nothing else if bitOverflow == 0
   * @post: the largest coefficient is less than 2^bitOverflow
   * @post: the degree and rhs are less than 2^bitOverflow * INF
   * @post: if overflow happened, all division until 2^bitReduce happened
   * @post: the constraint remains conflicting or propagating on asserting
   */
  void saturateAndFixOverflow(const IntVecIt& level, bool fullWeakening, int bitOverflow, int bitReduce, Lit asserting);
  /*
   * Fixes overflow for rationals
   * @post: saturated
   * @post: none of the coefficients, degree, or rhs exceed INFLPINT
   */
  void saturateAndFixOverflowRational(const std::vector<double>& lpSolution);
  bool fitsInDouble() const;
  bool largestCoefFitsIn(int bits) const;

  template <typename S, typename L>
  void addUp(CePtr<ConstrExp<S, L>> c, const SMALL& cmult = 1, const SMALL& thismult = 1) {
    assert(cmult >= 1);
    assert(thismult >= 1);
    if (plogger) proofBuffer << proofMult(thismult) << c->proofBuffer.str() << proofMult(cmult) << "+ ";
    if (thismult > 1) {
      degree *= thismult;
      rhs *= thismult;
      for (Var v : vars) coefs[v] *= thismult;
    }
    rhs += static_cast<LARGE>(cmult) * static_cast<LARGE>(c->rhs);
    degree += static_cast<LARGE>(cmult) * static_cast<LARGE>(c->degree);
    for (Var v : c->vars) {
      assert(v < (Var)coefs.size());
      assert(v > 0);
      SMALL val = cmult * static_cast<SMALL>(c->coefs[v]);
      if (!used[v]) {
        assert(coefs[v] == 0);
        vars.push_back(v);
        coefs[v] = val;
        used[v] = true;
      } else {
        if ((coefs[v] < 0) != (val < 0)) degree -= std::min(aux::abs(coefs[v]), aux::abs(val));
        coefs[v] += val;
      }
    }
  }

  void invert();
  void multiply(const SMALL& m);
  void divide(const LARGE& d);
  void divideRoundUp(const LARGE& d);
  void weakenDivideRound(const IntVecIt& level, Lit l, bool slackdiv, bool fullWeakening);
  void weakenNonDivisibleNonFalsifieds(const IntVecIt& level, const LARGE& div, bool fullWeakening, Lit asserting);
  void applyMIR(const LARGE& d, std::function<Lit(Var)> toLit);

  bool divideByGCD();
  // NOTE: only equivalence preserving operations!
  void postProcess(const IntVecIt& level, const std::vector<int>& pos, bool sortFirst, Stats& sts);
  AssertionStatus isAssertingBefore(const IntVecIt& level, int lvl) const;
  // @return: earliest decision level that propagates a variable
  int getAssertionLevel(const IntVecIt& level, const std::vector<int>& pos) const;
  // @post: preserves order after removeZeroes()
  void weakenNonImplied(const IntVecIt& level, const LARGE& slack, Stats& sts);
  // @post: preserves order after removeZeroes()
  bool weakenNonImplying(const IntVecIt& level, const SMALL& propCoef, const LARGE& slack, Stats& sts);
  // @post: preserves order after removeZeroes()
  void heuristicWeakening(const IntVecIt& level, const std::vector<int>& pos, Stats& sts);

  // @post: preserves order
  template <typename T>
  void weakenSmalls(const T& limit) {
    for (Var v : vars)
      if (aux::abs(coefs[v]) < limit) weaken(v);
    saturate();
  }

  LARGE absCoeffSum() const;

  // @post: preserves order of vars
  bool simplifyToCardinality(bool equivalencePreserving);
  bool isCardinality() const;
  int getCardinalityDegree() const;
  int getCardinalityDegreeWithZeroes();
  void simplifyToMinLengthCardinality();
  void simplifyToClause();
  bool isClause() const;

  void sortInDecreasingCoefOrder(const std::function<bool(Var, Var)>& tiebreaker = [](Var, Var) { return false; });
  bool isSortedInDecreasingCoefOrder() const;
  void sort(const std::function<bool(Var, Var)>& comp);

  ID logAsInput();
  ID logProofLine();
  ID logProofLineWithInfo([[maybe_unused]] std::string&& info, [[maybe_unused]] const Stats& sts);
  // @pre: reducible to unit over v
  void logUnit(const IntVecIt& level, const std::vector<int>& pos, Var v_unit, const Stats& sts);
  void logInconsistency(const IntVecIt& level, const std::vector<int>& pos, const Stats& sts);

  void toStreamAsOPB(std::ostream& o) const;
  void toStreamWithAssignment(std::ostream& o, const IntVecIt& level, const std::vector<int>& pos) const;

  void resolveWith(Ce32 c, Lit l, IntSet* actSet, const IntVecIt& Level, const std::vector<int>& Pos);
  void resolveWith(Ce64 c, Lit l, IntSet* actSet, const IntVecIt& Level, const std::vector<int>& Pos);
  void resolveWith(Ce96 c, Lit l, IntSet* actSet, const IntVecIt& Level, const std::vector<int>& Pos);
  void resolveWith(Ce128 c, Lit l, IntSet* actSet, const IntVecIt& Level, const std::vector<int>& Pos);
  void resolveWith(CeArb c, Lit l, IntSet* actSet, const IntVecIt& Level, const std::vector<int>& Pos);
  void resolveWith(const Clause& c, Lit l, IntSet* actSet, const IntVecIt& Level, const std::vector<int>& Pos);
  void resolveWith(const Cardinality& c, Lit l, IntSet* actSet, const IntVecIt& Level, const std::vector<int>& Pos);

 private:
  template <typename CF, typename DG>
  void genericResolve(CePtr<ConstrExp<CF, DG>> reason, Lit l, IntSet* actSet, const IntVecIt& Level,
                      const std::vector<int>& Pos) {
    assert(getCoef(-l) > 0);
    stats.NADDEDLITERALS += reason->vars.size();
    reason->removeUnitsAndZeroes(Level, Pos);
    if (options.weakenNonImplying) reason->weakenNonImplying(Level, reason->getCoef(l), reason->getSlack(Level), stats);
    reason->saturateAndFixOverflow(Level, (bool)options.weakenFull, options.bitsOverflow.get(),
                                   options.bitsReduced.get(), l);
    assert(reason->getCoef(l) > 0);
    assert(reason->getCoef(l) > reason->getSlack(Level));
    reason->weakenDivideRound(Level, l, (bool)options.slackdiv, (bool)options.weakenFull);
    assert(reason->getCoef(l) > 0);
    assert(reason->getSlack(Level) <= 0);
    if (actSet != nullptr) {
      for (Var v : reason->vars) {
        Lit l = reason->getLit(v);
        if (options.bumpLits) {
          actSet->add(l);
        } else {
          if (!options.bumpOnlyFalse || isFalse(Level, l)) actSet->add(v);
          if (options.bumpCanceling && getLit(v) == -l) actSet->add(-v);
        }
      }
    }
    SMALL reason_coef_l = static_cast<SMALL>(reason->getCoef(l));  // NOTE: SMALL >= CF
    SMALL confl_coef_l = getCoef(-l);
    SMALL gcd_coef_l = aux::gcd(reason_coef_l, confl_coef_l);
    addUp(reason, confl_coef_l / gcd_coef_l, reason_coef_l / gcd_coef_l);
    saturateAndFixOverflow(Level, (bool)options.weakenFull, options.bitsOverflow.get(), options.bitsReduced.get(), 0);
    assert(getCoef(-l) == 0);
    assert(hasNegativeSlack(Level));
  }

  template <typename S, typename L>
  void copyTo_(CePtr<ConstrExp<S, L>> out) const {
    // TODO: assert whether S/L can fit SMALL/LARGE? Not always possible.
    assert(out->isReset());
    out->degree = static_cast<L>(degree);
    out->rhs = static_cast<L>(rhs);
    out->orig = orig;
    out->vars = vars;
    assert(out->coefs.size() == coefs.size());
    for (Var v : vars) {
      out->coefs[v] = static_cast<S>(coefs[v]);
      assert(used[v] == true);
      assert(out->used[v] == false);
      out->used[v] = true;
    }
    if (plogger) {
      out->proofBuffer.str(std::string());
      out->proofBuffer << proofBuffer.str();
    }
  }

  template <typename S, typename L>
  std::unique_ptr<ConstrSimple<S, L>> toSimple_() const {
    std::unique_ptr<ConstrSimple<S, L>> result = std::make_unique<ConstrSimple<S, L>>();
    result->rhs = static_cast<L>(rhs);
    result->terms.reserve(vars.size());
    for (Var v : vars)
      if (coefs[v] != 0) result->terms.emplace_back(static_cast<S>(coefs[v]), v);
    if (plogger) result->proofLine = proofBuffer.str();
    result->orig = orig;
    return result;
  }
};

template <typename S, typename L>
std::ostream& operator<<(std::ostream& o, const ConstrExp<S, L>& C) {
  std::vector<Var> vars = C.vars;
  std::sort(vars.begin(), vars.end(), [](Var v1, Var v2) { return v1 < v2; });
  for (Var v : vars) {
    Lit l = C.coefs[v] < 0 ? -v : v;
    o << C.getCoef(l) << "x" << l << " ";
  }
  o << ">= " << C.degree;
  return o;
}

template <typename CE>
class ConstrExpPool {  // TODO: private constructor for ConstrExp, only accessible to ConstrExpPool?
  size_t n = 0;
  std::vector<CE*> ces;
  std::vector<CE*> availables;
  std::shared_ptr<Logger> plogger;

 public:
  ~ConstrExpPool() {
    for (CE* ce : ces) delete ce;
  }

  // add reset
  void reset() {
    n = 0;
    for (CE* ce : ces) delete ce;
    ces.clear();
    availables.clear();
    plogger.reset();	
  }

  void resize(size_t newn) {
    assert(n <= INF);
    n = newn;
    for (CE* ce : ces) ce->resize(n);
  }

  void initializeLogging(std::shared_ptr<Logger>& lgr) {
    plogger = lgr;
    for (CE* ce : ces) ce->initializeLogging(lgr);
  }

  CE* take() {
    assert(ces.size() < 20);  // Sanity check that no large amounts of ConstrExps are created
    if (availables.size() == 0) {
      ces.emplace_back(new CE(*this));
      ces.back()->resize(n);
      ces.back()->initializeLogging(plogger);
      availables.push_back(ces.back());
    }
    CE* result = availables.back();
    availables.pop_back();
    assert(result->isReset());
    assert(result->coefs.size() == n);
    return result;
  }

  void release(CE* ce) {
    assert(std::any_of(ces.begin(), ces.end(), [&](CE* i) { return i == ce; }));
    assert(std::none_of(availables.begin(), availables.end(), [&](CE* i) { return i == ce; }));
    ce->reset();
    availables.push_back(ce);
  }
};

class ConstrExpPools {
  ConstrExpPool<ConstrExp32> ce32s;
  ConstrExpPool<ConstrExp64> ce64s;
  ConstrExpPool<ConstrExp96> ce96s;
  ConstrExpPool<ConstrExp128> ce128s;
  ConstrExpPool<ConstrExpArb> ceArbs;

 public:
  void resize(size_t newn);
  void initializeLogging(std::shared_ptr<Logger> lgr);

  template <typename SMALL, typename LARGE>
  CePtr<ConstrExp<SMALL, LARGE>> take();  // NOTE: only call specializations

  Ce32 take32();
  Ce64 take64();
  Ce96 take96();
  Ce128 take128();
  CeArb takeArb();

  // add reset
  void reset() {
    ce32s.reset();
    ce64s.reset();
    ce96s.reset();
    ce128s.reset();
    ceArbs.reset();
  }
};

}  // namespace rs
