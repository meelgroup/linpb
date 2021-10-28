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
#include "globals.hpp"
#include "typedefs.hpp"

namespace rs {

enum class WatchStatus { DROPWATCH, KEEPWATCH, CONFLICTING };

class Solver;
struct Constr {  // internal solver constraint optimized for fast propagation
  virtual size_t getMemSize() const = 0;

  ID id;
  // NOTE: above attributes not strictly needed in cache-sensitive Constr, but it did not matter after testing
  struct {
    unsigned unused : 1;
    unsigned origin : 4;
    unsigned lbd : 27;
    unsigned markedfordel : 1;
    unsigned locked : 1;
    unsigned size : 30;
  } header;
  ActValC act;

  Constr(ID i, Origin o, bool lkd, unsigned int lngth) : id(i), act(0) {
    header = {0, (unsigned int)o, 0x07FFFFFF, 0, lkd, lngth};
  }
  virtual ~Constr() {}
  virtual void freeUp() = 0;  // poor man's destructor

  unsigned int size() const { return header.size; }
  void setLocked(bool lkd) { header.locked = lkd; }
  bool isLocked() { return header.locked; }
  Origin getOrigin() const { return (Origin)header.origin; }
  void setLBD(unsigned int lbd) { header.lbd = lbd; }
  unsigned int lbd() const { return header.lbd; }
  bool isMarkedForDelete() const { return header.markedfordel; }
  void markForDel() { header.markedfordel = 1; }

  virtual BigVal degree() const = 0;  // TODO: remove direct uses of these bigint methods, convert to ConstrExp instead
  virtual BigCoef coef(unsigned int i) const = 0;
  BigCoef largestCoef() const { return coef(0); };
  virtual Lit lit(unsigned int i) const = 0;

  virtual void initializeWatches(CRef cr, Solver& solver) = 0;
  virtual WatchStatus checkForPropagation(CRef cr, int& idx, Lit p, Solver& slvr) = 0;
  virtual void undoFalsified(int i) = 0;
  virtual void resolveWith(CeSuper confl, Lit l, IntSet* actSet, Solver& solver) = 0;

  virtual CeSuper toExpanded(ConstrExpPools& cePools) const = 0;

  std::ostream& operator<<(std::ostream& o) {
    for (size_t i = 0; i < size(); ++i) {
      o << coef(i) << "x" << lit(i) << " ";
    }
    o << ">= " << degree() << "\n";
    return o;
  }
  void print(const Solver& solver);

  bool isCorrectlyConflicting(const Solver& solver);
  bool isCorrectlyPropagating(const Solver& solver, int idx);
};

struct Clause final : public Constr {
  Lit data[];

  static size_t getMemSize(unsigned int length) { return (sizeof(Clause) + sizeof(Lit) * length) / sizeof(uint32_t); }
  size_t getMemSize() const { return getMemSize(size()); }

  BigVal degree() const { return 1; }
  BigCoef coef([[maybe_unused]] unsigned int i) const { return 1; }
  Lit lit(unsigned int i) const { return data[i]; }

  template <typename SMALL, typename LARGE>
  Clause(const ConstrExp<SMALL, LARGE>* constraint, bool locked, ID _id)
      : Constr(_id, constraint->orig, locked, constraint->vars.size()) {
    assert(_id > ID_Trivial);
    assert(constraint->vars.size() < INF);
    assert(constraint->getDegree() == 1);
    const unsigned int length = constraint->vars.size();

    for (unsigned int i = 0; i < length; ++i) {
      Var v = constraint->vars[i];
      assert(constraint->getLit(v) != 0);
      data[i] = constraint->getLit(v);
    }
  }
  void freeUp() {}

  void initializeWatches(CRef cr, Solver& solver);
  WatchStatus checkForPropagation(CRef cr, int& idx, Lit p, Solver& solver);
  void undoFalsified([[maybe_unused]] int i) { assert(false); }
  void resolveWith(CeSuper confl, Lit l, IntSet* actSet, Solver& solver);

  CeSuper toExpanded(ConstrExpPools& cePools) const;
};

struct Cardinality final : public Constr {
  unsigned int watchIdx;
  unsigned int degr;
  long long ntrailpops;
  Lit data[];

  static size_t getMemSize(unsigned int length) {
    return (sizeof(Cardinality) + sizeof(Lit) * length) / sizeof(uint32_t);
  }
  size_t getMemSize() const { return getMemSize(size()); }

  BigVal degree() const { return degr; }
  BigCoef coef([[maybe_unused]] unsigned int i) const { return 1; }
  Lit lit(unsigned int i) const { return data[i]; }

  template <typename SMALL, typename LARGE>
  Cardinality(const ConstrExp<SMALL, LARGE>* constraint, bool locked, ID _id)
      : Constr(_id, constraint->orig, locked, constraint->vars.size()),
        watchIdx(0),
        degr(static_cast<unsigned int>(constraint->getDegree())),
        ntrailpops(-1) {
    assert(_id > ID_Trivial);
    assert(constraint->vars.size() < INF);
    assert(aux::abs(constraint->coefs[constraint->vars[0]]) == 1);
    assert(constraint->getDegree() <= (LARGE)constraint->vars.size());
    const unsigned int length = constraint->vars.size();

    for (unsigned int i = 0; i < length; ++i) {
      Var v = constraint->vars[i];
      assert(constraint->getLit(v) != 0);
      data[i] = constraint->getLit(v);
    }
  }
  void freeUp() {}

  void initializeWatches(CRef cr, Solver& solver);
  WatchStatus checkForPropagation(CRef cr, int& idx, Lit p, Solver& solver);
  void undoFalsified([[maybe_unused]] int i) { assert(false); }
  void resolveWith(CeSuper confl, Lit l, IntSet* actSet, Solver& solver);

  CeSuper toExpanded(ConstrExpPools& cePools) const;
};

template <typename CF, typename DG>
struct Counting final : public Constr {
  unsigned int watchIdx;
  long long ntrailpops;
  DG degr;
  DG slack;
  Term<CF> data[];

  static size_t getMemSize(unsigned int length) {
    return (sizeof(Counting<CF, DG>) + sizeof(Term<CF>) * length) / sizeof(uint32_t);
  }
  size_t getMemSize() const { return getMemSize(size()); }

  BigVal degree() const { return degr; }
  BigCoef coef(unsigned int i) const { return data[i].c; }
  Lit lit(unsigned int i) const { return data[i].l; }

  template <typename SMALL, typename LARGE>
  Counting(const ConstrExp<SMALL, LARGE>* constraint, bool locked, ID _id)
      : Constr(_id, constraint->orig, locked, constraint->vars.size()),
        watchIdx(0),
        ntrailpops(-1),
        degr(static_cast<DG>(constraint->getDegree())),
        slack(0) {
    assert(_id > ID_Trivial);
    assert(aux::fitsIn<DG>(constraint->getDegree()));
    assert(aux::fitsIn<CF>(constraint->getLargestCoef()));
    ++stats.NCOUNTING;
    const unsigned int length = constraint->vars.size();

    for (unsigned int i = 0; i < length; ++i) {
      Var v = constraint->vars[i];
      assert(constraint->getLit(v) != 0);
      data[i] = {static_cast<CF>(aux::abs(constraint->coefs[v])), constraint->getLit(v)};
    }
  }
  void freeUp() {}

  void initializeWatches(CRef cr, Solver& solver);
  WatchStatus checkForPropagation(CRef cr, int& idx, Lit p, Solver& solver);
  void undoFalsified(int i);
  void resolveWith(CeSuper confl, Lit l, IntSet* actSet, Solver& solver);

  CePtr<ConstrExp<CF, DG>> expandTo(ConstrExpPools& cePools) const;
  CeSuper toExpanded(ConstrExpPools& cePools) const;

  bool hasCorrectSlack(const Solver& solver);
};

template <typename CF, typename DG>
struct Watched final : public Constr {
  unsigned int watchIdx;
  long long ntrailpops;
  DG degr;
  DG watchslack;
  Term<CF> data[];

  static size_t getMemSize(unsigned int length) {
    return (sizeof(Watched<CF, DG>) + sizeof(Term<CF>) * length) / sizeof(uint32_t);
  }
  size_t getMemSize() const { return getMemSize(size()); }

  BigVal degree() const { return degr; }
  BigCoef coef(unsigned int i) const { return aux::abs(data[i].c); }
  Lit lit(unsigned int i) const { return data[i].l; }

  template <typename SMALL, typename LARGE>
  Watched(const ConstrExp<SMALL, LARGE>* constraint, bool locked, ID _id)
      : Constr(_id, constraint->orig, locked, constraint->vars.size()),
        watchIdx(0),
        ntrailpops(-1),
        degr(static_cast<DG>(constraint->getDegree())),
        watchslack(0) {
    assert(_id > ID_Trivial);
    assert(aux::fitsIn<DG>(constraint->getDegree()));
    assert(aux::fitsIn<CF>(constraint->getLargestCoef()));
    ++stats.NWATCHED;
    const unsigned int length = constraint->vars.size();

    for (unsigned int i = 0; i < length; ++i) {
      Var v = constraint->vars[i];
      assert(constraint->getLit(v) != 0);
      data[i] = {static_cast<CF>(aux::abs(constraint->coefs[v])), constraint->getLit(v)};
    }
  }
  void freeUp() {}

  void initializeWatches(CRef cr, Solver& solver);
  WatchStatus checkForPropagation(CRef cr, int& idx, Lit p, Solver& solver);
  void undoFalsified(int i);
  void resolveWith(CeSuper confl, Lit l, IntSet* actSet, Solver& solver);

  CePtr<ConstrExp<CF, DG>> expandTo(ConstrExpPools& cePools) const;
  CeSuper toExpanded(ConstrExpPools& cePools) const;

  bool hasCorrectSlack(const Solver& solver);
  bool hasCorrectWatches(const Solver& solver);
};

template <typename CF, typename DG>
struct CountingSafe final : public Constr {
  unsigned int watchIdx;
  long long ntrailpops;
  DG* degr;
  DG* slack;
  Term<CF>* terms;  // array

  static size_t getMemSize([[maybe_unused]] unsigned int length) {
    return sizeof(CountingSafe<CF, DG>) / sizeof(uint32_t);
  }
  size_t getMemSize() const { return getMemSize(size()); }

  BigVal degree() const { return BigVal(*degr); }
  BigCoef coef(unsigned int i) const { return BigCoef(terms[i].c); }
  Lit lit(unsigned int i) const { return terms[i].l; }

  template <typename SMALL, typename LARGE>
  CountingSafe(const ConstrExp<SMALL, LARGE>* constraint, bool locked, ID _id)
      : Constr(_id, constraint->orig, locked, constraint->vars.size()),
        watchIdx(0),
        ntrailpops(-1),
        degr(new DG(static_cast<DG>(constraint->getDegree()))),
        slack(new DG(0)),
        terms(new Term<CF>[constraint->vars.size()]) {
    assert(_id > ID_Trivial);
    assert(aux::fitsIn<DG>(constraint->getDegree()));
    assert(aux::fitsIn<CF>(constraint->getLargestCoef()));
    ++stats.NCOUNTING;
    const unsigned int length = constraint->vars.size();

    for (unsigned int i = 0; i < length; ++i) {
      Var v = constraint->vars[i];
      assert(constraint->getLit(v) != 0);
      terms[i] = {static_cast<CF>(aux::abs(constraint->coefs[v])), constraint->getLit(v)};
    }
  }
  void freeUp() {
    delete degr;
    delete slack;
    delete[] terms;
  }

  void initializeWatches(CRef cr, Solver& solver);
  WatchStatus checkForPropagation(CRef cr, int& idx, Lit p, Solver& solver);
  void undoFalsified(int i);
  void resolveWith(CeSuper confl, Lit l, IntSet* actSet, Solver& solver);

  CePtr<ConstrExp<CF, DG>> expandTo(ConstrExpPools& cePools) const;
  CeSuper toExpanded(ConstrExpPools& cePools) const;

  bool hasCorrectSlack(const Solver& solver);
};

template <typename CF, typename DG>
struct WatchedSafe final : public Constr {
  unsigned int watchIdx;
  long long ntrailpops;
  DG* degr;
  DG* watchslack;
  Term<CF>* terms;  // array

  static size_t getMemSize([[maybe_unused]] unsigned int length) {
    return sizeof(WatchedSafe<CF, DG>) / sizeof(uint32_t);
  }
  size_t getMemSize() const { return getMemSize(size()); }

  BigVal degree() const { return BigVal(*degr); }
  BigCoef coef(unsigned int i) const { return BigCoef(aux::abs(terms[i].c)); }
  Lit lit(unsigned int i) const { return terms[i].l; }

  template <typename SMALL, typename LARGE>
  WatchedSafe(const ConstrExp<SMALL, LARGE>* constraint, bool locked, ID _id)
      : Constr(_id, constraint->orig, locked, constraint->vars.size()),
        watchIdx(0),
        ntrailpops(-1),
        degr(new DG(static_cast<DG>(constraint->getDegree()))),
        watchslack(new DG(0)),
        terms(new Term<CF>[constraint->vars.size()]) {
    assert(_id > ID_Trivial);
    assert(aux::fitsIn<DG>(constraint->getDegree()));
    assert(aux::fitsIn<CF>(constraint->getLargestCoef()));
    ++stats.NWATCHED;
    const unsigned int length = constraint->vars.size();

    for (unsigned int i = 0; i < length; ++i) {
      Var v = constraint->vars[i];
      assert(constraint->getLit(v) != 0);
      terms[i] = {static_cast<CF>(aux::abs(constraint->coefs[v])), constraint->getLit(v)};
    }
  }
  void freeUp() {
    delete degr;
    delete watchslack;
    delete[] terms;
  }

  void initializeWatches(CRef cr, Solver& solver);
  WatchStatus checkForPropagation(CRef cr, int& idx, Lit p, Solver& solver);
  void undoFalsified(int i);
  void resolveWith(CeSuper confl, Lit l, IntSet* actSet, Solver& solver);

  CePtr<ConstrExp<CF, DG>> expandTo(ConstrExpPools& cePools) const;
  CeSuper toExpanded(ConstrExpPools& cePools) const;

  bool hasCorrectSlack(const Solver& solver);
  bool hasCorrectWatches(const Solver& solver);
};

}  // namespace rs
