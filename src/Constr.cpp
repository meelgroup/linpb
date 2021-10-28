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

#include "Constr.hpp"
#include "Solver.hpp"

namespace rs {

void Clause::initializeWatches(CRef cr, Solver& solver) {
  auto& Level = solver.Level;
  auto& adj = solver.adj;

  unsigned int length = size();
  assert(length >= 1);
  if (length == 1) {
    assert(solver.decisionLevel() == 0);
    assert(isCorrectlyPropagating(solver, 0));
    solver.propagate(data[0], cr);
    return;
  }

  unsigned int watch = 0;
  for (unsigned int i = 0; i < length && watch <= 1; ++i) {
    Lit l = data[i];
    if (!isFalse(Level, l)) {
      data[i] = data[watch];
      data[watch++] = l;
    }
  }
  assert(watch >= 1);  // we found enough watches to satisfy the constraint
  assert((watch == 1) == isFalse(Level, data[1]));
  if (watch == 1) {
    assert(!isFalse(Level, data[0]));
    if (!isTrue(Level, data[0])) {
      assert(isCorrectlyPropagating(solver, 0));
      solver.propagate(data[0], cr);
    }
    for (unsigned int i = 2; i < length; ++i) {  // ensure last watch is last falsified literal
      Lit l = data[i];
      assert(isFalse(Level, l));
      if (Level[-l] > Level[-data[1]]) {
        data[i] = data[1];
        data[1] = l;
      }
    }
  }
  for (unsigned int i = 0; i < 2; ++i) adj[data[i]].emplace_back(cr, data[1 - i] - INF);  // add blocked literal
}

WatchStatus Clause::checkForPropagation(CRef cr, int& idx, Lit p, Solver& solver) {
  auto& Level = solver.Level;
  auto& adj = solver.adj;

  assert(idx < 0);
  assert(p == data[0] || p == data[1]);
  assert(size() > 1);
  int widx = 0;
  Lit watch = data[0];
  Lit otherwatch = data[1];
  if (p == data[1]) {
    widx = 1;
    watch = data[1];
    otherwatch = data[0];
  }
  assert(p == watch);
  assert(p != otherwatch);
  if (isTrue(Level, otherwatch)) {
    idx = otherwatch - INF;         // set new blocked literal
    return WatchStatus::KEEPWATCH;  // constraint is satisfied
  }

  const unsigned int length = size();
  for (unsigned int i = 2; i < length; ++i) {
    Lit l = data[i];
    if (!isFalse(Level, l)) {
      unsigned int mid = i / 2 + 1;
      data[i] = data[mid];
      data[mid] = watch;
      data[widx] = l;
      adj[l].emplace_back(cr, otherwatch - INF);
      stats.NWATCHCHECKS += i - 1;
      return WatchStatus::DROPWATCH;
    }
  }
  stats.NWATCHCHECKS += length - 2;

  assert(isFalse(Level, watch));
  for (unsigned int i = 2; i < length; ++i) assert(isFalse(Level, data[i]));
  if (isFalse(Level, otherwatch)) {
    assert(isCorrectlyConflicting(solver));
    return WatchStatus::CONFLICTING;
  } else {
    assert(!isTrue(Level, otherwatch));
    ++stats.NPROPCLAUSE;
    assert(isCorrectlyPropagating(solver, 1 - widx));
    solver.propagate(otherwatch, cr);
  }
  ++stats.NPROPCHECKS;
  return WatchStatus::KEEPWATCH;
}

void Clause::resolveWith(CeSuper confl, Lit l, IntSet* actSet, Solver& solver) {
  confl->resolveWith(*this, l, actSet, solver.getLevel(), solver.getPos());
}

CeSuper Clause::toExpanded(ConstrExpPools& cePools) const {
  Ce32 result = cePools.take32();
  result->addRhs(1);
  for (size_t i = 0; i < size(); ++i) {
    result->addLhs(1, data[i]);
  }
  result->orig = getOrigin();
  if (result->plogger) result->resetBuffer(id);
  return result;
}

void Cardinality::initializeWatches(CRef cr, Solver& solver) {
  auto& Level = solver.Level;
  [[maybe_unused]] auto& Pos = solver.Pos;
  auto& adj = solver.adj;

  unsigned int length = size();
  if (degr >= length) {
    assert(solver.decisionLevel() == 0);
    for (unsigned int i = 0; i < length; ++i) {
      assert(isUnknown(Pos, data[i]));
      assert(isCorrectlyPropagating(solver, i));
      solver.propagate(data[i], cr);
    }
    return;
  }

  unsigned int watch = 0;
  for (unsigned int i = 0; i < length && watch <= degr; ++i) {
    Lit l = data[i];
    if (!isFalse(Level, l)) {
      data[i] = data[watch];
      data[watch++] = l;
    }
  }
  assert(watch >= degr);  // we found enough watches to satisfy the constraint
  if (isFalse(Level, data[degr])) {
    for (unsigned int i = 0; i < degr; ++i) {
      assert(!isFalse(Level, data[i]));
      if (!isTrue(Level, data[i])) {
        assert(isCorrectlyPropagating(solver, i));
        solver.propagate(data[i], cr);
      }
    }
    for (unsigned int i = degr + 1; i < size(); ++i) {  // ensure last watch is last falsified literal
      Lit l = data[i];
      assert(isFalse(Level, l));
      if (Level[-l] > Level[-data[degr]]) {
        data[i] = data[degr];
        data[degr] = l;
      }
    }
  }
  for (unsigned int i = 0; i <= degr; ++i) adj[data[i]].emplace_back(cr, i);  // add watch index
}

WatchStatus Cardinality::checkForPropagation(CRef cr, int& idx, [[maybe_unused]] Lit p, Solver& solver) {
  auto& Level = solver.Level;
  auto& adj = solver.adj;

  assert(idx >= 0);
  assert(idx < INF);
  assert(data[idx] == p);
  const unsigned int length = size();
  if (!options.propIdx || ntrailpops < stats.NTRAILPOPS) {
    ntrailpops = stats.NTRAILPOPS;
    watchIdx = degr + 1;
  }
  assert(watchIdx > degr);
  stats.NWATCHCHECKS -= watchIdx;
  for (; watchIdx < length; ++watchIdx) {
    Lit l = data[watchIdx];
    if (!isFalse(Level, l)) {
      unsigned int mid = (watchIdx + degr + 1) / 2;
      assert(mid <= watchIdx);
      assert(mid > degr);
      data[watchIdx] = data[mid];
      data[mid] = data[idx];
      data[idx] = l;
      adj[l].emplace_back(cr, idx);
      stats.NWATCHCHECKS += watchIdx + 1;
      return WatchStatus::DROPWATCH;
    }
  }
  stats.NWATCHCHECKS += watchIdx;
  assert(isFalse(Level, data[idx]));
  for (unsigned int i = degr + 1; i < length; ++i) assert(isFalse(Level, data[i]));
  for (int i = 0; i <= (int)degr; ++i)
    if (i != idx && isFalse(Level, data[i])) {
      assert(isCorrectlyConflicting(solver));
      return WatchStatus::CONFLICTING;
    }
  for (int i = 0; i <= (int)degr; ++i) {
    Lit l = data[i];
    if (i != idx && !isTrue(Level, l)) {
      ++stats.NPROPCARD;
      assert(isCorrectlyPropagating(solver, i));
      solver.propagate(l, cr);
    }
  }
  stats.NPROPCHECKS += degr + 1;
  return WatchStatus::KEEPWATCH;
}

void Cardinality::resolveWith(CeSuper confl, Lit l, IntSet* actSet, Solver& solver) {
  confl->resolveWith(*this, l, actSet, solver.getLevel(), solver.getPos());
}

CeSuper Cardinality::toExpanded(ConstrExpPools& cePools) const {
  Ce32 result = cePools.take32();
  result->addRhs(degr);
  for (size_t i = 0; i < size(); ++i) {
    result->addLhs(1, data[i]);
  }
  result->orig = getOrigin();
  if (result->plogger) result->resetBuffer(id);
  return result;
}

template <typename CF, typename DG>
void Counting<CF, DG>::initializeWatches(CRef cr, Solver& solver) {
  auto& Level = solver.Level;
  auto& Pos = solver.Pos;
  auto& adj = solver.adj;
  auto& qhead = solver.qhead;

  slack = -degr;
  unsigned int length = size();
  for (unsigned int i = 0; i < length; ++i) {
    Lit l = data[i].l;
    adj[l].emplace_back(cr, i + INF);
    if (!isFalse(Level, l) || Pos[toVar(l)] >= qhead) slack += data[i].c;
  }

  assert(slack >= 0);
  assert(hasCorrectSlack(solver));
  if (slack < data[0].c) {  // propagate
    for (unsigned int i = 0; i < length && data[i].c > slack; ++i)
      if (isUnknown(Pos, data[i].l)) {
        assert(isCorrectlyPropagating(solver, i));
        solver.propagate(data[i].l, cr);
      }
  }
}

template <typename CF, typename DG>
WatchStatus Counting<CF, DG>::checkForPropagation(CRef cr, int& idx, Lit p, Solver& solver) {
  _unused(p);
  auto& Pos = solver.Pos;

  assert(idx >= INF);
  assert(data[idx - INF].l == p);
  const unsigned int length = size();
  const CF& lrgstCf = data[0].c;
  const CF& c = data[idx - INF].c;

  slack -= c;
  assert(hasCorrectSlack(solver));

  if (slack < 0) {
    assert(isCorrectlyConflicting(solver));
    return WatchStatus::CONFLICTING;
  }
  if (slack < lrgstCf) {
    if (!options.propIdx || ntrailpops < stats.NTRAILPOPS) {
      ntrailpops = stats.NTRAILPOPS;
      watchIdx = 0;
    }
    stats.NPROPCHECKS -= watchIdx;
    for (; watchIdx < length && data[watchIdx].c > slack; ++watchIdx) {
      const Lit l = data[watchIdx].l;
      if (isUnknown(Pos, l)) {
        stats.NPROPCLAUSE += (degr == 1);
        stats.NPROPCARD += (degr != 1 && lrgstCf == 1);
        ++stats.NPROPCOUNTING;
        assert(isCorrectlyPropagating(solver, watchIdx));
        solver.propagate(l, cr);
      }
    }
    stats.NPROPCHECKS += watchIdx;
  }
  return WatchStatus::KEEPWATCH;
}

template <typename CF, typename DG>
void Counting<CF, DG>::undoFalsified(int i) {
  assert(i >= INF);
  slack += data[i - INF].c;
  ++stats.NWATCHLOOKUPSBJ;
}

template <typename CF, typename DG>
void Counting<CF, DG>::resolveWith(CeSuper confl, Lit l, IntSet* actSet, Solver& solver) {
  confl->resolveWith(expandTo(solver.cePools), l, actSet, solver.getLevel(), solver.getPos());
}

template <typename CF, typename DG>
CePtr<ConstrExp<CF, DG>> Counting<CF, DG>::expandTo(ConstrExpPools& cePools) const {
  CePtr<ConstrExp<CF, DG>> result = cePools.take<CF, DG>();
  result->addRhs(degr);
  for (size_t i = 0; i < size(); ++i) {
    result->addLhs(data[i].c, data[i].l);
  }
  result->orig = getOrigin();
  if (result->plogger) result->resetBuffer(id);
  return result;
}

template <typename CF, typename DG>
CeSuper Counting<CF, DG>::toExpanded(ConstrExpPools& cePools) const {
  return expandTo(cePools);
}

template <typename CF, typename DG>
void Watched<CF, DG>::initializeWatches(CRef cr, Solver& solver) {
  auto& Level = solver.Level;
  auto& Pos = solver.Pos;
  auto& adj = solver.adj;
  auto& qhead = solver.qhead;

  watchslack = -degr;
  unsigned int length = size();
  const CF lrgstCf = aux::abs(data[0].c);
  for (unsigned int i = 0; i < length && watchslack < lrgstCf; ++i) {
    Lit l = data[i].l;
    if (!isFalse(Level, l) || Pos[toVar(l)] >= qhead) {
      assert(data[i].c > 0);
      watchslack += data[i].c;
      data[i].c = -data[i].c;
      adj[l].emplace_back(cr, i + INF);
    }
  }
  assert(watchslack >= 0);
  assert(hasCorrectSlack(solver));
  if (watchslack < lrgstCf) {
    // set sufficient falsified watches
    std::vector<unsigned int> falsifiedIdcs;
    falsifiedIdcs.reserve(length);
    for (unsigned int i = 0; i < length; ++i)
      if (isFalse(Level, data[i].l) && Pos[toVar(data[i].l)] < qhead) falsifiedIdcs.push_back(i);
    std::sort(falsifiedIdcs.begin(), falsifiedIdcs.end(),
              [&](unsigned int i1, unsigned int i2) { return Pos[toVar(data[i1].l)] > Pos[toVar(data[i2].l)]; });
    DG diff = lrgstCf - watchslack;
    for (unsigned int i : falsifiedIdcs) {
      assert(data[i].c > 0);
      diff -= data[i].c;
      data[i].c = -data[i].c;
      adj[data[i].l].emplace_back(cr, i + INF);
      if (diff <= 0) break;
    }
    // perform initial propagation
    for (unsigned int i = 0; i < length && aux::abs(data[i].c) > watchslack; ++i)
      if (isUnknown(Pos, data[i].l)) {
        assert(isCorrectlyPropagating(solver, i));
        solver.propagate(data[i].l, cr);
      }
  }
}

template <typename CF, typename DG>
WatchStatus Watched<CF, DG>::checkForPropagation(CRef cr, int& idx, Lit p, Solver& solver) {
  _unused(p);
  auto& Level = solver.Level;
  auto& Pos = solver.Pos;
  auto& adj = solver.adj;

  assert(idx >= INF);
  assert(data[idx - INF].l == p);
  const unsigned int length = size();
  const CF lrgstCf = aux::abs(data[0].c);
  CF& c = data[idx - INF].c;

  if (!options.propIdx || ntrailpops < stats.NTRAILPOPS) {
    ntrailpops = stats.NTRAILPOPS;
    watchIdx = 0;
  }

  assert(c < 0);
  watchslack += c;
  if (!options.propSup ||
      watchslack - c >= lrgstCf) {  // look for new watches if previously, slack was at least lrgstCf
    stats.NWATCHCHECKS -= watchIdx;
    for (; watchIdx < length && watchslack < lrgstCf; ++watchIdx) {
      const CF& cf = data[watchIdx].c;
      const Lit l = data[watchIdx].l;
      if (cf > 0 && !isFalse(Level, l)) {
        watchslack += cf;
        data[watchIdx].c = -cf;
        adj[l].emplace_back(cr, watchIdx + INF);
      }
    }  // NOTE: first innermost loop of RoundingSat
    stats.NWATCHCHECKS += watchIdx;
    if (watchslack < lrgstCf) {
      assert(watchIdx == length);
      watchIdx = 0;
    }
  }  // else we did not find enough watches last time, so we can skip looking for them now

  assert(hasCorrectSlack(solver));
  assert(hasCorrectWatches(solver));

  if (watchslack >= lrgstCf) {
    c = -c;
    return WatchStatus::DROPWATCH;
  }
  if (watchslack < 0) {
    assert(isCorrectlyConflicting(solver));
    return WatchStatus::CONFLICTING;
  }
  // keep the watch, check for propagation
  stats.NPROPCHECKS -= watchIdx;
  for (; watchIdx < length && aux::abs(data[watchIdx].c) > watchslack; ++watchIdx) {
    const Lit l = data[watchIdx].l;
    if (isUnknown(Pos, l)) {
      stats.NPROPCLAUSE += (degr == 1);
      stats.NPROPCARD += (degr != 1 && lrgstCf == 1);
      ++stats.NPROPWATCH;
      assert(isCorrectlyPropagating(solver, watchIdx));
      solver.propagate(l, cr);
    }  // NOTE: second innermost loop of RoundingSat
  }
  stats.NPROPCHECKS += watchIdx;
  return WatchStatus::KEEPWATCH;
}

template <typename CF, typename DG>
void Watched<CF, DG>::undoFalsified(int i) {
  assert(i >= INF);
  assert(data[i - INF].c < 0);
  watchslack -= data[i - INF].c;
  ++stats.NWATCHLOOKUPSBJ;
}

template <typename CF, typename DG>
void Watched<CF, DG>::resolveWith(CeSuper confl, Lit l, IntSet* actSet, Solver& solver) {
  confl->resolveWith(expandTo(solver.cePools), l, actSet, solver.getLevel(), solver.getPos());
}

template <typename CF, typename DG>
CePtr<ConstrExp<CF, DG>> Watched<CF, DG>::expandTo(ConstrExpPools& cePools) const {
  CePtr<ConstrExp<CF, DG>> result = cePools.take<CF, DG>();
  result->addRhs(degr);
  for (size_t i = 0; i < size(); ++i) {
    result->addLhs(aux::abs(data[i].c), data[i].l);
  }
  result->orig = getOrigin();
  if (result->plogger) result->resetBuffer(id);
  return result;
}

template <typename CF, typename DG>
CeSuper Watched<CF, DG>::toExpanded(ConstrExpPools& cePools) const {
  return expandTo(cePools);
}

template <typename CF, typename DG>
void CountingSafe<CF, DG>::initializeWatches(CRef cr, Solver& solver) {
  auto& Level = solver.Level;
  auto& Pos = solver.Pos;
  auto& adj = solver.adj;
  auto& qhead = solver.qhead;

  DG& slk = *slack;
  slk = -*degr;
  unsigned int length = size();
  for (unsigned int i = 0; i < length; ++i) {
    Lit l = terms[i].l;
    adj[l].emplace_back(cr, i + INF);
    if (!isFalse(Level, l) || Pos[toVar(l)] >= qhead) slk += terms[i].c;
  }

  assert(slk >= 0);
  assert(hasCorrectSlack(solver));
  if (slk < terms[0].c) {  // propagate
    for (unsigned int i = 0; i < length && terms[i].c > slk; ++i) {
      Lit l = terms[i].l;
      if (isUnknown(Pos, l)) {
        assert(isCorrectlyPropagating(solver, i));
        solver.propagate(l, cr);
      }
    }
  }
}

template <typename CF, typename DG>
WatchStatus CountingSafe<CF, DG>::checkForPropagation(CRef cr, int& idx, [[maybe_unused]] Lit p, Solver& solver) {
  auto& Pos = solver.Pos;

  assert(idx >= INF);
  assert(terms[idx - INF].l == p);
  const unsigned int length = size();
  const CF& lrgstCf = terms[0].c;
  const CF& c = terms[idx - INF].c;

  DG& slk = *slack;
  slk -= c;
  assert(hasCorrectSlack(solver));

  if (slk < 0) {
    assert(isCorrectlyConflicting(solver));
    return WatchStatus::CONFLICTING;
  }
  if (slk < lrgstCf) {
    if (!options.propIdx || ntrailpops < stats.NTRAILPOPS) {
      ntrailpops = stats.NTRAILPOPS;
      watchIdx = 0;
    }
    stats.NPROPCHECKS -= watchIdx;
    for (; watchIdx < length && terms[watchIdx].c > slk; ++watchIdx) {
      const Lit l = terms[watchIdx].l;
      if (isUnknown(Pos, l)) {
        stats.NPROPCLAUSE += (*degr == 1);
        stats.NPROPCARD += (*degr != 1 && lrgstCf == 1);
        ++stats.NPROPCOUNTING;
        assert(isCorrectlyPropagating(solver, watchIdx));
        solver.propagate(l, cr);
      }
    }
    stats.NPROPCHECKS += watchIdx;
  }
  return WatchStatus::KEEPWATCH;
}

template <typename CF, typename DG>
void CountingSafe<CF, DG>::undoFalsified(int i) {
  assert(i >= INF);
  DG& slk = *slack;
  slk += terms[i - INF].c;
  ++stats.NWATCHLOOKUPSBJ;
}

template <typename CF, typename DG>
void CountingSafe<CF, DG>::resolveWith(CeSuper confl, Lit l, IntSet* actSet, Solver& solver) {
  confl->resolveWith(expandTo(solver.cePools), l, actSet, solver.getLevel(), solver.getPos());
}

template <typename CF, typename DG>
CePtr<ConstrExp<CF, DG>> CountingSafe<CF, DG>::expandTo(ConstrExpPools& cePools) const {
  CePtr<ConstrExp<CF, DG>> result = cePools.take<CF, DG>();
  result->addRhs(*degr);
  for (size_t i = 0; i < size(); ++i) {
    result->addLhs(terms[i].c, terms[i].l);
  }
  result->orig = getOrigin();
  if (result->plogger) result->resetBuffer(id);
  return result;
}

template <typename CF, typename DG>
CeSuper CountingSafe<CF, DG>::toExpanded(ConstrExpPools& cePools) const {
  return expandTo(cePools);
}

template <typename CF, typename DG>
void WatchedSafe<CF, DG>::initializeWatches(CRef cr, Solver& solver) {
  auto& Level = solver.Level;
  auto& Pos = solver.Pos;
  auto& adj = solver.adj;
  auto& qhead = solver.qhead;

  DG& wslk = *watchslack;
  wslk = -*degr;
  unsigned int length = size();
  const CF lrgstCf = aux::abs(terms[0].c);
  for (unsigned int i = 0; i < length && wslk < lrgstCf; ++i) {
    Lit l = terms[i].l;
    if (!isFalse(Level, l) || Pos[toVar(l)] >= qhead) {
      assert(terms[i].c > 0);
      wslk += terms[i].c;
      terms[i].c = -terms[i].c;
      adj[l].emplace_back(cr, i + INF);
    }
  }
  assert(wslk >= 0);
  assert(hasCorrectSlack(solver));
  if (wslk < lrgstCf) {
    // set sufficient falsified watches
    std::vector<unsigned int> falsifiedIdcs;
    falsifiedIdcs.reserve(length);
    for (unsigned int i = 0; i < length; ++i)
      if (isFalse(Level, terms[i].l) && Pos[toVar(terms[i].l)] < qhead) falsifiedIdcs.push_back(i);
    std::sort(falsifiedIdcs.begin(), falsifiedIdcs.end(),
              [&](unsigned int i1, unsigned int i2) { return Pos[toVar(terms[i1].l)] > Pos[toVar(terms[i2].l)]; });
    DG diff = lrgstCf - wslk;
    for (unsigned int i : falsifiedIdcs) {
      assert(terms[i].c > 0);
      diff -= terms[i].c;
      terms[i].c = -terms[i].c;
      adj[terms[i].l].emplace_back(cr, i + INF);
      if (diff <= 0) break;
    }
    // perform initial propagation
    for (unsigned int i = 0; i < length && aux::abs(terms[i].c) > wslk; ++i)
      if (isUnknown(Pos, terms[i].l)) {
        assert(isCorrectlyPropagating(solver, i));
        solver.propagate(terms[i].l, cr);
      }
  }
}

template <typename CF, typename DG>
WatchStatus WatchedSafe<CF, DG>::checkForPropagation(CRef cr, int& idx, Lit p, Solver& solver) {
  _unused(p);
  auto& Level = solver.Level;
  auto& Pos = solver.Pos;
  auto& adj = solver.adj;

  assert(idx >= INF);
  assert(terms[idx - INF].l == p);
  const unsigned int length = size();
  const CF lrgstCf = aux::abs(terms[0].c);
  CF& c = terms[idx - INF].c;

  if (!options.propIdx || ntrailpops < stats.NTRAILPOPS) {
    ntrailpops = stats.NTRAILPOPS;
    watchIdx = 0;
  }

  assert(c < 0);
  DG& wslk = *watchslack;
  wslk += c;
  if (!options.propSup || wslk - c >= lrgstCf) {  // look for new watches if previously, slack was at least lrgstCf
    stats.NWATCHCHECKS -= watchIdx;
    for (; watchIdx < length && wslk < lrgstCf; ++watchIdx) {
      const CF& cf = terms[watchIdx].c;
      const Lit l = terms[watchIdx].l;
      if (cf > 0 && !isFalse(Level, l)) {
        wslk += cf;
        terms[watchIdx].c = -cf;
        adj[l].emplace_back(cr, watchIdx + INF);
      }
    }  // NOTE: first innermost loop of RoundingSat
    stats.NWATCHCHECKS += watchIdx;
    if (wslk < lrgstCf) {
      assert(watchIdx == length);
      watchIdx = 0;
    }
  }  // else we did not find enough watches last time, so we can skip looking for them now

  assert(hasCorrectSlack(solver));
  assert(hasCorrectWatches(solver));

  if (wslk >= lrgstCf) {
    c = -c;
    return WatchStatus::DROPWATCH;
  }
  if (wslk < 0) {
    assert(isCorrectlyConflicting(solver));
    return WatchStatus::CONFLICTING;
  }
  // keep the watch, check for propagation
  stats.NPROPCHECKS -= watchIdx;
  for (; watchIdx < length && aux::abs(terms[watchIdx].c) > wslk; ++watchIdx) {
    const Lit l = terms[watchIdx].l;
    if (isUnknown(Pos, l)) {
      stats.NPROPCLAUSE += (*degr == 1);
      stats.NPROPCARD += (*degr != 1 && lrgstCf == 1);
      ++stats.NPROPWATCH;
      assert(isCorrectlyPropagating(solver, watchIdx));
      solver.propagate(l, cr);
    }  // NOTE: second innermost loop of RoundingSat
  }
  stats.NPROPCHECKS += watchIdx;
  return WatchStatus::KEEPWATCH;
}

template <typename CF, typename DG>
void WatchedSafe<CF, DG>::undoFalsified(int i) {
  assert(i >= INF);
  assert(terms[i - INF].c < 0);
  DG& wslk = *watchslack;
  wslk -= terms[i - INF].c;
  ++stats.NWATCHLOOKUPSBJ;
}

template <typename CF, typename DG>
void WatchedSafe<CF, DG>::resolveWith(CeSuper confl, Lit l, IntSet* actSet, Solver& solver) {
  confl->resolveWith(expandTo(solver.cePools), l, actSet, solver.getLevel(), solver.getPos());
}

template <typename CF, typename DG>
CePtr<ConstrExp<CF, DG>> WatchedSafe<CF, DG>::expandTo(ConstrExpPools& cePools) const {
  CePtr<ConstrExp<CF, DG>> result = cePools.take<CF, DG>();
  result->addRhs(*degr);
  for (size_t i = 0; i < size(); ++i) {
    result->addLhs(aux::abs(terms[i].c), terms[i].l);
  }
  result->orig = getOrigin();
  if (result->plogger) result->resetBuffer(id);
  return result;
}

template <typename CF, typename DG>
CeSuper WatchedSafe<CF, DG>::toExpanded(ConstrExpPools& cePools) const {
  return expandTo(cePools);
}

// TODO: keep below test methods?

bool Constr::isCorrectlyConflicting(const Solver& solver) {
  return true;  // comment to run check
  BigVal slack = -degree();
  for (int i = 0; i < (int)size(); ++i) {
    slack += isFalse(solver.getLevel(), lit(i)) ? 0 : coef(i);
  }
  return slack < 0;
}

bool Constr::isCorrectlyPropagating(const Solver& solver, int idx) {
  return true;  // comment to run check
  assert(isUnknown(solver.getPos(), lit(idx)));
  BigVal slack = -degree();
  for (int i = 0; i < (int)size(); ++i) {
    slack += isFalse(solver.getLevel(), lit(i)) ? 0 : coef(i);
  }
  return slack < coef(idx);
}

void Constr::print(const Solver& solver) {
  for (size_t i = 0; i < size(); ++i) {
    int pos = solver.getPos()[toVar(lit(i))];
    std::cout << coef(i) << "x" << lit(i)
              << (pos < solver.qhead ? (isTrue(solver.getLevel(), lit(i)) ? "t" : "f") : "u") << (pos == INF ? -1 : pos)
              << " ";
  }
  std::cout << ">= " << degree() << std::endl;
}

template <typename CF, typename DG>
bool Counting<CF, DG>::hasCorrectSlack(const Solver& solver) {
  return true;  // comment to run check
  DG slk = -degr;
  for (int i = 0; i < (int)size(); ++i) {
    if (solver.getPos()[toVar(lit(i))] >= solver.qhead || !isFalse(solver.getLevel(), lit(i))) {
      slk += data[i].c;
    }
  }
  return (slk == slack);
}

template <typename CF, typename DG>
bool Watched<CF, DG>::hasCorrectSlack(const Solver& solver) {
  return true;  // comment to run check
  DG slk = -degr;
  for (int i = 0; i < (int)size(); ++i) {
    if (data[i].c < 0 && (solver.getPos()[toVar(lit(i))] >= solver.qhead || !isFalse(solver.getLevel(), lit(i))))
      slk += aux::abs(data[i].c);
  }
  return (slk == watchslack);
}

template <typename CF, typename DG>
bool CountingSafe<CF, DG>::hasCorrectSlack(const Solver& solver) {
  return true;  // comment to run check
  DG slk = -*degr;
  for (int i = 0; i < (int)size(); ++i) {
    if (solver.getPos()[toVar(lit(i))] >= solver.qhead || !isFalse(solver.getLevel(), lit(i))) {
      slk += terms[i].c;
    }
  }
  return (slk == *slack);
}

template <typename CF, typename DG>
bool WatchedSafe<CF, DG>::hasCorrectSlack(const Solver& solver) {
  return true;  // comment to run check
  DG slk = -*degr;
  for (int i = 0; i < (int)size(); ++i) {
    if (terms[i].c < 0 && (solver.getPos()[toVar(lit(i))] >= solver.qhead || !isFalse(solver.getLevel(), lit(i))))
      slk += aux::abs(terms[i].c);
  }
  return (slk == *watchslack);
}

template <typename CF, typename DG>
bool Watched<CF, DG>::hasCorrectWatches(const Solver& solver) {
  return true;  // comment to run check
  if (watchslack >= aux::abs(data[0].c)) return true;
  for (int i = 0; i < (int)watchIdx; ++i) assert(!isUnknown(solver.getPos(), lit(i)));
  for (int i = 0; i < (int)size(); ++i) {
    if (!(data[i].c < 0 || isFalse(solver.getLevel(), data[i].l))) {
      std::cout << i << " " << data[i].c << " " << isFalse(solver.getLevel(), data[i].l) << std::endl;
      print(solver);
    }
    assert(data[i].c < 0 || isFalse(solver.getLevel(), data[i].l));
  }
  return true;
}

template <typename CF, typename DG>
bool WatchedSafe<CF, DG>::hasCorrectWatches(const Solver& solver) {
  return true;  // comment to run check
  if (*watchslack >= aux::abs(terms[0].c)) return true;
  for (int i = 0; i < (int)watchIdx; ++i) assert(!isUnknown(solver.getPos(), lit(i)));
  for (int i = 0; i < (int)size(); ++i) {
    if (!(terms[i].c < 0 || isFalse(solver.getLevel(), terms[i].l))) {
      std::cout << i << " " << terms[i].c << " " << isFalse(solver.getLevel(), terms[i].l) << std::endl;
      print(solver);
    }
    assert(terms[i].c < 0 || isFalse(solver.getLevel(), terms[i].l));
  }
  return true;
}

template struct Counting<int, long long>;

template struct Watched<int, long long>;

template struct CountingSafe<long long, int128>;
template struct CountingSafe<int128, int128>;
template struct CountingSafe<bigint, bigint>;

template struct WatchedSafe<long long, int128>;
template struct WatchedSafe<int128, int128>;
template struct WatchedSafe<bigint, bigint>;

}  // namespace rs
