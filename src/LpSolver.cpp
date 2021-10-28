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

#include "LpSolver.hpp"
#include <queue>
#include "Solver.hpp"

namespace rs {

#if WITHSOPLEX

CandidateCut::CandidateCut(CeSuper in, const std::vector<double>& sol) {
  assert(in->isSaturated());
  in->saturateAndFixOverflowRational(sol);
  in->toSimple()->copyTo(simpcons);
  // NOTE: simpcons is already in var-normal form
  initialize(sol);
}

CandidateCut::CandidateCut(const Constr& in, CRef cref, const std::vector<double>& sol, ConstrExpPools& pools)
    : cr(cref) {
  assert(in.degree() > 0);
  CeSuper tmp = in.toExpanded(pools);
  tmp->saturateAndFixOverflowRational(sol);
  if (tmp->isTautology()) {
    return;
  }
  tmp->toSimple()->copyTo(simpcons);
  // NOTE: simpcons is already in var-normal form
  initialize(sol);
  assert(cr != CRef_Undef);
}

void CandidateCut::initialize(const std::vector<double>& sol) {
  std::sort(simpcons.terms.begin(), simpcons.terms.end(),
            [](const Term<long long>& t1, const Term<long long>& t2) { return t1.l < t2.l; });
  assert(norm == 1);
  norm = 0;
  for (const Term<long long>& p : simpcons.terms) norm += (double)p.c * (double)p.c;
  norm = std::sqrt(norm);
  ratSlack = -static_cast<double>(simpcons.rhs);
  for (Term<long long>& p : simpcons.terms) {
    assert(p.l > 0);  // simpcons is in var-normal form
    ratSlack += (double)p.c * sol[p.l];
  }
  assert(norm >= 0);
  if (norm == 0) norm = 1;
  ratSlack /= norm;
}

// @pre: simpcons is ordered and norm is calculated
double CandidateCut::cosOfAngleTo(const CandidateCut& other) const {
  assert(norm != 0);
  assert(other.norm != 0);
  double cos = 0;
  int i = 0;
  int j = 0;
  while (i < (int)simpcons.terms.size() && j < (int)other.simpcons.terms.size()) {
    int x = simpcons.terms[i].l;
    int y = other.simpcons.terms[j].l;
    if (x < y)
      ++i;
    else if (x > y)
      ++j;
    else {  // x==y
      cos += (double)simpcons.terms[i].c * (double)other.simpcons.terms[j].c;
      ++i;
      ++j;
    }
  }
  return cos / (norm * other.norm);
}

std::ostream& operator<<(std::ostream& o, const CandidateCut& cc) {
  return o << cc.simpcons << " norm " << cc.norm << " ratSlack " << cc.ratSlack;
}

LpSolver::LpSolver(Solver& slvr, const CeArb obj) : solver(slvr) {
  assert(INFTY == lp.realParam(lp.INFTY));

  if (options.verbosity.get() > 1) std::cout << "c Initializing LP" << std::endl;
  setNbVariables(slvr.getNbVars() + 1);
  lp.setIntParam(soplex::SoPlex::SYNCMODE, soplex::SoPlex::SYNCMODE_ONLYREAL);
  lp.setIntParam(soplex::SoPlex::SOLVEMODE, soplex::SoPlex::SOLVEMODE_REAL);
  lp.setIntParam(soplex::SoPlex::CHECKMODE, soplex::SoPlex::CHECKMODE_REAL);
  lp.setIntParam(soplex::SoPlex::SIMPLIFIER, soplex::SoPlex::SIMPLIFIER_OFF);
  lp.setIntParam(soplex::SoPlex::OBJSENSE, soplex::SoPlex::OBJSENSE_MINIMIZE);
  lp.setIntParam(soplex::SoPlex::VERBOSITY, options.verbosity.get());
  lp.setRandomSeed(0);

  // add two empty rows for objective bound constraints
  while (row2data.size() < 2) {
    soplex::DSVectorReal row(0);
    lp.addRowReal(soplex::LPRowReal(row, soplex::LPRowReal::Type::GREATER_EQUAL, 0));
    row2data.emplace_back(ID_Trivial, false);
  }

  // add all formula constraints
  for (CRef cr : solver.constraints)
    if (solver.ca[cr].getOrigin() == Origin::FORMULA) addConstraint(cr, false);

  // NOTE: scaling objective is not needed, as if it does not fit in double (i.e. >1e300), it will still be sound.
  soplex::DVectorReal objective;
  objective.reDim(getNbVariables());  // NOTE: automatically set to zero
  if (obj->vars.size() > 0)
    for (Var v : obj->vars) objective[v] = static_cast<double>(obj->coefs[v]);
  else
    for (int v = 1; v < getNbVariables(); ++v) objective[v] = 1;  // add default objective function
  lp.changeObjReal(objective);

  if (options.verbosity.get() > 1) std::cout << "c Finished initializing LP" << std::endl;
}

void LpSolver::setNbVariables(int n) {
  if (n <= getNbVariables()) return;
  soplex::LPColSetReal allCols;
  allCols.reMax(n - getNbVariables());
  soplex::DSVectorReal dummycol(0);
  for (Var v = getNbVariables(); v < n; ++v) {
    allCols.add(soplex::LPColReal(0, dummycol, 1, 0));
  }
  lp.addColsReal(allCols);

  lpSol.reDim(n);
  lpSolution.resize(n, 0);
  lowerBounds.reDim(n);
  upperBounds.reDim(n);
  assert(getNbVariables() == n);
}

int LpSolver::getNbVariables() const { return lp.numCols(); }
int LpSolver::getNbRows() const { return lp.numRows(); }

CeSuper LpSolver::createLinearCombinationFarkas(soplex::DVectorReal& mults) {
  double scale = getScaleFactor(mults, true);
  if (scale == 0) return CeNull();
  assert(scale > 0);

  CeArb out = solver.cePools.takeArb();
  for (int r = 0; r < mults.dim(); ++r) {
    bigint factor = static_cast<bigint>(mults[r] * scale);
    if (factor <= 0) continue;
    assert(lp.lhsReal(r) != INFTY);
    out->addUp(rowToConstraint(r), factor);
  }
  out->removeUnitsAndZeroes(solver.getLevel(), solver.getPos());
  assert(out->hasNoZeroes());
  out->weakenSmalls(out->absCoeffSum() / static_cast<bigint>((double)out->vars.size() / options.lpIntolerance.get()));
  out->saturateAndFixOverflow(solver.getLevel(), (bool)options.weakenFull, options.bitsOverflow.get(),
                              options.bitsReduced.get(), 0);
  return out;
}

CandidateCut LpSolver::createLinearCombinationGomory(soplex::DVectorReal& mults) {
  double scale = getScaleFactor(mults, false);
  if (scale == 0) return CandidateCut();
  assert(scale > 0);
  CeArb lcc = solver.cePools.takeArb();

  std::vector<std::pair<BigCoef, int>> slacks;
  for (int r = 0; r < mults.dim(); ++r) {
    bigint factor = static_cast<bigint>(mults[r] * scale);
    if (factor == 0) continue;
    Ce64 ce = rowToConstraint(r);
    if (factor < 0) ce->invert();
    lcc->addUp(ce, aux::abs(factor));
    slacks.emplace_back(-factor, r);
  }

  BigVal b = lcc->getRhs();
  for (Var v : lcc->vars)
    if (lpSolution[v] > 0.5) b -= lcc->coefs[v];
  if (b == 0) {
    return CandidateCut();
  }

  assert(scale > 0);
  bigint divisor = static_cast<bigint>(std::ceil(scale));
  while ((b % divisor) == 0) ++divisor;
  lcc->applyMIR(divisor, [&](Var v) -> Lit { return lpSolution[v] <= 0.5 ? v : -v; });

  // round up the slack variables MIR style and cancel out the slack variables
  bigint bmodd = aux::mod_safe(b, divisor);
  for (unsigned i = 0; i < slacks.size(); ++i) {
    bigint factor =
        bmodd * aux::floordiv_safe(slacks[i].first, divisor) + std::min(aux::mod_safe(slacks[i].first, divisor), bmodd);
    if (factor == 0) continue;
    Ce64 ce = rowToConstraint(slacks[i].second);
    if (factor < 0) ce->invert();
    lcc->addUp(ce, aux::abs(factor));
  }
  if (lcc->plogger) lcc->logAsInput();
  // TODO: fix logging for Gomory cuts

  lcc->removeUnitsAndZeroes(solver.getLevel(), solver.getPos());
  if (lcc->isTautology())
    lcc->reset();
  else {
    assert(lcc->hasNoZeroes());
    lcc->weakenSmalls(lcc->absCoeffSum() / static_cast<bigint>((double)lcc->vars.size() / options.lpIntolerance.get()));
  }
  CandidateCut result(lcc, lpSolution);
  return result;
}

void LpSolver::constructGomoryCandidates() {
  std::vector<int> indices;
  indices.resize(getNbRows());
  lp.getBasisInd(indices.data());

  assert(lpSlackSolution.dim() == getNbRows());
  std::vector<std::pair<double, int>> fracrowvec;
  for (int row = 0; row < getNbRows(); ++row) {
    if (asynch_interrupt) throw asynchInterrupt;
    double fractionality = 0;
    if (indices[row] >= 0) {  // basic original variable / column
      assert(indices[row] < (int)lpSolution.size());
      fractionality = nonIntegrality(lpSolution[indices[row]]);
    } else {  // basic slack variable / row
      assert(-indices[row] - 1 < lpSlackSolution.dim());
      fractionality = nonIntegrality(lpSlackSolution[-indices[row] - 1]);
    }
    assert(fractionality >= 0);
    if (fractionality > 0) fracrowvec.emplace_back(fractionality, row);
  }
  std::priority_queue<std::pair<double, int>> fracrows(std::less<std::pair<double, int>>(), fracrowvec);

  [[maybe_unused]] double last = 0.5;
  for (int i = 0; i < options.gomoryCutLimit.get() && !fracrows.empty(); ++i) {
    assert(last >= fracrows.top().first);
    last = fracrows.top().first;
    int row = fracrows.top().second;
    fracrows.pop();

    assert(lpMultipliers.dim() == getNbRows());
    lpMultipliers.clear();
    lp.getBasisInverseRowReal(row, lpMultipliers.get_ptr());
    candidateCuts.push_back(createLinearCombinationGomory(lpMultipliers));
    if (candidateCuts.back().ratSlack >= -options.lpIntolerance.get()) candidateCuts.pop_back();
    for (int i = 0; i < lpMultipliers.dim(); ++i) lpMultipliers[i] = -lpMultipliers[i];
    candidateCuts.push_back(createLinearCombinationGomory(lpMultipliers));
    if (candidateCuts.back().ratSlack >= -options.lpIntolerance.get()) candidateCuts.pop_back();
  }
}

void LpSolver::constructLearnedCandidates() {
  for (CRef cr : solver.constraints) {
    if (asynch_interrupt) throw asynchInterrupt;
    const Constr& c = solver.ca[cr];
    if (c.getOrigin() == Origin::LEARNED || c.getOrigin() == Origin::LEARNEDFARKAS || c.getOrigin() == Origin::GOMORY) {
      bool containsNewVars = false;
      for (unsigned int i = 0; i < c.size() && !containsNewVars; ++i) {
        containsNewVars = toVar(c.lit(i)) >= getNbVariables();
        assert((toVar(c.lit(i)) > solver.getNbOrigVars()) == containsNewVars);
        // for now, getNbVariables() == solver.getNbOrigVars().nbOrigVars+1
      }
      if (containsNewVars) continue;
      candidateCuts.emplace_back(c, cr, lpSolution, solver.cePools);
      if (candidateCuts.back().ratSlack >= -options.lpIntolerance.get()) candidateCuts.pop_back();
    }
  }
}

void LpSolver::addFilteredCuts() {
  for ([[maybe_unused]] const CandidateCut& cc : candidateCuts) {
    assert(cc.norm != 0);
  }
  std::sort(candidateCuts.begin(), candidateCuts.end(), [](const CandidateCut& x1, const CandidateCut& x2) {
    return x1.ratSlack > x2.ratSlack ||
           (x1.ratSlack == x2.ratSlack && x1.simpcons.terms.size() < x2.simpcons.terms.size());
  });

  // filter the candidate cuts
  std::vector<int> keptCuts;  // indices
  for (unsigned int i = 0; i < candidateCuts.size(); ++i) {
    bool parallel = false;
    for (unsigned int j = 0; j < keptCuts.size() && !parallel; ++j) {
      if (asynch_interrupt) throw asynchInterrupt;
      parallel = candidateCuts[keptCuts[j]].cosOfAngleTo(candidateCuts[i]) > options.maxCutCos.get();
    }
    if (!parallel) keptCuts.push_back(i);
  }

  for (int i : keptCuts) {
    CandidateCut& cc = candidateCuts[i];
    CeSuper ce = cc.simpcons.toExpanded(solver.cePools);
    ce->postProcess(solver.getLevel(), solver.getPos(), true, stats);
    assert(ce->fitsInDouble());
    assert(!ce->isTautology());
    if (cc.cr == CRef_Undef) {  // Gomory cut
      solver.learnConstraint(ce, Origin::GOMORY);
    } else {  // learned cut
      ++stats.NLPLEARNEDCUTS;
    }
    addConstraint(ce, true);
  }
}

void LpSolver::pruneCuts() {
  assert(getNbRows() == (int)row2data.size());
  lpMultipliers.clear();
  if (!lp.getDual(lpMultipliers)) return;
  for (int r = 0; r < getNbRows(); ++r)
    if (row2data[r].removable && lpMultipliers[r] == 0) {
      ++stats.NLPDELETEDCUTS;
      toRemove.push_back(r);
    }
}

// NOTE: it is possible that mults are negative (e.g., when calculating Gomory cuts)
double LpSolver::getScaleFactor(soplex::DVectorReal& mults, bool removeNegatives) {
  double largest = 0;
  for (int i = 0; i < mults.dim(); ++i) {
    if (std::isnan(mults[i]) || std::isinf(mults[i]) || (removeNegatives && mults[i] < 0)) mults[i] = 0;
    largest = std::max(aux::abs(mults[i]), largest);
  }
  if (largest == 0) return 0;
  return maxMult / largest;
}

Ce64 LpSolver::rowToConstraint(int row) {
  Ce64 ce = solver.cePools.take64();
  double rhs = lp.lhsReal(row);
  assert(aux::abs(rhs) != INFTY);
  assert(validVal(rhs));
  ce->addRhs((long long)rhs);

  lpRow.clear();
  lp.getRowVectorReal(row, lpRow);
  for (int i = 0; i < lpRow.size(); ++i) {
    const soplex::Nonzero<double>& el = lpRow.element(i);
    assert(validVal(el.val));
    assert(el.val != 0);
    ce->addLhs((long long)el.val, el.idx);
  }
  if (ce->plogger) ce->resetBuffer(row2data[row].id);
  return ce;
}

std::pair<LpStatus, CeSuper> LpSolver::checkFeasibility(bool inProcessing) {
  if (solver.logger) solver.logger->logComment("Checking LP", stats);
  if (options.lpPivotRatio.get() < 0)
    lp.setIntParam(soplex::SoPlex::ITERLIMIT, -1);  // no pivot limit
  else if (options.lpPivotRatio.get() * stats.NCONFL < (inProcessing ? stats.NLPPIVOTSROOT : stats.NLPPIVOTSINTERNAL))
    return {LpStatus::PIVOTLIMIT, CeNull()};  // pivot ratio exceeded
  else
    lp.setIntParam(soplex::SoPlex::ITERLIMIT, options.lpPivotBudget.get() * lpPivotMult);
  flushConstraints();

  // Set the  LP's bounds based on the current trail
  for (Var v = 1; v < getNbVariables(); ++v) {
    lowerBounds[v] = isTrue(solver.getLevel(), v);
    upperBounds[v] = !isFalse(solver.getLevel(), v);
  }
  lp.changeBoundsReal(lowerBounds, upperBounds);

  // Run the LP
  soplex::SPxSolver::Status stat;
  stat = lp.optimize();
  ++stats.NLPCALLS;
  if (inProcessing)
    stats.NLPPIVOTSROOT += lp.numIterations();
  else
    stats.NLPPIVOTSINTERNAL += lp.numIterations();
  stats.LPSOLVETIME += lp.solveTime();

  if (options.verbosity.get() > 1)
    std::cout << "c " << (inProcessing ? "root" : "internal") << " LP status: " << stat << std::endl;
  assert(stat != soplex::SPxSolver::Status::NO_PROBLEM);
  assert(stat <= soplex::SPxSolver::Status::OPTIMAL_UNSCALED_VIOLATIONS);
  assert(stat >= soplex::SPxSolver::Status::ERROR);

  if (stat == soplex::SPxSolver::Status::ABORT_ITER) {
    lpPivotMult *= 2;  // increase pivot budget when calling the LP solver
    return {LpStatus::PIVOTLIMIT, CeNull()};
  }

  if (stat == soplex::SPxSolver::Status::OPTIMAL) {
    ++stats.NLPOPTIMAL;
    if (!lp.hasSol()) {
      ++stats.NLPNOPRIMAL;
      resetBasis();
    }
    if (lp.numIterations() == 0) ++stats.NLPNOPIVOT;
    return {LpStatus::OPTIMAL, CeNull()};
  }

  if (stat == soplex::SPxSolver::Status::ABORT_CYCLING) {
    ++stats.NLPCYCLING;
    resetBasis();
    return {LpStatus::UNDETERMINED, CeNull()};
  }
  if (stat == soplex::SPxSolver::Status::SINGULAR) {
    ++stats.NLPSINGULAR;
    resetBasis();
    return {LpStatus::UNDETERMINED, CeNull()};
  }
  if (stat != soplex::SPxSolver::Status::INFEASIBLE) {
    ++stats.NLPOTHER;
    resetBasis();
    return {LpStatus::UNDETERMINED, CeNull()};
  }

  // Infeasible LP :)
  ++stats.NLPINFEAS;

  // To prove that we have an inconsistency, let's build the Farkas proof
  if (!lp.getDualFarkas(lpMultipliers)) {
    ++stats.NLPNOFARKAS;
    resetBasis();
    return {LpStatus::UNDETERMINED, CeNull()};
  }

  CeSuper confl = createLinearCombinationFarkas(lpMultipliers);
  if (!confl) return {LpStatus::UNDETERMINED, CeNull()};
  solver.learnConstraint(confl, Origin::FARKAS);
  if (confl->hasNegativeSlack(solver.getLevel())) return {LpStatus::INFEASIBLE, confl};
  return {LpStatus::UNDETERMINED, CeNull()};
}

void LpSolver::inProcess() {
  assert(solver.decisionLevel() == 0);
  std::pair<LpStatus, CeSuper> lpResult = checkFeasibility(true);
  LpStatus lpstat = lpResult.first;
  [[maybe_unused]] CeSuper confl = lpResult.second;
  assert((lpstat == LpStatus::INFEASIBLE) == (confl && confl->hasNegativeSlack(solver.getLevel())));
  // NOTE: we don't handle confl here, as it is added as a learned constraint already.
  if (lpstat != LpStatus::OPTIMAL) return;
  if (!lp.hasSol()) return;
  lp.getPrimal(lpSol);
  assert(lpSol.dim() == (int)lpSolution.size());
  for (int i = 0; i < lpSol.dim(); ++i) lpSolution[i] = lpSol[i];
  lp.getSlacksReal(lpSlackSolution);
  assert((int)solver.phase.size() >= getNbVariables());
  for (Var v = 1; v < getNbVariables(); ++v) solver.phase[v] = (lpSolution[v] <= 0.5) ? -v : v;
  if (options.verbosity.get() > 0) std::cout << "c rational objective " << lp.objValueReal() << std::endl;
  candidateCuts.clear();
  if (solver.logger && (options.addGomoryCuts || options.addLearnedCuts)) solver.logger->logComment("cutting", stats);
  if (options.addLearnedCuts) constructLearnedCandidates();  // first to avoid adding gomory cuts twice
  if (options.addGomoryCuts) constructGomoryCandidates();
  addFilteredCuts();
  pruneCuts();
}

void LpSolver::resetBasis() {
  ++stats.NLPRESETBASIS;
  lp.clearBasis();  // and hope next iteration works fine
}

void LpSolver::convertConstraint(const ConstrSimple64& c, soplex::DSVectorReal& row, double& rhs) {
  assert(row.max() >= (int)c.terms.size());
  for (auto& t : c.terms) {
    if (t.c == 0) continue;
    assert(t.l > 0);
    assert(t.l < lp.numColsReal());
    assert(t.c < INFLPINT);
    row.add(t.l, t.c);
  }
  rhs = static_cast<double>(c.rhs);
  assert(validVal(rhs));
}

void LpSolver::addConstraint(CeSuper c, bool removable, bool upperbound, bool lowerbound) {
  assert(!upperbound || c->orig == Origin::UPPERBOUND);
  c->saturateAndFixOverflowRational(lpSolution);
  ID id =
      solver.logger ? c->logProofLineWithInfo("LP", stats) : ++solver.crefID;  // TODO: fix this kind of logger check
  if (upperbound || lowerbound) {
    boundsToAdd[lowerbound].id = id;
    c->toSimple()->copyTo(boundsToAdd[lowerbound].cs);
  } else {
    toAdd[id] = {ConstrSimple64(), removable};
    c->toSimple()->copyTo(toAdd[id].cs);
  }
}

void LpSolver::addConstraint(CRef cr, bool removable, bool upperbound, bool lowerbound) {
  assert(cr != CRef_Undef);
  addConstraint(solver.ca[cr].toExpanded(solver.cePools), removable, upperbound, lowerbound);
}

void LpSolver::flushConstraints() {
  if (toRemove.size() > 0) {  // first remove rows
    std::vector<int> rowsToRemove(getNbRows(), 0);
    for (int row : toRemove) {
      stats.NLPDELETEDROWS += (rowsToRemove[row] == 0);
      assert(row < (int)rowsToRemove.size());
      rowsToRemove[row] = -1;
    }
    lp.removeRowsReal(rowsToRemove.data());  // TODO: use other removeRowsReal method?
    for (int r = 0; r < (int)rowsToRemove.size(); ++r) {
      int newrow = rowsToRemove[r];
      if (newrow < 0 || newrow == r) continue;
      row2data[newrow] = row2data[r];
    }
    row2data.resize(getNbRows());
    toRemove.clear();
  }

  if (toAdd.size() > 0) {  // then add rows
    soplex::LPRowSetReal rowsToAdd;
    rowsToAdd.reMax(toAdd.size());
    row2data.reserve(row2data.size() + toAdd.size());
    for (auto& p : toAdd) {
      double rhs;
      soplex::DSVectorReal row(p.second.cs.terms.size());
      convertConstraint(p.second.cs, row, rhs);
      rowsToAdd.add(soplex::LPRowReal(row, soplex::LPRowReal::Type::GREATER_EQUAL, rhs));
      row2data.emplace_back(p.first, p.second.removable);
      ++stats.NLPADDEDROWS;
    }
    lp.addRowsReal(rowsToAdd);
    toAdd.clear();
  }

  for (int i = 0; i < 2; ++i) {
    if (boundsToAdd[i].id == row2data[i].id) continue;
    double rhs;
    soplex::DSVectorReal row(boundsToAdd[i].cs.terms.size());
    convertConstraint(boundsToAdd[i].cs, row, rhs);
    lp.changeRowReal(i, soplex::LPRowReal(row, soplex::LPRowReal::Type::GREATER_EQUAL, rhs));
    row2data[i] = {boundsToAdd[i].id, false};  // so upper bound resides in row[0]
  }

  lpSlackSolution.reDim(getNbRows());
  lpMultipliers.reDim(getNbRows());
  assert((int)row2data.size() == getNbRows());
}

#endif  // WITHSOPLEX

}  // namespace rs
