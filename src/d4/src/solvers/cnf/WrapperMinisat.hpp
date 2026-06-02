/*
 * d4
 * Copyright (C) 2020  Univ. Artois & CNRS
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this library; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA
 */
#pragma once

#include "../WrapperSolver.hpp"
#include "minisat/Solver.hpp"
#include "src/problem/ProblemManager.hpp"
#include "src/problem/ProblemTypes.hpp"
#include "../../../../d4_smt.hpp"

namespace d4 {
class WrapperMinisat : public WrapperSolver {
 protected:
  minisat::Solver m_solver;
  minisat::vec<minisat::Var> m_setOfVar_m;
  int m_smtTrail = 0;
  // Variables of the currently solved connected component. Used to
  // determine which projection variables must be assigned before invoking
  // SMT consistency checks and to count decompositions.
  std::vector<Var> m_currentComponent;
  // Snapshot of the Minisat trail that has already been translated into
  // theory assertions.
  std::vector<minisat::Lit> m_syncedTrailLits;
  // Number of assertions currently pushed on the SMT solver stack.
  int m_smtContextLevel = 0;
  // Indicates whether newly asserted theory literals require an SMT check.
  bool m_pendingTheoryCheck = false;

  using WrapperSolver::m_isInAssumption;

 public:
  ~WrapperMinisat() override {}

  void initSolver(ProblemManager &p) override;
  bool solve(std::vector<Var> &setOfVar) override;
  bool solve() override;
  void uncheckedEnqueue(Lit l) override;
  bool varIsAssigned(Var v) override;
  bool getPolarity(Var v) override;
  bool decideAndComputeUnit(Lit l, std::vector<Lit> &units) override;
  bool failedLiteralProbing(Lit l) override;
  void whichAreUnits(std::vector<Var> &component,
                     std::vector<Lit> &units) override;
  void restart() override;
  void setAssumption(std::vector<Lit> &assums) override;
  std::vector<Lit> &getAssumption() override;
  void pushAssumption(Lit l) override;
  void popAssumption(unsigned count) override;
  void displayAssumption(std::ostream &out) override;
  void setNeedModel(bool b) override;
  void showTrail() override;
  std::vector<lbool> &getModel() override;
  lbool getModelVar(Var v) override;
  void getUnits(std::vector<Lit> &units) override;
  bool propagateAssumption() override;
  void addClause(const std::vector<Lit> &clause) override;

  double getActivity(Var v) override;
  double getCountConflict(Var v) override;
  void setCountConflict(Var v, double count) override;
  unsigned getNbConflict() override;
  void setReversePolarity(bool value) override;
  void decayCountConflict() override;
  bool isUnsat() override;

  void getCore() override;
  void getLastIUP(Lit l) override;

 private:
  void syncTrail();
  bool allProjAssigned();
};
}  // namespace d4
