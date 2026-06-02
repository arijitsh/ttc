  #include "ResidualSimplifier.hpp"

  #include <chrono>
  #include <unordered_set>

  #include "../../../../ddnnf/profiler.hpp"
  #include "../../../../features.hpp"

  namespace d4 {

  ResidualSimplifier::ResidualSimplifier(cvc5::Solver& solver)
      : d_solver(solver),
        d_trueTerm(ttc::getTermBuilder(solver).mkBoolean(true)),
        d_falseTerm(ttc::getTermBuilder(solver).mkBoolean(false)) {}

  bool ResidualSimplifier::extractEqualityAssignment(const cvc5::Term& varCandidate,
                                                    const cvc5::Term& valueCandidate,
                                                    cvc5::Term& term,
                                                    cvc5::Term& value) const {
    const cvc5::Sort sort = varCandidate.getSort();
    if (!(sort.isInteger() || sort.isReal())) {
      return false;
    }

    if ((sort.isInteger() && varCandidate.isIntegerValue()) ||
        (sort.isReal() && varCandidate.isRealValue())) {
      return false;
    }

    const cvc5::Sort valueSort = valueCandidate.getSort();
    if (valueSort != sort) {
      if (!(sort.isReal() && valueSort.isInteger())) {
        return false;
      }
    }

    if (sort.isInteger()) {
      if (!valueCandidate.isIntegerValue()) {
        return false;
      }
    } else {
      if (!(valueCandidate.isRealValue() || valueCandidate.isIntegerValue())) {
        return false;
      }
    }

    if (valueCandidate.isIntegerValue() || valueCandidate.isRealValue()) {
      term = varCandidate;
      value = valueCandidate;
      return true;
    }

    return false;
  }

  bool ResidualSimplifier::extractAssignment(const cvc5::Term& literal,
                                            cvc5::Term& term,
                                            cvc5::Term& value) const {
    Trace("extracta") << "c Trying to extract assignment from literal (Sort:" << literal.getSort()
                      << " , Kind: " << literal.getKind()  << ") \n";
    if (literal.getNumChildren() == 2)
      Trace("extracta") << literal[0].getKind() << " -> " << literal[1].getKind()
                        << "  \n\t" << literal
                      << std::endl;
    if (literal.getSort().isBoolean() )  {
      if (literal.isBooleanValue()) {
        Trace("extracta") << "c Extracted assignment Z " << literal << " -> " << literal.getBooleanValue() << "\n";
        term = literal;
        value = d_trueTerm;
        return true;
      }
      if (literal.getKind() == cvc5::Kind::EQUAL ) {
        auto leftKind = literal[0].getKind();
        auto rightKind = literal[1].getKind();
        if (leftKind == cvc5::Kind::CONSTANT && (rightKind == cvc5::Kind::CONST_RATIONAL || rightKind == cvc5::Kind::CONST_INTEGER)) {
          Trace("extracta") << "c Extracted assignment Y " << literal[0] << " -> " << literal[1] << "\n";
          term = literal[0]  ;
          value = literal[1];
          return true;
        }
        Trace("extracta") << "c [fail] Left or right side is not a variable or constant\n";
        return false;
      }
      if (literal.getKind() == cvc5::Kind::CONSTANT) {
        Trace("extracta") << "c Extracted assignment W " << literal << " -> " << d_trueTerm << "\n";
        term = literal;
        value = d_trueTerm;
        return true;
      }

      if (literal.getKind() == cvc5::Kind::NOT) {
        Trace("extracta") << "c [fail] NOT literal does not have exactly one child\n";
        return false;
      }

      Trace("extracta") << "c [fail] Literal is not a variable, equality, or NOT\n";
      return false;




      if (literal.getKind() == cvc5::Kind::NOT && literal.getNumChildren() == 1) {
        const cvc5::Term& child = literal[0];
        if (!child.getSort().isBoolean()) {
          Trace("extracta") << "c [fail] Child is not boolean\n";
          return false;
        }
        term = child;
        value = d_falseTerm;
        Trace("extracta") << "c Extracted assignment A " << term << " -> " << value << "\n";
        return true;
      }

      term = literal;
      value = d_trueTerm;
      Trace("extracta") << "c Extracted assignment B " << term << " -> " << value << "\n";
      Trace("extracta") << "c Literal is a boolean variable\n";
      if (term.getKind() == cvc5::Kind::EQUAL && term[0].getNumChildren() == 1 && term[1].getNumChildren() == 1) {
        Trace("extracta") << "c Extracting equality assignment from equality term\n";
        auto literal = term;
        return extractAssignment(literal, term, value);
      }
      return true;
    }
    Trace("extracta") << "c Literal is not boolean\n";

    if (literal.getKind() == cvc5::Kind::EQUAL && literal.getNumChildren() == 2) {
      const cvc5::Term& lhs = literal[0];
      const cvc5::Term& rhs = literal[1];
      if (extractEqualityAssignment(lhs, rhs, term, value) ||
          extractEqualityAssignment(rhs, lhs, term, value)) {
        Trace("extracta") << "c Extracted assignment C " << term << " -> " << value << "\n";
        return true;
      }
    }
    Trace("extracta") << "c Failed to extract assignment from literal " << literal << "\n";
    return false;
  }

  bool ResidualSimplifier::updateTrailSubstitutions() {
    const cvc5::modes::LearnedLitType types[] = {
        cvc5::modes::LearnedLitType::PREPROCESS_SOLVED,
        cvc5::modes::LearnedLitType::PREPROCESS,
        cvc5::modes::LearnedLitType::INPUT,
        cvc5::modes::LearnedLitType::SOLVABLE,
        cvc5::modes::LearnedLitType::CONSTANT_PROP,
        cvc5::modes::LearnedLitType::INTERNAL};

    std::unordered_set<cvc5::Term> seen;
    bool updated = false;
    bool foundAny = false;

    auto stream = Trace("extract");
    stream << "c Extracting assignments from learned literals of types ";

    for (const auto type : types) {
      d_solver.checkSat();
      std::vector<cvc5::Term> literals = d_solver.getLearnedLiterals(type);
      for (const auto& literal : literals) {
        if (!seen.insert(literal).second) {
          continue;
        }
        cvc5::Term term;
        cvc5::Term value;
        if (!extractAssignment(literal, term, value)) {
          Trace("extract") << "c Failed to extract assignment from literal " << literal << "\n";
          continue;
        }
        stream << "[" << literal << ", " << type << "] " << term << " -> " << value << "\n";
        foundAny = true;
        auto [it, inserted] = d_termToIndex.emplace(term, d_substTerms.size());
        if (inserted) {
          d_substTerms.push_back(term);
          d_substValues.push_back(value);
          updated = true;
        } else if (d_substValues[it->second] != value) {
          d_substValues[it->second] = value;
          updated = true;
        }
      }
    }
    stream << std::endl;

    if (!foundAny && !d_substTerms.empty()) {
      d_termToIndex.clear();
      d_substTerms.clear();
      d_substValues.clear();
      updated = true;
    }

    return updated;
  }

  void ResidualSimplifier::flattenAnd(const cvc5::Term& t,
                                      std::vector<cvc5::Term>& out) {
    if (t.getKind() == cvc5::Kind::AND) {
      for (size_t i = 0, n = t.getNumChildren(); i < n; ++i) {
        flattenAnd(t[i], out);
      }
      return;
    }
    if (t.isBooleanValue() && t.getBooleanValue()) return;
    out.push_back(t);
  }

  SimplifyResult ResidualSimplifier::simplify(
      const std::vector<cvc5::Term>& assertions) {
    SimplifyResult result;
    if (assertions.empty()) {
      return result;
    }

    Trace("simplify") << "c Residual SMT assertions before simplification:\n";
    for (const auto& a : assertions) {
      Trace("simplify") << "c   " << a << "\n";
    }

    updateTrailSubstitutions();

    double substitutionSeconds = 0.0;
    std::vector<cvc5::Term> substitutedAssertions;
    substitutedAssertions.reserve(assertions.size());
    auto stream = Trace("substitute");
    stream << "c Residual SMT assertions after substitution of "
          << d_substTerms << "\n";



    for (const auto& a : assertions) {
      cvc5::Term substituted = a;
      if (!d_substTerms.empty()) {
        const auto substStart = std::chrono::steady_clock::now();
        substituted = substituted.substitute(d_substTerms, d_substValues);
        const auto substEnd = std::chrono::steady_clock::now();
        substitutionSeconds +=
            std::chrono::duration_cast<std::chrono::duration<double>>(substEnd -
                                                                      substStart)
                .count();
        stream << a << " => " << substituted
              << std::endl;
      }
      // stream << std::endl;

      if (substituted.isBooleanValue()) {
        if (!substituted.getBooleanValue()) {
          result.unsat = true;
          substitutedAssertions.clear();
          break;
        }
        continue;
      }

      substitutedAssertions.push_back(substituted);
    }

    Profile.addResidualSubstitute(substitutionSeconds);

    if (result.unsat) {
      Profile.addResidualSimplify(0.0);
      return result;
    }

    const auto simplifyStart = std::chrono::steady_clock::now();

    std::vector<cvc5::Term> simplifiedAssertions;
    simplifiedAssertions.reserve(substitutedAssertions.size());

    for (const auto& term : substitutedAssertions) {
      auto simp = d_solver.simplify(term);
      Trace("simplify") << term << " => " << simp << "\n";

      if (simp.isBooleanValue()) {
        if (!simp.getBooleanValue()) {
          result.unsat = true;
          simplifiedAssertions.clear();
          break;
        }
        continue;
      }

      simplifiedAssertions.push_back(simp);
    }

    const auto simplifyEnd = std::chrono::steady_clock::now();
    const double simplifySeconds =
        std::chrono::duration_cast<std::chrono::duration<double>>(simplifyEnd -
                                                                  simplifyStart)
            .count();
    Profile.addResidualSimplify(simplifySeconds);

    if (!result.unsat) {
      result.assertions.insert(result.assertions.end(),
                              simplifiedAssertions.begin(),
                              simplifiedAssertions.end());
    }

    return result;
  }

  }  // namespace d4
