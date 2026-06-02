#ifndef REWRITE_HPP
#define REWRITE_HPP

#include <cvc5/cvc5.h>
#include <vector>

#include "features.hpp"

inline cvc5::Term simpleRewrite(cvc5::Solver& solver, const cvc5::Term& term)
{
    using namespace cvc5;
    auto& tm = ttc::getTermBuilder(solver);
    switch (term.getKind())
    {
        case Kind::NOT:
        {
            Term child = simpleRewrite(solver, term[0]);
            if (child.isBooleanValue())
            {
                return tm.mkBoolean(!child.getBooleanValue());
            }
            return tm.mkTerm(Kind::NOT, {child});
        }
        case Kind::AND:
        {
            std::vector<Term> children;
            for (size_t i = 0, n = term.getNumChildren(); i < n; ++i)
            {
                Term c = simpleRewrite(solver, term[i]);
                if (c.isBooleanValue())
                {
                    if (!c.getBooleanValue())
                    {
                        return tm.mkBoolean(false);
                    }
                }
                else
                {
                    children.push_back(c);
                }
            }
            if (children.empty())
            {
                return tm.mkBoolean(true);
            }
            if (children.size() == 1)
            {
                return children[0];
            }
            return tm.mkTerm(Kind::AND, children);
        }
        case Kind::OR:
        {
            std::vector<Term> children;
            for (size_t i = 0, n = term.getNumChildren(); i < n; ++i)
            {
                Term c = simpleRewrite(solver, term[i]);
                if (c.isBooleanValue())
                {
                    if (c.getBooleanValue())
                    {
                        return tm.mkBoolean(true);
                    }
                }
                else
                {
                    children.push_back(c);
                }
            }
            if (children.empty())
            {
                return tm.mkBoolean(false);
            }
            if (children.size() == 1)
            {
                return children[0];
            }
            return tm.mkTerm(Kind::OR, children);
        }
        default:
            return term;
    }
}

#endif // REWRITE_HPP
