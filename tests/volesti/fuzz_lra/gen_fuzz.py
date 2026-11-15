#!/usr/bin/env python3
"""Seeded generator of tiny satisfiable QF_LRA formulas that exercise SMT-LIB syntax.

Every formula has 1-4 Real variables, is conjoined with a bounding box (so its
solution set is bounded) and is kept only when a Monte-Carlo estimate with the
QF_LRA evaluator of run_union_decomp.py says that at least --min-frac of the
box is covered (so the volume is positive and the MC reference is accurate).

Each file starts with a `; fuzz-meta: {json}` line: index, seed, variables
(SMT-LIB symbols as written in the declarations), the bounding box (exact
rationals as strings), the syntax features the formula actually uses, and the
generator's MC volume estimate.

Usage: python3 tests/volesti/fuzz_lra/gen_fuzz.py [--count 500] [--seed 1]
           [--out tests/volesti/fuzz_lra/gen] [--min-frac 0.03]
"""
import argparse
import json
import math
import os
import random
import re
import sys
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))          # tests/volesti

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np                                  # noqa: E402
import run_union_decomp as rud                      # noqa: E402

# ---------------------------------------------------------------------------
# Evaluator: run_union_decomp.Formula + SMT-LIB scoping (push/pop), define-fun
# with arguments (Real- or Bool-valued), :named references, Real ite.
# ---------------------------------------------------------------------------


class FuzzFormula(rud.Formula):
    def __init__(self, text):                       # noqa: D401 (replaces base init)
        levels = [[]]
        for cmd in rud.parse_sexprs(text):
            if not isinstance(cmd, list) or not cmd:
                continue
            head = cmd[0]
            if head in ("push", "pop"):
                k = int(cmd[1]) if len(cmd) > 1 else 1
                for _ in range(k):
                    if head == "push":
                        levels.append([])
                    else:
                        if len(levels) == 1:
                            raise ValueError("pop below level 0")
                        levels.pop()
            elif head == "exit":
                break
            elif head in ("declare-const", "declare-fun", "define-fun", "assert"):
                levels[-1].append(cmd)
        self.reals, self.bools, self.globals, self.funs = [], [], {}, {}
        active = [c for lvl in levels for c in lvl]
        for cmd in active:
            if cmd[0] in ("declare-const", "declare-fun"):
                sort = cmd[2] if cmd[0] == "declare-const" else cmd[3]
                if cmd[0] == "declare-fun" and cmd[2]:
                    raise ValueError("uninterpreted functions are not supported")
                (self.reals if sort == "Real" else self.bools).append(cmd[1])
        self.n = len(self.reals)
        self.index = {v: i for i, v in enumerate(self.reals)}
        roots = []
        for cmd in active:
            if cmd[0] == "define-fun":
                if cmd[2]:
                    self.funs[cmd[1]] = ([p[0] for p in cmd[2]], cmd[4], dict(self.globals))
                else:
                    self.globals[cmd[1]] = self.compile(cmd[4], self.globals)
            elif cmd[0] == "assert":
                roots.append(self.compile(cmd[1], self.globals))
        self.root = rud.b_and(roots) if roots else rud.B("const", val=True)

    def compile(self, e, env):
        if isinstance(e, list) and e and e[0] == "!":
            out = self.compile(e[1], env)
            for i in range(2, len(e) - 1):
                if e[i] == ":named":
                    self.globals[e[i + 1]] = out    # later asserts may refer to it
            return out
        if isinstance(e, list) and e and isinstance(e[0], str) and e[0] in self.funs \
                and e[0] not in env:
            params, body, genv = self.funs[e[0]]
            xs = [self.compile(a, env) for a in e[1:]]
            new = dict(genv)
            new.update(zip(params, xs))
            return self.compile(body, new)
        return super().compile(e, env)


def mc_volume(text, box, n_points, seed, chunk=20000):
    """(volume, standard error, fraction) of F inside `box` by uniform sampling."""
    F = FuzzFormula(text)
    lo = np.array([float(a) for a, _ in box])
    hi = np.array([float(b) for _, b in box])
    if F.n != len(box):
        raise ValueError("evaluator sees %d Real variables, box has %d" % (F.n, len(box)))
    rng = np.random.default_rng(seed)
    hits = 0
    done = 0
    while done < n_points:
        m = min(chunk, n_points - done)
        X = lo + (hi - lo) * rng.random((m, F.n))
        hits += int(F.evaluate(X).sum())
        done += m
    p = hits / n_points
    bv = float(np.prod(hi - lo))
    se = bv * math.sqrt(max(p * (1 - p), 1.0 / n_points) / n_points)
    return bv * p, se, p


# ---------------------------------------------------------------------------
# Generator
# ---------------------------------------------------------------------------
QUOTED_NAMES = ["|x 0|", "|v(1)|", "|a;b|", "|#y|", "|µ z|".encode("ascii", "replace").decode(),
                "|x)(|", "|sp ace|", "|tab\tq|", "|p\"q|"]
SPECIAL_NAMES = ["x.k", "y$2", "z->w", "<w>", "a+b", "x!0", "_u", "q?r", "~t", "%m"]

# Optional features and their probability of being switched on per formula.
OPTIONAL = {
    "hdr_set_info": .35, "hdr_set_option": .35, "option_after_logic": .1,
    "print_success": .1, "no_set_logic": .08, "string_literal_info": .15,
    "quoted_symbols": .25, "quoted_alias": .15, "special_symbols": .2,
    "multi_assert": .6, "comments": .3, "comment_tricky": .15,
    "num_int": .5, "num_int_as_decimal": .4, "num_decimal": .5, "rational_div": .35,
    "neg_unary_const": .35, "large_consts": .1, "small_consts": .1,
    "mul_const_left": .5, "mul_const_right": .4, "div_var": .2, "unary_minus_var": .3,
    "nary_plus": .5, "nary_minus": .3, "const_lhs": .3, "const_moved": .2,
    "chain_cmp": .3, "chain_vars": .15, "bbox_chain": .4, "bbox_not": .15,
    "and": .5, "or": .5, "not": .4, "nested_not": .2, "implies": .3, "implies_nary": .1,
    "xor": .2, "bool_ite": .25, "bool_eq": .15, "bool_eq_nary": .1,
    "real_ite": .15, "let": .3, "let_shadow": .2, "let_bool": .2, "let_parallel": .1,
    "define0_real": .1, "define0_bool": .1, "define_args_real": .1,
    "define_args_bool": .1, "define_param_shadow": .08, "define_nested": .08,
    "named": .25, "named_ref": .1,
    "eq_nary": .2, "eq_nary_neg": .1, "distinct_nary": .2, "distinct2": .15,
    "const_atoms": .25, "true_false": .25, "dup_assert": .2, "dup_atom": .15,
    "redundant_atom": .2,
    "push_pop": .2, "push_declare": .1, "push_unclosed": .1, "push_n": .08,
    "mid_check_sat": .1, "get_model": .3, "get_model_noopt": .08, "get_value": .1,
    "exit": .4, "no_check_sat": .05, "to_real": .08, "to_real_in_qflra": .03,
}
FEATURES = sorted(OPTIONAL) + ["decl_const", "decl_fun"]


def is_dyadic_decimal(q):
    d = q.denominator
    while d % 2 == 0:
        d //= 2
    while d % 5 == 0:
        d //= 5
    return d == 1


def decimal_str(q):
    """Exact decimal rendering of a non-negative Fraction with denominator 2^a 5^b."""
    k = 0
    while (q * 10 ** k).denominator != 1:
        k += 1
    k = max(k, 1)
    n = int(q * 10 ** k)
    s = str(n).rjust(k + 1, "0")
    return s[:-k] + "." + s[-k:]


class Gen:
    def __init__(self, idx, seed, dims=None):
        self.idx = idx
        self.rng = random.Random(seed)
        self.seed = seed
        r = self.rng
        self.on = {f for f, p in OPTIONAL.items() if r.random() < p}
        self.used = set()
        self.n = dims or r.randint(1, 4)
        self.scale = Fraction(1)
        if "large_consts" in self.on:
            self.scale = Fraction(10) ** r.choice([4, 6])
        elif "small_consts" in self.on:
            self.scale = Fraction(1, 10 ** r.choice([3, 5]))
        if self.scale != 1:
            self.use("large_consts" if self.scale > 1 else "small_consts")
        self.names, self.refs, self.decls = [], [], []
        pool_q = list(QUOTED_NAMES)
        pool_s = list(SPECIAL_NAMES)
        r.shuffle(pool_q)
        r.shuffle(pool_s)
        for i in range(self.n):
            name = "x%d" % i
            decl = ref = name
            c = r.random()
            if self.feat("quoted_symbols", .5) and pool_q:
                decl = ref = pool_q.pop()
            elif self.feat("special_symbols", .5) and pool_s:
                decl = ref = pool_s.pop()
            elif self.feat("quoted_alias", .6):
                if c < .5:
                    decl, ref = "|%s|" % name, name
                else:
                    decl, ref = name, "|%s|" % name
            self.names.append(decl)
            self.refs.append(ref)
        self.box = []
        for i in range(self.n):
            lo = Fraction(r.randint(-40, 40), 10)
            w = Fraction(r.randint(5, 40), 10)
            self.box.append((lo * self.scale, (lo + w) * self.scale))
        self.preamble, self.defs, self.cmds = [], [], []
        self.counter = 0

    # -- helpers ------------------------------------------------------------
    def use(self, f):
        self.used.add(f)

    def feat(self, f, p=1.0):
        if f in self.on and self.rng.random() < p:
            self.used.add(f)
            return True
        return False

    def fresh(self, base):
        self.counter += 1
        return "%s%d" % (base, self.counter)

    def point(self):
        """A random point in the middle 80% of the box (exact rationals)."""
        return [lo + (hi - lo) * Fraction(self.rng.randint(1, 9), 10) for lo, hi in self.box]

    def const(self, q):
        q = Fraction(q)
        neg = q < 0
        a = -q if neg else q
        r = self.rng
        if a.denominator == 1:
            if ("to_real" in self.on or "to_real_in_qflra" in self.on) and r.random() < .4:
                self.use("to_real")
                s = "(to_real %d)" % a.numerator
            elif self.feat("num_int_as_decimal", .5):
                s = str(a.numerator) + ".0"
            else:
                self.use("num_int")
                s = str(a.numerator)
        elif is_dyadic_decimal(a) and not self.feat("rational_div", .4):
            self.use("num_decimal")
            s = decimal_str(a)
        else:
            self.use("rational_div")
            num = str(a.numerator) + (".0" if self.feat("num_int_as_decimal", .3) else "")
            s = "(/ %s %d)" % (num, a.denominator)
        if neg:
            if a.denominator != 1 and s.startswith("(/") and r.random() < .5:
                s = "(/ (- %d) %d)" % (a.numerator, a.denominator)
            else:
                s = "(- %s)" % s
            self.use("neg_unary_const")
        return s

    def term(self, i, k):
        """k * x_i."""
        x = self.refs[i]
        r = self.rng
        if k == 1:
            return x
        if k == -1 and self.feat("unary_minus_var", .7):
            return "(- %s)" % x
        if self.feat("div_var", .3) and k.numerator == 1 and k > 0:
            return "(/ %s %s)" % (x, self.const(Fraction(k.denominator)))
        if self.feat("mul_const_right", .5):
            return "(* %s %s)" % (x, self.const(k))
        self.use("mul_const_left")
        return "(* %s %s)" % (self.const(k), x)

    def linexpr(self, coefs, const=Fraction(0)):
        """Render sum_i coefs[i] x_i + const (coefs: dict i -> Fraction)."""
        items = [(i, Fraction(k)) for i, k in coefs.items() if k != 0]
        self.rng.shuffle(items)
        if const != 0:
            items.append((None, Fraction(const)))
        if not items:
            return self.const(0)

        def piece(i, k):
            return self.const(k) if i is None else self.term(i, k)

        if len(items) == 1:
            return piece(*items[0])
        if self.feat("nary_minus", .4):
            first = piece(*items[0])
            rest = [piece(i, -k) for i, k in items[1:]]
            return "(- %s %s)" % (first, " ".join(rest))
        if len(items) > 2 or self.feat("nary_plus", .8):
            if len(items) > 2:
                self.use("nary_plus")
            return "(+ %s)" % " ".join(piece(i, k) for i, k in items)
        return "(+ %s %s)" % (piece(*items[0]), piece(*items[1]))

    def rand_coefs(self, nvars=None):
        r = self.rng
        k = nvars or r.randint(1, self.n)
        idx = r.sample(range(self.n), k)
        coefs = {}
        for i in idx:
            c = Fraction(r.choice([-3, -2, -1, 1, 1, 2, 3]))
            if r.random() < .2:
                c /= r.choice([2, 3, 4])
            coefs[i] = c
        return coefs

    # -- atoms ----------------------------------------------------------------
    def atom(self, depth=0):
        """A half-space atom through a random interior point (sometimes fancier)."""
        r = self.rng
        p = self.point()
        coefs = self.rand_coefs()
        c = sum(k * p[i] for i, k in coefs.items())
        op = r.choice(["<", "<=", ">", ">="])
        if self.feat("chain_vars", .15) and self.n >= 2:
            i, j = r.sample(range(self.n), 2)
            return "(%s %s %s)" % (r.choice(["<", "<="]), self.refs[i], self.refs[j])
        if self.feat("chain_cmp", .3):
            # c1 < e < c2 with c1 < c < c2 (half-width relative to the box scale).
            w = sum(abs(k) * (hi - lo) for (i, k), (lo, hi) in
                    zip(coefs.items(), [self.box[i] for i in coefs])) * Fraction(r.randint(1, 4), 10)
            o1, o2 = r.choice([("<", "<"), ("<=", "<="), ("<", "<=")])
            if r.random() < .5:
                return "(%s %s %s %s)" % (o1, self.const(c - w), self.linexpr(coefs),
                                          self.const(c + w))
            return "(%s %s %s %s)" % (o1.replace("<", ">"), self.const(c + w),
                                      self.linexpr(coefs), self.const(c - w))
        if self.feat("real_ite", .25) and depth < 2:
            # Both branches pass through p, so the atom keeps p on its boundary.
            coefs2 = self.rand_coefs()
            c2 = sum(k * p[i] for i, k in coefs2.items())
            cond = self.atom(depth + 1)
            e1 = self.linexpr(coefs, -c)
            e2 = self.linexpr(coefs2, -c2)
            return "(%s (ite %s %s %s) %s)" % (op, cond, e1, e2, self.const(0))
        if self.feat("define_args_real", .3) and len(coefs) <= 3:
            vs = sorted(coefs)
            fname = self.fresh("lf")
            params = ["u%d" % j for j in range(len(vs))]
            if self.feat("define_param_shadow", .5):
                # Parameters named like the global variables (in reverse order).
                params = [self.refs[v] for v in reversed(vs)]
            body = self.linexpr_named(params, [coefs[v] for v in vs])
            if self.feat("define_nested", .4):
                inner = self.fresh("id")
                self.defs.append("(define-fun %s ((t Real)) Real t)" % inner)
                body = "(%s %s)" % (inner, body)
            self.defs.append("(define-fun %s (%s) Real %s)" % (
                fname, " ".join("(%s Real)" % q for q in params), body))
            return "(%s (%s %s) %s)" % (op, fname, " ".join(self.refs[v] for v in vs),
                                        self.const(c))
        if self.feat("let", .3):
            v = self.fresh("s")
            lhs = self.linexpr(coefs)
            if self.feat("let_shadow", .5):
                return "(let ((%s %s)) (let ((%s (+ %s 1))) (%s %s %s)))" % (
                    v, lhs, v, v, op, v, self.const(c + 1))
            return "(let ((%s %s)) (%s %s %s))" % (v, lhs, op, v, self.const(c))
        if self.feat("let_parallel", .4) and self.n >= 2:
            i, j = r.sample(range(self.n), 2)
            # Inside (let ((xi xj) (xj xi)) body) the name xi denotes x_j.
            saved = list(self.refs)
            self.refs[i], self.refs[j] = saved[j], saved[i]
            body = "(%s %s %s)" % (op, self.linexpr(coefs), self.const(c))
            self.refs = saved
            return "(let ((%s %s) (%s %s)) %s)" % (saved[i], saved[j], saved[j], saved[i], body)
        if self.feat("define0_real", .3):
            k = self.fresh("K")
            self.defs.append("(define-fun %s () Real %s)" % (k, self.const(c)))
            rhs = k
        else:
            rhs = self.const(c)
        if self.feat("const_moved", .3):
            return "(%s %s %s)" % (op, self.linexpr(coefs, -c), self.const(0))
        if self.feat("const_lhs", .3):
            flip = {"<": ">", "<=": ">=", ">": "<", ">=": "<="}[op]
            return "(%s %s %s)" % (flip, rhs, self.linexpr(coefs))
        return "(%s %s %s)" % (op, self.linexpr(coefs), rhs)

    def linexpr_named(self, names, ks):
        saved = self.refs
        self.refs = names
        out = self.linexpr({j: k for j, k in enumerate(ks)})
        self.refs = saved
        return out

    def const_atom(self, value):
        """A variable-free atom with the given truth value."""
        r = self.rng
        if self.feat("true_false", .5):
            return "true" if value else "false"
        self.use("const_atoms")
        a, b = Fraction(r.randint(0, 9)), Fraction(r.randint(10, 19), 2)
        if a == b:
            b += 1
        lo, hi = min(a, b), max(a, b)
        return ("(< %s %s)" if value else "(> %s %s)") % (self.const(lo), self.const(hi))

    def bexpr(self, depth=0):
        r = self.rng
        if depth >= 2 or r.random() < .55:
            a = self.atom()
            if self.feat("nested_not", .2):
                return "(not (not %s))" % a
            if self.feat("not", .2):
                return "(not %s)" % a
            return a
        ops = [o for o in ("and", "or", "implies", "xor", "bool_ite", "bool_eq") if o in self.on]
        if not ops:
            return self.atom()
        op = r.choice(ops)
        self.use(op)
        k = 2 if r.random() < .7 else 3
        sub = [self.bexpr(depth + 1) for _ in range(k)]
        if op in ("and", "or"):
            if self.feat("const_atoms", .3) or self.feat("true_false", .3):
                sub.append(self.const_atom(op == "and"))
            if self.feat("dup_atom", .3):
                sub.append(sub[0])
            return "(%s %s)" % (op, " ".join(sub))
        if op == "implies":
            if k == 3:
                self.use("implies_nary") if "implies_nary" in self.on else sub.pop()
            return "(=> %s)" % " ".join(sub)
        if op == "xor":
            return "(xor %s %s)" % (sub[0], sub[1])
        if op == "bool_ite":
            while len(sub) < 3:
                sub.append(self.bexpr(depth + 1))
            return "(ite %s)" % " ".join(sub)
        if k == 3 and "bool_eq_nary" in self.on:
            self.use("bool_eq_nary")
            return "(= %s)" % " ".join(sub)
        return "(= %s %s)" % (sub[0], sub[1])

    # -- whole formula ----------------------------------------------------------
    def bbox_constraints(self):
        out = []
        r = self.rng
        argfun = None
        if self.feat("define_args_bool", .7):
            argfun = self.fresh("inb")
            self.defs.append("(define-fun %s ((u Real) (lo Real) (hi Real)) Bool "
                             "(and (<= lo u) (<= u hi)))" % argfun)
        for i, (lo, hi) in enumerate(self.box):
            x = self.refs[i]
            if argfun and r.random() < .6:
                out.append("(%s %s %s %s)" % (argfun, x, self.const(lo), self.const(hi)))
            elif self.feat("bbox_chain", .6):
                o = r.choice(["<", "<="])
                if r.random() < .5:
                    out.append("(%s %s %s %s)" % (o, self.const(lo), x, self.const(hi)))
                else:
                    out.append("(%s %s %s %s)" % (o.replace("<", ">"), self.const(hi), x,
                                                  self.const(lo)))
            elif self.feat("bbox_not", .5):
                out.append("(not (< %s %s))" % (x, self.const(lo)))
                out.append("(not (> %s %s))" % (x, self.const(hi)))
            else:
                out.append("(>= %s %s)" % (x, self.const(lo)))
                out.append("(<= %s %s)" % (x, self.const(hi)))
        return out

    def comment(self):
        if self.feat("comment_tricky", .5):
            return self.rng.choice(["; ) unbalanced ( paren | pipe \" quote",
                                    ";(assert false)", ";; |unterminated",
                                    "; (check-sat) (exit)"])
        self.use("comments")
        return "; " + self.rng.choice(["constraint", "box", "cut", "x < y"])

    def build(self, ncuts):
        r = self.rng
        lines = []
        logic = "(set-logic QF_LRA)"
        if self.feat("hdr_set_option", 1):
            if self.feat("print_success", .5):
                lines.append("(set-option :print-success true)")
            lines.append("(set-option :produce-models true)")
        if self.feat("hdr_set_info", 1):
            lines.append("(set-info :smt-lib-version 2.6)")
            if self.feat("string_literal_info", 1):
                lines.append('(set-info :source "fuzz) ""generated"" (seed %d) ; not a comment")'
                             % self.seed)
            else:
                lines.append("(set-info :source |fuzz generator\n multi-line ( source|)")
            lines.append("(set-info :category \"crafted\")")
            lines.append("(set-info :status sat)")
        if not self.feat("no_set_logic", 1):
            lines.append(logic)
        if self.feat("option_after_logic", 1):
            lines.append("(set-option :produce-models true)")
        for i in range(self.n):
            if r.random() < .5:
                self.use("decl_const")
                lines.append("(declare-const %s Real)" % self.names[i])
            else:
                self.use("decl_fun")
                lines.append("(declare-fun %s () Real)" % self.names[i])
        body = []                           # list of assertion bodies (conjuncts)
        bb = self.bbox_constraints()
        for c in range(ncuts):
            body.append(self.bexpr())
        if self.feat("define0_bool", 1):
            b = self.fresh("B")
            self.defs.append("(define-fun %s () Bool %s)" % (b, self.bexpr(1)))
            body.append("(or %s %s)" % (b, self.atom()) if r.random() < .5 else b)
        if self.feat("let_bool", 1):
            b1, b2 = self.fresh("b"), self.fresh("b")
            body.append("(let ((%s %s) (%s %s)) (or %s %s))" % (b1, self.atom(), b2, self.atom(),
                                                               b1, b2))
        if self.feat("eq_nary", 1):
            vs = [self.refs[i] for i in range(self.n)] + [self.const(self.point()[0])]
            while len(vs) < 3:
                vs.append(vs[-1])
            body.append("(or (= %s) %s)" % (" ".join(vs), self.atom()))
        if self.feat("eq_nary_neg", 1):
            vs = [self.refs[i] for i in range(self.n)] + [self.const(self.box[0][0])]
            if len(vs) < 3:
                vs.append(self.const(self.box[0][1]))
            body.append("(not (= %s))" % " ".join(vs))
        if self.feat("distinct_nary", 1):
            vs = [self.refs[i] for i in range(self.n)] + [self.const(self.point()[0])]
            if len(vs) < 3:
                vs.append(self.const(self.box[0][1] + 1))
            r.shuffle(vs)
            body.append("(distinct %s)" % " ".join(vs))
        if self.feat("distinct2", 1):
            p = self.point()
            coefs = self.rand_coefs()
            body.append("(distinct %s %s)" % (self.linexpr(coefs),
                                              self.const(sum(k * p[i] for i, k in coefs.items()))))
        if self.feat("redundant_atom", 1):
            i = r.randrange(self.n)
            body.append("(< %s %s)" % (self.refs[i], self.const(self.box[i][1] + self.scale)))
        if self.feat("const_atoms", .5) or self.feat("true_false", .5):
            body.append(self.const_atom(True))
        # Assemble asserts.
        if self.feat("dup_assert", 1):
            body.append(r.choice(bb + body))
        named = []
        if self.feat("named", 1):
            k = r.randrange(len(body)) if body else None
            if k is not None:
                nm = self.fresh("n")
                body[k] = "(! %s :named %s)" % (body[k], nm)
                named.append(nm)
        conj = bb + body
        r.shuffle(conj)
        asserts = []
        if self.feat("multi_assert", 1) or len(conj) == 1:
            groups = [conj]
            if len(conj) > 1:
                cut = sorted(r.sample(range(1, len(conj)), r.randint(1, len(conj) - 1)))
                groups = [conj[a:b] for a, b in zip([0] + cut, cut + [len(conj)])]
        else:
            groups = [conj]
        for g in groups:
            asserts.append("(assert %s)" % (g[0] if len(g) == 1 else "(and %s)" % "\n    ".join(g)))
        if named and self.feat("named_ref", 1):
            asserts.append("(assert %s)" % named[0])
        # Scoping junk.
        if self.feat("push_pop", 1):
            pos = r.randrange(len(asserts) + 1)
            junk = ["(push 1)", "(assert false)" if r.random() < .5 else
                    "(assert (< %s %s))" % (self.refs[0], self.const(self.box[0][0] - 1)), "(pop 1)"]
            asserts[pos:pos] = junk
        if self.feat("push_n", 1):
            asserts[0:0] = ["(push 2)", "(assert (> %s %s))" % (self.refs[0],
                                                                  self.const(self.box[0][1] + 5)),
                            "(pop 2)"]
        if self.feat("push_declare", 1):
            asserts += ["(push 1)", "(declare-const tmpv Real)",
                        "(assert (and (< tmpv 0) (> tmpv 1)))", "(pop 1)"]
        if self.feat("push_unclosed", 1):
            asserts.insert(r.randrange(len(asserts) + 1), "(push 1)")
        if self.feat("mid_check_sat", 1) and len(asserts) > 1:
            asserts.insert(r.randrange(1, len(asserts)), "(check-sat)")
        lines += self.defs
        lines += asserts
        if not self.feat("no_check_sat", 1):
            lines.append("(check-sat)")
            if self.feat("get_model_noopt", 1):
                lines = [l for l in lines if l != "(set-option :produce-models true)"]
                lines.append("(get-model)")
            elif "(set-option :produce-models true)" in lines:
                if self.feat("get_model", 1):
                    lines.append("(get-model)")
                if self.feat("get_value", 1):
                    lines.append("(get-value (%s))" % self.refs[0])
        if self.feat("exit", 1):
            lines.append("(exit)")
        if "comments" in self.on or "comment_tricky" in self.on:
            for _ in range(r.randint(1, 3)):
                pos = r.randrange(len(lines) + 1)
                lines.insert(pos, self.comment())
            # an inline comment inside a multi-line assert
            for k, l in enumerate(lines):
                if "\n    " in l and r.random() < .5:
                    lines[k] = l.replace("\n    ", "  " + self.comment() + "\n    ", 1)
                    break
        if "to_real" in self.used:
            if "to_real_in_qflra" in self.on:
                self.use("to_real_in_qflra")
            else:
                lines = [("(set-logic QF_LIRA)" if l == logic else l) for l in lines]
        return "\n".join(lines) + "\n"


def generate_one(idx, seed, min_frac, mc_points=20000, dims=None):
    for attempt in range(40):
        g = Gen(idx, seed * 1000003 + idx * 101 + attempt, dims)
        ncuts = g.rng.randint(1, 2) if attempt < 20 else 1
        text = g.build(ncuts)
        try:
            vol, se, frac = mc_volume(text, g.box, mc_points, seed=idx)
        except Exception as e:  # generator/evaluator mismatch: report, retry
            sys.stderr.write("gen %d attempt %d: evaluator error %r\n" % (idx, attempt, e))
            continue
        if frac >= min_frac:
            meta = {"index": idx, "seed": g.seed, "vars": g.names, "refs": g.refs,
                    "box": [[str(a), str(b)] for a, b in g.box], "features": sorted(g.used),
                    "mc_volume": vol, "mc_frac": frac}
            return "; fuzz-meta: %s\n%s" % (json.dumps(meta), text), meta
    raise RuntimeError("could not generate formula %d" % idx)


def read_meta(text):
    m = re.search(r"^; fuzz-meta: (\{.*\})$", text, re.M)
    if not m:
        return None
    meta = json.loads(m.group(1))
    meta["box"] = [(Fraction(a), Fraction(b)) for a, b in meta["box"]]
    return meta


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--count", type=int, default=500)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--start", type=int, default=0)
    ap.add_argument("--min-frac", type=float, default=0.03)
    ap.add_argument("--out", default=os.path.join(HERE, "gen"))
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    stats = {}
    for idx in range(args.start, args.start + args.count):
        text, meta = generate_one(idx, args.seed, args.min_frac)
        with open(os.path.join(args.out, "fuzz_%04d.smt2" % idx), "w") as f:
            f.write(text)
        for ft in meta["features"]:
            stats[ft] = stats.get(ft, 0) + 1
    print("wrote %d formulas to %s" % (args.count, args.out))
    missing = [f for f in FEATURES if f not in stats]
    print("feature coverage: %s" % ", ".join("%s=%d" % kv for kv in sorted(stats.items())))
    if missing:
        print("features never used: %s" % ", ".join(missing))


if __name__ == "__main__":
    main()
