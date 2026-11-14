#!/usr/bin/env python3
"""T6 / component C1: is ttc's polytope decomposition exact?

thm:decomposition (KR25 manuscript, Alg. 1 lines 1-2): the AllSAT cover of the
Boolean abstraction of F gives cubes whose polytopes satisfy

    (completeness)  sol(F) = U_i P(cube_i)      (up to a measure-zero set)
    (soundness)     P(cube_i) is a subset of sol(F) for every i.

ttc writes every canonicalized polytope with `--dump-ine PREFIX`
(PREFIX_cubeN.ine, cdd H-format).  The .ine files carry no variable names: ttc
numbers its columns by first occurrence of each Real variable in a DFS of the
raw asserted term (src/parser.cpp `collect`), NOT by declaration order; the
evaluator reproduces that order.  This script runs ttc on each formula, stops it once phase 1a
(the dump) is done, and checks the decomposition against an independent
SMT-LIB QF_LRA evaluator written here:

  * points: uniform in a bounding box; hit-and-run samples inside every dumped
    P_i (soundness); hit-and-run samples inside every full-dimensional disjunct
    of an independent DNF of F (completeness); points at distance 1e-6..1e-2
    (relative to the box) on both sides of every facet of every P_i and of every
    atom hyperplane of F;
  * for every point x: F(x) and c(x) = #{i : x in P_i} (closed P_i).  A point
    with F(x) != (c(x) > 0) is a mismatch.  It is an INTERIOR mismatch if x is
    farther than tau = 1e-9 * max|bbox coordinate| (a floating-point guard) from
    every atom hyperplane of F and every facet hyperplane of every P_i; then F
    and the cover are both constant on the ball B(x, dist), so the error has
    positive measure -- a decomposition
    bug, no statistics needed.  Other mismatches are BOUNDARY (measure zero,
    e.g. strict vs non-strict) and only reported.  Interior mismatches are
    re-checked with an exact cvc5 query when a cvc5 binary is available.

Also reported: AllSAT cube count, dumped / zero-volume polytopes, cover
overlap ratio sum_i vol(P_i) / vol(sol F) (vol(P_i) exact via qhull when scipy
is available; vol(sol F) by box Monte Carlo, and by Karp-Luby over the P_i
samples), and the evaluator sanity check MC vol(F) vs `; exact-volume:` (z-test,
|z| < 4).

Header keys in the .smt2 inputs:  ; bbox: lo hi  (or lo1 hi1 lo2 hi2 ...),
; exact-volume: v, ; report-only: true.  Extra inputs (e.g. syntheticLRA
.smt2.xz files) can be passed on the command line; .xz is decompressed into
the work directory; a missing bbox is derived from the DNF of F (scipy LP) or
from the dumped polytopes.

Exit status is non-zero iff a non-report-only input has an interior mismatch,
ttc fails / does not run its volume engine, or the evaluator sanity check fails.

Usage: python3 tests/volesti/run_union_decomp.py [--ttc ./build/ttc]
          [--filter RE] [--seed S] [--points N] [--keep] [--workdir DIR]
          [--json OUT] [extra.smt2[.xz] ...]
Requires numpy; scipy is optional (LPs for the DNF regions / exact volumes).
"""
import argparse
import glob
import json
import lzma
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile
import threading
import time
from fractions import Fraction

# One compute thread (the user's CPU budget): set before numpy loads its BLAS.
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np  # noqa: E402

try:
    import scipy.optimize
    import scipy.spatial
    HAVE_SCIPY = True
except Exception:  # pragma: no cover
    HAVE_SCIPY = False

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
INPUT_DIR = os.path.join(HERE, "union_decomp")
DNF_CAP = 4096

# ---------------------------------------------------------------------------
# S-expressions
# ---------------------------------------------------------------------------
TOKEN_RE = re.compile(r'\s+|;[^\n]*|\(|\)|\|[^|]*\||"(?:[^"]|"")*"|[^\s();|"]+')


def parse_sexprs(text):
    stack = [[]]
    for m in TOKEN_RE.finditer(text):
        tok = m.group(0)
        if tok[0].isspace() or tok[0] == ';':
            continue
        if tok == '(':
            stack.append([])
        elif tok == ')':
            if len(stack) == 1:
                raise ValueError("unbalanced ')'")
            done = stack.pop()
            stack[-1].append(done)
        else:
            if tok[0] == '|':
                tok = tok[1:-1]
            stack[-1].append(tok)
    if len(stack) != 1:
        raise ValueError("unbalanced '('")
    return stack[0]


def numeral(tok):
    if isinstance(tok, str) and re.fullmatch(r'[0-9]+(\.[0-9]+)?', tok):
        return Fraction(tok)
    return None


# ---------------------------------------------------------------------------
# Formula DAG.  Real nodes: Lin (affine), RIte, RAdd, RScale.  Bool nodes below.
# ---------------------------------------------------------------------------
class Lin:
    __slots__ = ("a", "c")

    def __init__(self, a, c):
        self.a, self.c = a, float(c)

    def is_const(self):
        return not np.any(self.a)


class RIte:
    __slots__ = ("cond", "t", "e")

    def __init__(self, cond, t, e):
        self.cond, self.t, self.e = cond, t, e


class RAdd:
    __slots__ = ("terms",)

    def __init__(self, terms):
        self.terms = terms


class RScale:
    __slots__ = ("k", "x")

    def __init__(self, k, x):
        self.k, self.x = float(k), x


REAL = (Lin, RIte, RAdd, RScale)


class B:
    """Bool node: kind in const/var/atom/not/and/or/xor/iff/ite."""
    __slots__ = ("kind", "args", "val")

    def __init__(self, kind, args=(), val=None):
        self.kind, self.args, self.val = kind, tuple(args), val


def b_and(xs):
    return xs[0] if len(xs) == 1 else B("and", xs)


def b_or(xs):
    return xs[0] if len(xs) == 1 else B("or", xs)


def b_not(x):
    return B("not", [x])


class Formula:
    def __init__(self, text):
        self.reals, self.bools = [], []
        self.globals = {}
        asserts = []
        for cmd in parse_sexprs(text):
            if not isinstance(cmd, list) or not cmd:
                continue
            head = cmd[0]
            if head in ("declare-const", "declare-fun"):
                name = cmd[1]
                sort = cmd[2] if head == "declare-const" else cmd[3]
                if head == "declare-fun" and cmd[2]:
                    raise ValueError("uninterpreted functions are not supported")
                if sort == "Real":
                    self.reals.append(name)
                elif sort == "Bool":
                    self.bools.append(name)
                else:
                    raise ValueError("unsupported sort %r" % (sort,))
            elif head == "define-fun":
                if cmd[2]:
                    raise ValueError("define-fun with arguments is not supported")
                asserts.append(("define", cmd[1], cmd[4]))
            elif head == "assert":
                asserts.append(("assert", cmd[1]))
        self.n = len(self.reals)
        self.index = {v: i for i, v in enumerate(self.reals)}
        self.ttc_order = self._first_occurrence([a[1] for a in asserts if a[0] == "assert"])
        roots = []
        for item in asserts:
            if item[0] == "define":
                self.globals[item[1]] = self.compile(item[2], self.globals)
            else:
                roots.append(self.compile(item[1], self.globals))
        self.root = b_and(roots) if roots else B("const", val=True)
        self._collect_atoms()

    def _first_occurrence(self, exprs):
        """Real variables in the order ttc numbers its columns: first occurrence
        in a DFS (children left to right, let-bound names expanded where used)
        of the raw asserted term (src/parser.cpp, `collect`)."""
        order, seen, done = [], set(), set()

        def visit(e, env):
            key = (id(e), id(env))
            if key in done:
                return
            done.add(key)
            if isinstance(e, str):
                if e in env:
                    visit(*env[e])
                elif e in self.index and e not in seen:
                    seen.add(e)
                    order.append(e)
                return
            if not e:
                return
            if e[0] == "let":
                new = dict(env)
                for name, term in e[1]:
                    new[name] = (term, env)
                visit(e[2], new)
                return
            if e[0] == "!":
                visit(e[1], env)
                return
            for x in e[1:]:
                visit(x, env)

        for e in exprs:
            visit(e, {})
        return order

    # ---- compilation ------------------------------------------------------
    def lin_const(self, v):
        return Lin(np.zeros(self.n), v)

    def compile(self, e, env):
        if isinstance(e, str):
            if e in env:
                return env[e]
            if e == "true" or e == "false":
                return B("const", val=(e == "true"))
            if e in self.index:
                a = np.zeros(self.n)
                a[self.index[e]] = 1.0
                return Lin(a, 0.0)
            if e in self.bools:
                return B("var", val=e)
            v = numeral(e)
            if v is not None:
                return self.lin_const(float(v))
            raise ValueError("unknown symbol %r" % e)
        head, args = e[0], e[1:]
        if head == "let":
            new = dict(env)
            for name, term in args[0]:
                new[name] = self.compile(term, env)
            return self.compile(args[1], new)
        if head == "!":
            return self.compile(args[0], env)
        if isinstance(head, list):
            raise ValueError("unsupported application")
        xs = [self.compile(a, env) for a in args]
        if head in ("+", "-", "*", "/", "to_real", "abs"):
            return self.arith(head, xs)
        if head == "not":
            return b_not(xs[0])
        if head == "and":
            return b_and(xs) if xs else B("const", val=True)
        if head == "or":
            return b_or(xs) if xs else B("const", val=False)
        if head == "=>":
            out = xs[-1]
            for x in reversed(xs[:-1]):
                out = B("or", [b_not(x), out])
            return out
        if head == "xor":
            return B("xor", xs)
        if head == "ite":
            c, t, f = xs
            if isinstance(t, REAL):
                if c.kind == "const":
                    return t if c.val else f
                return RIte(c, t, f)
            return B("ite", [c, t, f])
        if head in ("=", "distinct"):
            if isinstance(xs[0], REAL):
                if head == "=":
                    return b_and([B("atom", [self.sub(xs[i], xs[i + 1])], "=")
                                  for i in range(len(xs) - 1)])
                return b_and([b_not(B("atom", [self.sub(xs[i], xs[j])], "="))
                              for i in range(len(xs)) for j in range(i + 1, len(xs))])
            if head == "=":
                return B("iff", xs)
            if len(xs) == 2:
                return B("xor", xs)
            return B("const", val=False)
        if head in ("<", "<=", ">", ">="):
            return b_and([B("atom", [self.sub(xs[i], xs[i + 1])], head)
                          for i in range(len(xs) - 1)])
        raise ValueError("unsupported operator %r" % head)

    def add(self, xs):
        if all(isinstance(x, Lin) for x in xs):
            return Lin(sum((x.a for x in xs), np.zeros(self.n)), sum(x.c for x in xs))
        return RAdd(xs)

    def scale(self, k, x):
        if isinstance(x, Lin):
            return Lin(k * x.a, k * x.c)
        return RScale(k, x)

    def sub(self, x, y):
        return self.add([x, self.scale(-1.0, y)])

    def arith(self, head, xs):
        if head == "to_real":
            return xs[0]
        if head == "abs":
            z = self.lin_const(0.0)
            return RIte(B("atom", [self.sub(xs[0], z)], ">="), xs[0], self.scale(-1.0, xs[0]))
        if head == "+":
            return self.add(xs)
        if head == "-":
            if len(xs) == 1:
                return self.scale(-1.0, xs[0])
            return self.add([xs[0]] + [self.scale(-1.0, x) for x in xs[1:]])
        if head == "*":
            k, rest = 1.0, []
            for x in xs:
                if isinstance(x, Lin) and x.is_const():
                    k *= x.c
                else:
                    rest.append(x)
            if not rest:
                return self.lin_const(k)
            if len(rest) > 1:
                raise ValueError("non-linear multiplication")
            return self.scale(k, rest[0])
        if head == "/":
            out = xs[0]
            for d in xs[1:]:
                if not (isinstance(d, Lin) and d.is_const()) or d.c == 0.0:
                    raise ValueError("division by a non-constant")
                out = self.scale(1.0 / d.c, out)
            return out
        raise AssertionError(head)

    # ---- symbolic branches of a real node: list of (conds, a, c) -----------
    def branches(self, r):
        if isinstance(r, Lin):
            return [((), r.a, r.c)]
        if isinstance(r, RScale):
            return [(cs, r.k * a, r.k * c) for cs, a, c in self.branches(r.x)]
        if isinstance(r, RIte):
            return ([(((r.cond, True),) + cs, a, c) for cs, a, c in self.branches(r.t)]
                    + [(((r.cond, False),) + cs, a, c) for cs, a, c in self.branches(r.e)])
        out = [((), np.zeros(self.n), 0.0)]
        for t in r.terms:
            out = [(c1 + c2, a1 + a2, k1 + k2) for c1, a1, k1 in out
                   for c2, a2, k2 in self.branches(t)]
            if len(out) > DNF_CAP:
                raise ValueError("too many ite branches")
        return out

    def _collect_atoms(self):
        seen, stack, hyps = set(), [self.root], []
        while stack:
            node = stack.pop()
            if id(node) in seen:
                continue
            seen.add(id(node))
            if isinstance(node, B):
                if node.kind == "atom":
                    for cs, a, c in self.branches(node.args[0]):
                        if np.any(a):
                            hyps.append((a, c))
                    stack.append(node.args[0])
                else:
                    stack.extend(node.args)
            elif isinstance(node, RIte):
                stack.extend([node.cond, node.t, node.e])
            elif isinstance(node, RAdd):
                stack.extend(node.terms)
            elif isinstance(node, RScale):
                stack.append(node.x)
        # Hyperplanes a.x + c = 0, deduplicated after normalisation.
        uniq = {}
        for a, c in hyps:
            nrm = np.linalg.norm(a)
            a2, c2 = a / nrm, c / nrm
            key = tuple(np.round(np.append(a2, c2), 12))
            uniq.setdefault(key, (a2, c2))
        self.hyperplanes = list(uniq.values())

    # ---- vectorised evaluation ---------------------------------------------
    def evaluate(self, X):
        """F(x) for every row of X, Boolean variables existentially quantified."""
        k = len(self.bools)
        if k > 14:
            raise ValueError("too many free Boolean variables")
        res = np.zeros(len(X), dtype=bool)
        for bits in range(1 << k):
            assign = {v: bool((bits >> i) & 1) for i, v in enumerate(self.bools)}
            res |= self._eval_bool(self.root, X, {}, assign)
        return res

    def _eval_real(self, r, X, memo, assign):
        key = id(r)
        if key in memo:
            return memo[key]
        if isinstance(r, Lin):
            out = [(None, r.a, r.c)]
        elif isinstance(r, RScale):
            out = [(m, r.k * a, r.k * c) for m, a, c in self._eval_real(r.x, X, memo, assign)]
        elif isinstance(r, RIte):
            cm = self._eval_bool(r.cond, X, memo, assign)
            out = []
            for m, a, c in self._eval_real(r.t, X, memo, assign):
                out.append((cm if m is None else (m & cm), a, c))
            for m, a, c in self._eval_real(r.e, X, memo, assign):
                out.append((~cm if m is None else (m & ~cm), a, c))
        else:
            out = [(None, np.zeros(self.n), 0.0)]
            for t in r.terms:
                nxt = []
                for m1, a1, c1 in out:
                    for m2, a2, c2 in self._eval_real(t, X, memo, assign):
                        m = m1 if m2 is None else (m2 if m1 is None else m1 & m2)
                        nxt.append((m, a1 + a2, c1 + c2))
                out = nxt
        memo[key] = out
        return out

    def _eval_bool(self, node, X, memo, assign):
        key = id(node)
        if key in memo:
            return memo[key]
        k = node.kind
        if k == "const":
            out = np.full(len(X), node.val)
        elif k == "var":
            out = np.full(len(X), assign[node.val])
        elif k == "atom":
            out = np.zeros(len(X), dtype=bool)
            for m, a, c in self._eval_real(node.args[0], X, memo, assign):
                v = X @ a + c
                op = node.val
                r = (v < 0) if op == "<" else (v <= 0) if op == "<=" else \
                    (v > 0) if op == ">" else (v >= 0) if op == ">=" else (v == 0)
                out |= r if m is None else (r & m)
        else:
            vals = [self._eval_bool(a, X, memo, assign) for a in node.args]
            if k == "not":
                out = ~vals[0]
            elif k == "and":
                out = np.logical_and.reduce(vals)
            elif k == "or":
                out = np.logical_or.reduce(vals)
            elif k == "xor":
                out = np.logical_xor.reduce(vals)
            elif k == "iff":
                out = np.logical_and.reduce([vals[i] == vals[i + 1]
                                             for i in range(len(vals) - 1)])
            elif k == "ite":
                out = np.where(vals[0], vals[1], vals[2])
            else:
                raise AssertionError(k)
        memo[key] = out
        return out

    # ---- independent DNF of F (list of half-space systems) -------------------
    def dnf(self):
        """Disjuncts of F as lists of (a, c, strict): a.x + c (<|<=) 0.

        Negated equalities are dropped (true almost everywhere); Boolean
        variables are projected out.  Returns None when the DNF exceeds
        DNF_CAP terms."""
        memo = {}

        def consistent(c):
            vs = {(l[1], l[2]) for l in c if l[0] == "v"}
            return not any((v, not p) in vs for v, p in vs)

        def conj_and(lists):
            out = [frozenset()]
            for L in lists:
                out = [x | y for x in out for y in L]
                out = [c for c in out if consistent(c)]
                if len(out) > DNF_CAP:
                    raise OverflowError
            return list(set(out))

        def conj_or(lists):
            out = list({c for L in lists for c in L})
            if len(out) > DNF_CAP:
                raise OverflowError
            return out

        def lit(a, c, strict):
            return ("h", tuple(np.round(a, 15)), round(c, 15), strict)

        def d(node, pol):
            key = (id(node), pol)
            if key in memo:
                return memo[key]
            k = node.kind
            if k == "const":
                out = [frozenset()] if node.val == pol else []
            elif k == "var":
                out = [frozenset([("v", node.val, pol)])]
            elif k == "not":
                out = d(node.args[0], not pol)
            elif k in ("and", "or"):
                parts = [d(a, pol) for a in node.args]
                out = conj_and(parts) if (k == "and") == pol else conj_or(parts)
            elif k == "ite":
                c, t, e = node.args
                out = conj_or([conj_and([d(c, True), d(t, pol)]),
                               conj_and([d(c, False), d(e, pol)])])
            elif k == "xor" or k == "iff":
                args = list(node.args)
                if k == "iff":
                    if len(args) > 2:
                        out = conj_and([d(B("iff", [args[i], args[i + 1]]), pol)
                                        for i in range(len(args) - 1)]) if pol else \
                            conj_or([d(B("iff", [args[i], args[i + 1]]), False)
                                     for i in range(len(args) - 1)])
                        memo[key] = out
                        return out
                    pol2 = not pol
                else:
                    pol2 = pol
                # parity of args == pol2
                acc = {True: d(args[0], True), False: d(args[0], False)}
                for a in args[1:]:
                    at, af = d(a, True), d(a, False)
                    acc = {p: conj_or([conj_and([acc[p], af]), conj_and([acc[not p], at])])
                           for p in (True, False)}
                out = acc[pol2]
            elif k == "atom":
                op = node.val
                out = []
                for conds, a, c in self.branches(node.args[0]):
                    cl = [d(cn, cp) for cn, cp in conds]
                    if op == "=":
                        if pol:
                            hl = [frozenset([lit(a, c, False), lit(-a, -c, False)])]
                        else:
                            hl = [frozenset()]
                    else:
                        if not pol:
                            op2 = {"<": ">=", "<=": ">", ">": "<=", ">=": "<"}[op]
                        else:
                            op2 = op
                        if op2 in ("<", "<="):
                            hl = [frozenset([lit(a, c, op2 == "<")])]
                        else:
                            hl = [frozenset([lit(-a, -c, op2 == ">")])]
                    out = conj_or([out, conj_and(cl + [hl])])
            else:
                raise AssertionError(k)
            memo[key] = out
            return out

        try:
            terms = d(self.root, True)
        except (OverflowError, ValueError):
            return None
        systems = []
        for t in terms:
            hs = [(np.array(l[1]), l[2], l[3]) for l in t if l[0] == "h"]
            systems.append(hs)
        return systems


# ---------------------------------------------------------------------------
# Geometry helpers
# ---------------------------------------------------------------------------
def read_ine(path, n, perm=None):
    """A x <= b in our column order; perm[j] = our index of ttc's column j."""
    with open(path) as f:
        lines = [l.split() for l in f if l.strip()]
    i = next(k for k, l in enumerate(lines) if l[0] == "begin")
    rows, cols = int(lines[i + 1][0]), int(lines[i + 1][1])
    if cols != n + 1:
        raise ValueError("%s: %d columns, expected %d" % (path, cols, n + 1))
    M = np.array([[float(Fraction(v)) for v in l] for l in lines[i + 2:i + 2 + rows]])
    A = np.zeros((rows, n))
    A[:, perm if perm is not None else np.arange(n)] = -M[:, 1:]
    return A, M[:, 0]                  # A x <= b


class Frame:
    """Affine frame x = mid + s*y.  LPs, qhull and hit-and-run work in y with unit
    row norms, so their absolute tolerances do not depend on the formula's units."""

    def __init__(self, mid, s):
        self.mid, self.s = np.asarray(mid, float), float(s)

    def unit(self, A, b):
        A2 = A * self.s
        b2 = b - A @ self.mid
        nr = np.linalg.norm(A2, axis=1)
        nr[nr == 0] = 1.0
        return A2 / nr[:, None], b2 / nr

    def to_x(self, Y):
        return self.mid + self.s * Y

    def to_y(self, X):
        return (X - self.mid) / self.s


def chebyshev(A, b, fr):
    """Chebyshev centre and radius of {A x <= b} (scipy LP), or (None, -1)."""
    if not HAVE_SCIPY or len(A) == 0:
        return None, -1.0
    A, b = fr.unit(A, b)
    n = A.shape[1]
    cost = np.zeros(n + 1)
    cost[-1] = -1.0
    res = scipy.optimize.linprog(cost, A_ub=np.hstack([A, np.ones((len(A), 1))]), b_ub=b,
                                 bounds=[(None, None)] * n + [(0, None)], method="highs")
    if res.status == 3:
        return None, math.inf
    if res.status != 0:
        return None, -1.0
    return fr.to_x(res.x[:n]), res.x[n] * fr.s


def lp_bbox(A, b, fr):
    A, b = fr.unit(A, b)
    n = A.shape[1]
    lo, hi = np.empty(n), np.empty(n)
    for j in range(n):
        e = np.zeros(n)
        e[j] = 1.0
        r1 = scipy.optimize.linprog(e, A_ub=A, b_ub=b, bounds=[(None, None)] * n, method="highs")
        r2 = scipy.optimize.linprog(-e, A_ub=A, b_ub=b, bounds=[(None, None)] * n, method="highs")
        if r1.status != 0 or r2.status != 0:
            return None
        lo[j], hi[j] = r1.fun, -r2.fun
    return fr.to_x(lo), fr.to_x(hi)


def hit_and_run(A, b, x0, count, rng, fr, chains=64):
    """~uniform points of {A x <= b} from an interior x0 (vectorised chains)."""
    A, b = fr.unit(A, b)
    n = A.shape[1]
    per = max(1, -(-count // chains))
    X = np.tile(fr.to_y(x0), (chains, 1)).astype(float)
    out = []
    burn, thin = 20 * n + 50, n + 2
    for step in range(burn + per * thin):
        D = rng.standard_normal((chains, n))
        D /= np.linalg.norm(D, axis=1, keepdims=True)
        AD = D @ A.T
        slack = np.maximum(b[None, :] - X @ A.T, 0.0)
        with np.errstate(divide="ignore", invalid="ignore"):
            t = slack / AD
        tmax = np.where(AD > 1e-300, t, np.inf).min(axis=1)
        tmin = np.where(AD < -1e-300, t, -np.inf).max(axis=1)
        ok = np.isfinite(tmax) & np.isfinite(tmin)
        tmax, tmin = np.where(ok, tmax, 0.0), np.where(ok, tmin, 0.0)
        X = X + (tmin + (tmax - tmin) * rng.random(chains))[:, None] * D
        if step >= burn and (step - burn) % thin == thin - 1:
            out.append(X.copy())
    return fr.to_x(np.vstack(out)[:count])


def polytope_volume(A, b, center, fr):
    if not HAVE_SCIPY or A.shape[1] > 8:
        return None
    A2, b2 = fr.unit(A, b)
    n = A.shape[1]
    try:
        hs = scipy.spatial.HalfspaceIntersection(np.hstack([A2, -b2[:, None]]), fr.to_y(center))
        if n == 1:
            v = float(np.ptp(hs.intersections))
        else:
            v = float(scipy.spatial.ConvexHull(hs.intersections).volume)
        return v * fr.s ** n
    except Exception:
        return None


def near_plane_points(base, a, c, widths, rng, count):
    """Points at distance w*10^U(-6,-2) on both sides of {a.x + c = 0}, where w is
    the width of the bounding box along the unit normal."""
    if len(base) == 0 or count <= 0:
        return np.empty((0, len(a)))
    nrm = np.linalg.norm(a)
    a, c = a / nrm, c / nrm
    w = float(np.abs(a) @ widths)
    P = base[rng.integers(0, len(base), size=count)]
    P = P - (P @ a + c)[:, None] * a[None, :]
    s = w * 10.0 ** rng.uniform(-6, -2, size=count)
    s *= np.where(np.arange(count) % 2 == 0, 1.0, -1.0)
    return P + s[:, None] * a[None, :]


# ---------------------------------------------------------------------------
# ttc driver
# ---------------------------------------------------------------------------
def run_ttc(ttc, path, prefix, full, timeout):
    cmd = ["nice", "-n", "19", ttc, path, "--dump-ine", prefix, "-v", "2"]
    info = {"cubes": None, "skipped": 0, "volume_mode": False, "error": None,
            "returncode": None, "s_volume": None, "log": []}
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                            text=True, errors="replace")
    timer = threading.Timer(timeout, proc.kill)
    timer.start()
    stopped = False
    try:
        for line in proc.stdout:
            info["log"].append(line.rstrip("\n"))
            m = re.match(r"c polytopes: (\d+)", line)
            if m:
                info["cubes"] = int(m.group(1))
            m = re.match(r"Skipping (\d+) polytopes out of (\d+)", line)
            if m:
                info["skipped"] = int(m.group(1))
            if line.startswith("c counting: volume"):
                info["volume_mode"] = True
            m = re.match(r"s volume (\S+)", line)
            if m:
                info["s_volume"] = float(m.group(1))
            if re.search(r"[Ee]rror|exception", line) and info["error"] is None:
                info["error"] = line.strip()
            if not full and "[ volume computation ]" in line:
                stopped = True       # phase 1a (the dump) is complete
                proc.terminate()
                break
    finally:
        proc.wait()
        timer.cancel()
    info["returncode"] = 0 if stopped else proc.returncode
    return info


def cvc5_check(cvc5, text, reals, x):
    """Exact F(x) via cvc5 (None if unavailable)."""
    if not cvc5:
        return None
    body = re.sub(r"\(\s*(check-sat|exit|get-model)\s*\)", "", text)
    eqs = []
    for v, val in zip(reals, x):
        q = Fraction(float(val))
        num = "(- %d)" % -q.numerator if q.numerator < 0 else "%d" % q.numerator
        eqs.append("(= |%s| (/ %s %d))" % (v, num, q.denominator))
    query = body + "\n(assert (and %s))\n(check-sat)\n" % " ".join(eqs)
    try:
        r = subprocess.run(["nice", "-n", "19", cvc5, "--lang=smt2"], input=query,
                           capture_output=True, text=True, timeout=60)
    except Exception:
        return None
    out = r.stdout.strip().splitlines()
    return {"sat": True, "unsat": False}.get(out[-1] if out else "", None)


# ---------------------------------------------------------------------------
# Checking one formula
# ---------------------------------------------------------------------------
def headers(text):
    meta = {}
    for line in text.splitlines():
        m = re.match(r"\s*;\s*([a-z-]+):\s*(.*)$", line)
        if m:
            meta[m.group(1)] = m.group(2).strip()
    return meta


def check(path, args, workdir):
    name = os.path.basename(path)
    for suf in (".xz", ".smt2"):
        if name.endswith(suf):
            name = name[: -len(suf)]
    if path.endswith(".xz"):
        text = lzma.open(path, "rt").read()
        path = os.path.join(workdir, name + ".smt2")
        with open(path, "w") as f:
            f.write(text)
    else:
        text = open(path).read()
    meta = headers(text)
    rep = {"name": name, "report_only": meta.get("report-only", "") == "true",
           "fails": [], "notes": []}
    rng = np.random.default_rng(args.seed)
    F = Formula(text)
    n = F.n
    rep["dim"] = n
    prefix = os.path.join(workdir, name)
    for old in glob.glob(prefix + "_cube*.ine"):
        os.remove(old)
    t0 = time.time()
    info = run_ttc(args.ttc, path, prefix, args.full, args.timeout)
    rep["ttc_time"] = time.time() - t0
    rep["cubes"], rep["skipped"] = info["cubes"], info["skipped"]
    if not info["volume_mode"]:
        rep["fails"].append("ttc did not run the LRA volume engine")
    if info["returncode"] != 0 or (info["error"] and info["volume_mode"]):
        rep["fails"].append("ttc error (rc=%s): %s" % (info["returncode"], info["error"]))
    files = sorted(glob.glob(prefix + "_cube*.ine"),
                   key=lambda p: int(re.search(r"_cube(\d+)\.ine$", p).group(1)))
    if len(F.ttc_order) != n:
        rep["fails"].append("declared but unused Real variables: %s (ttc drops them)"
                            % sorted(set(F.reals) - set(F.ttc_order)))
        return rep
    perm = np.array([F.index[v] for v in F.ttc_order])
    rep["ttc_column_order"] = F.ttc_order
    polys = [read_ine(p, n, perm) for p in files]
    rep["dumped"] = len(polys)

    # Independent DNF of F (LPs in the header-bbox frame when there is one).
    hdr_box = None
    if "bbox" in meta:
        v = [float(Fraction(t)) for t in meta["bbox"].split()]
        hdr_box = (np.full(n, v[0]), np.full(n, v[1])) if len(v) == 2 else \
            (np.array(v[0::2]), np.array(v[1::2]))
        fr0 = Frame((hdr_box[0] + hdr_box[1]) / 2, np.max(hdr_box[1] - hdr_box[0]) / 2)
    else:
        fr0 = Frame(np.zeros(n), 1.0)
    dnf = F.dnf()
    regions = []
    if dnf is not None and HAVE_SCIPY:
        for hs in dnf:
            if not hs:
                continue
            A = np.array([h[0] for h in hs])
            b = -np.array([h[1] for h in hs])
            ctr, r = chebyshev(A, b, fr0)
            if ctr is not None and r > 1e-12 * fr0.s:
                regions.append((A, b, ctr, r))
    rep["dnf_terms"] = None if dnf is None else len(dnf)
    rep["dnf_regions"] = len(regions)

    # Bounding box.
    if hdr_box is not None:
        lo, hi = hdr_box
    else:
        boxes = [lp_bbox(A, b, fr0) for A, b, _, _ in regions]
        if not regions:
            for A, b in polys:
                ctr, r = chebyshev(A, b, fr0)
                if ctr is not None and r > 0:
                    boxes.append(lp_bbox(A, b, fr0))
        boxes = [bx for bx in boxes if bx is not None]
        if not boxes:
            rep["fails"].append("no bounding box")
            return rep
        lo = np.min([bx[0] for bx in boxes], axis=0)
        hi = np.max([bx[1] for bx in boxes], axis=0)
        w = hi - lo
        lo, hi = lo - 0.1 * w, hi + 0.1 * w
    fr = Frame((lo + hi) / 2, np.max(hi - lo) / 2)

    # Per-polytope interior point, radius, exact volume.
    pinfo = []
    for A, b in polys:
        ctr, r = chebyshev(A, b, fr)
        vol = polytope_volume(A, b, ctr, fr) if (ctr is not None and 0 < r < math.inf) else None
        pinfo.append({"center": ctr, "radius": r, "vol": vol})
    scale = float(np.max(hi - lo))
    # Floating-point guard only: any distance > 0 already gives a ball of mismatch.
    tau = 1e-9 * float(np.max(np.maximum(np.abs(lo), np.abs(hi))))
    box_vol = float(np.prod(hi - lo))

    # ---- point sets --------------------------------------------------------
    sets = {}
    sets["uniform"] = lo + (hi - lo) * rng.random((args.points, n))
    pin = []
    for (A, b), pi in zip(polys, pinfo):
        if pi["center"] is not None and 0 < pi["radius"] < math.inf:
            pin.append(hit_and_run(A, b, pi["center"], args.inner, rng, fr))
        else:
            pin.append(np.empty((0, n)))
    sets["in_P"] = np.vstack(pin) if pin else np.empty((0, n))
    rin = [hit_and_run(A, b, ctr, args.inner, rng, fr) for A, b, ctr, r in regions
           if r < math.inf]
    sets["in_dnf"] = np.vstack(rin) if rin else np.empty((0, n))
    near = []
    for (A, b), S in zip(polys, pin):
        base = S if len(S) else sets["uniform"]
        for k in range(len(A)):
            near.append(near_plane_points(base, A[k], -b[k], hi - lo, rng, args.near))
    sets["near_P_facets"] = np.vstack(near) if near else np.empty((0, n))
    base = np.vstack([sets["uniform"], sets["in_P"], sets["in_dnf"]])
    near = [near_plane_points(base, a, c, hi - lo, rng, args.near) for a, c in F.hyperplanes]
    sets["near_F_atoms"] = np.vstack(near) if near else np.empty((0, n))

    # Hyperplanes for the boundary classification (unit normals).
    H, h0 = [], []
    for a, c in F.hyperplanes:
        H.append(a)
        h0.append(c)
    for A, b in polys:
        nr = np.linalg.norm(A, axis=1)
        H.extend(A / nr[:, None])
        h0.extend(-b / nr)
    H, h0 = (np.array(H), np.array(h0)) if H else (np.zeros((0, n)), np.zeros(0))

    def classify(X):
        f = np.zeros(len(X), dtype=bool)
        cnt = np.zeros(len(X), dtype=int)
        inside = np.zeros((len(X), len(polys)), dtype=bool)
        dist = np.full(len(X), np.inf)
        for s in range(0, len(X), 4096):
            Y = X[s:s + 4096]
            f[s:s + 4096] = F.evaluate(Y)
            for i, (A, b) in enumerate(polys):
                inside[s:s + 4096, i] = np.all(Y @ A.T <= b, axis=1)
            if len(H):
                dist[s:s + 4096] = np.abs(Y @ H.T + h0).min(axis=1)
        cnt = inside.sum(axis=1)
        return f, cnt, inside, dist

    rep["sets"] = {}
    examples = []
    tot_int = tot_bnd = 0
    for key, X in sets.items():
        if len(X) == 0:
            rep["sets"][key] = {"n": 0}
            continue
        f, cnt, inside, dist = classify(X)
        mis = f != (cnt > 0)
        interior = mis & (dist > tau)
        inc = interior & f            # F true, no polytope covers x
        uns = interior & ~f           # some P_i contains x, F false
        d = {"n": int(len(X)), "F_true": int(f.sum()), "incomplete": int(inc.sum()),
             "unsound": int(uns.sum()), "boundary_mismatch": int((mis & ~interior).sum())}
        if key == "in_dnf":
            d["F_false_in_dnf"] = int(((~f) & (dist > tau)).sum())
            if d["F_false_in_dnf"]:
                rep["fails"].append("evaluator/DNF disagree on %d points" % d["F_false_in_dnf"])
        rep["sets"][key] = d
        tot_int += d["incomplete"] + d["unsound"]
        tot_bnd += d["boundary_mismatch"]
        for idx in np.flatnonzero(interior)[: 2]:
            examples.append({"set": key, "x": X[idx].tolist(), "F": bool(f[idx]),
                             "in_polytopes": [int(i) + 1 for i in np.flatnonzero(inside[idx])],
                             "dist": float(dist[idx])})
        if key == "uniform":
            rep["F_true_uniform"] = int(f.sum())
            p = f.mean()
            rep["volF_mc"] = box_vol * p
            rep["volF_se"] = box_vol * math.sqrt(max(p * (1 - p), 1e-300) / len(X))
            q = (cnt > 0).mean()
            rep["volCover_mc"] = box_vol * q
            # 95% one-sided upper bound on the mismatch-measure fraction.
            k = int(mis.sum())
            rep["mismatch_frac_ub95"] = (3.0 / len(X)) if k == 0 else \
                (k + 1.645 * math.sqrt(k) + 1) / len(X)
    rep["interior_mismatches"] = tot_int
    rep["boundary_mismatches"] = tot_bnd
    rep["examples"] = examples
    if tot_int:
        rep["fails"].append("%d interior decomposition mismatches" % tot_int)
        for ex in examples[: args.confirm]:
            ex["cvc5_F"] = cvc5_check(args.cvc5, text, F.reals, ex["x"])
            if ex["cvc5_F"] is not None and ex["cvc5_F"] != ex["F"]:
                rep["fails"].append("cvc5 disagrees with the evaluator at %s" % ex["x"])

    # ---- overlap ratio ---------------------------------------------------
    vols = [pi["vol"] for pi in pinfo if pi["radius"] > 0]
    rep["degenerate_dumped"] = sum(1 for pi in pinfo if pi["radius"] <= 1e-12 * scale)
    if vols and all(v is not None for v in vols):
        sv = sum(vols)
        rep["sum_vol_P"] = sv
        if rep.get("volF_mc"):
            rep["overlap_box"] = sv / rep["volF_mc"]
        # Karp-Luby: vol(U P) = sum_i v_i E_{x~P_i}[1/c(x)].
        kl, used = 0.0, 0
        for i, (S, pi) in enumerate(zip(pin, pinfo)):
            if pi["radius"] <= 0 or len(S) == 0:
                continue
            _, cnt, _, _ = classify(S)
            kl += pi["vol"] * np.mean(1.0 / np.maximum(cnt, 1))
            used += 1
        if kl > 0:
            rep["overlap_kl"] = sv / kl
            rep["volCover_kl"] = kl
    if "exact-volume" in meta:
        ex = float(meta["exact-volume"])
        rep["exact"] = ex
        # Evaluator sanity: k = #F-true uniform points ~ Binomial(N, p0), p0 = exact/box.
        N, k = args.points, rep.get("F_true_uniform", 0)
        p0 = min(ex / box_vol, 1.0)
        mu = N * p0
        if mu >= 20 and p0 < 1:
            z = (k - mu) / math.sqrt(mu * (1 - p0))
            rep["z_exact"] = z
            if abs(z) > 4:
                rep["fails"].append("evaluator sanity: MC vol(F) %.4g vs exact %.4g (z=%.1f)"
                                    % (rep["volF_mc"], ex, z))
        elif k == 0 and mu > 0 and math.exp(-mu) < 1e-4:
            rep["fails"].append("evaluator sanity: no F-true samples, expected %.1f" % mu)
    rep["s_volume"] = info["s_volume"]
    if args.verbose:
        rep["log"] = info["log"][-40:]
    return rep


def fmt(v, f="%.4g"):
    return "-" if v is None else (f % v)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("inputs", nargs="*", help="extra .smt2 / .smt2.xz files")
    ap.add_argument("--ttc", default=os.path.join(REPO, "build", "ttc"))
    ap.add_argument("--cvc5", default=shutil.which("cvc5"))
    ap.add_argument("--filter", default=None)
    ap.add_argument("--no-default", action="store_true", help="skip tests/volesti/union_decomp")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--points", type=int, default=20000, help="uniform box points")
    ap.add_argument("--inner", type=int, default=512, help="points per polytope / DNF region")
    ap.add_argument("--near", type=int, default=64, help="points per hyperplane (both sides)")
    ap.add_argument("--confirm", type=int, default=3, help="cvc5 checks per failing input")
    ap.add_argument("--full", action="store_true", help="let ttc finish (prints s volume)")
    ap.add_argument("--timeout", type=float, default=300.0)
    ap.add_argument("--workdir", default=None)
    ap.add_argument("--keep", action="store_true")
    ap.add_argument("--json", default=None)
    ap.add_argument("-v", "--verbose", action="store_true")
    args = ap.parse_args()
    try:
        os.nice(19)
    except OSError:
        pass
    files = [] if args.no_default else sorted(glob.glob(os.path.join(INPUT_DIR, "*.smt2")))
    files += args.inputs
    if args.filter:
        files = [f for f in files if re.search(args.filter, os.path.basename(f))]
    workdir = args.workdir or tempfile.mkdtemp(prefix="union_decomp_")
    os.makedirs(workdir, exist_ok=True)
    if not HAVE_SCIPY:
        print("note: scipy missing -> no DNF-region sampling, no exact vol(P_i)")
    hdr = ("%-44s %3s %5s %4s %4s %4s %7s %6s %6s %6s %9s %9s %9s %7s %7s  %s" %
           ("formula", "n", "cubes", "dump", "skip", "dnf", "points", "incmp", "unsnd",
            "bndry", "vol(F)MC", "exact", "sumVolP", "ovl_box", "ovl_KL", "status"))
    print(hdr)
    print("-" * len(hdr))
    reports, bad = [], 0
    for path in files:
        try:
            rep = check(path, args, workdir)
        except Exception as e:  # evaluator limitation or crash
            rep = {"name": os.path.basename(path), "fails": ["harness error: %r" % (e,)],
                   "report_only": False}
        reports.append(rep)
        status = "PASS" if not rep["fails"] else ("REPORT" if rep["report_only"] else "FAIL")
        if status == "FAIL":
            bad += 1
        s = rep.get("sets", {})
        npts = sum(v.get("n", 0) for v in s.values())
        inc = sum(v.get("incomplete", 0) for v in s.values())
        uns = sum(v.get("unsound", 0) for v in s.values())
        print("%-44s %3s %5s %4s %4s %4s %7d %6d %6d %6s %9s %9s %9s %7s %7s  %s" % (
            rep["name"][:44], rep.get("dim", "-"), fmt(rep.get("cubes"), "%d"),
            fmt(rep.get("dumped"), "%d"), fmt(rep.get("skipped"), "%d"),
            fmt(rep.get("dnf_regions"), "%d"), npts, inc, uns,
            fmt(rep.get("boundary_mismatches"), "%d"), fmt(rep.get("volF_mc")),
            fmt(rep.get("exact")), fmt(rep.get("sum_vol_P")), fmt(rep.get("overlap_box"), "%.3f"),
            fmt(rep.get("overlap_kl"), "%.3f"), status))
        for msg in rep["fails"]:
            print("    - " + msg)
        for ex in rep.get("examples", [])[:3]:
            print("      e.g. [%s] x=%s F=%s in P%s dist=%.2g%s" % (
                ex["set"], np.array2string(np.array(ex["x"]), precision=6, separator=","),
                ex["F"], ex["in_polytopes"], ex["dist"],
                "" if ex.get("cvc5_F") is None else " cvc5:F=%s" % ex["cvc5_F"]))
        sys.stdout.flush()
    if args.json:
        with open(args.json, "w") as f:
            json.dump(reports, f, indent=1, default=float)
    if not args.keep and not args.workdir:
        shutil.rmtree(workdir, ignore_errors=True)
    print("\n%d input(s), %d FAIL" % (len(files), bad))
    return 1 if bad else 0


if __name__ == "__main__":
    sys.setrecursionlimit(200000)
    threading.stack_size(512 * 1024 * 1024)
    result = []
    t = threading.Thread(target=lambda: result.append(main()))
    t.start()
    t.join()
    sys.exit(result[0] if result else 2)
