#!/usr/bin/env python3
"""Random QF_LRA formulas with EXACTLY known volume, for union fuzzing.

A formula is a random Boolean combination (and/or/not/xor/=>/ite) of
axis-aligned atoms (y_i < c, y_i <= c, y_i > c, y_i >= c) inside a bounding
box.  Every atom's breakpoint is a coordinate of a compressed grid, so the
formula is constant on each grid cell and its exact volume is the sum of the
volumes of the cells whose centre satisfies it (fractions, exact).

Optionally the formula is written in sheared coordinates x = M y with a random
integer matrix M (|det M| >= 1): atoms a.y < c become (a M^-1).x < c with
rational coefficients, and vol_x = |det M| vol_y.  This gives non-axis-aligned
unions (skinny, rotated-ish) with exact volume.

Usage: gen_bool_unions.py --out DIR --count N --seed S [--dims 2,3,4] ...
Writes DIR/<name>.smt2 with '; exact-volume:' and '; dim:' headers.
"""
import argparse
import itertools
import os
import random
from fractions import Fraction


def frac_str(q):
    q = Fraction(q)
    if q.denominator == 1:
        return str(q.numerator) if q >= 0 else "(- %d)" % -q.numerator
    s = "(/ %d %d)" % (abs(q.numerator), q.denominator)
    return s if q >= 0 else "(- %s)" % s


def rand_breaks(rng, lo, hi, k):
    """k distinct rational breakpoints strictly inside (lo, hi)."""
    pts = set()
    while len(pts) < k:
        pts.add(lo + (hi - lo) * Fraction(rng.randint(1, 999), 1000))
    return sorted(pts)


class Gen:
    def __init__(self, rng, n, scale, atoms_per_dim, depth):
        self.rng, self.n = rng, n
        self.lo = [Fraction(0)] * n
        self.hi = [scale[i] for i in range(n)]
        self.breaks = [rand_breaks(rng, self.lo[i], self.hi[i], atoms_per_dim) for i in range(n)]
        self.depth = depth

    def atom(self):
        i = self.rng.randrange(self.n)
        c = self.rng.choice(self.breaks[i])
        op = self.rng.choice(["<", "<=", ">", ">="])
        return ("atom", i, op, c)

    def tree(self, d):
        if d == 0 or self.rng.random() < 0.25:
            return self.atom()
        op = self.rng.choices(["and", "or", "not", "xor", "=>", "ite"], [4, 4, 1, 1, 1, 1])[0]
        if op == "not":
            return ("not", self.tree(d - 1))
        if op == "ite":
            return ("ite", self.tree(d - 1), self.tree(d - 1), self.tree(d - 1))
        k = 2 if op in ("xor", "=>") else self.rng.randint(2, 3)
        return (op,) + tuple(self.tree(d - 1) for _ in range(k))


def evaluate(t, y):
    tag = t[0]
    if tag == "atom":
        _, i, op, c = t
        v = y[i]
        return {"<": v < c, "<=": v <= c, ">": v > c, ">=": v >= c}[op]
    if tag == "not":
        return not evaluate(t[1], y)
    if tag == "and":
        return all(evaluate(s, y) for s in t[1:])
    if tag == "or":
        return any(evaluate(s, y) for s in t[1:])
    if tag == "xor":
        return evaluate(t[1], y) != evaluate(t[2], y)
    if tag == "=>":
        return (not evaluate(t[1], y)) or evaluate(t[2], y)
    if tag == "ite":
        return evaluate(t[2], y) if evaluate(t[1], y) else evaluate(t[3], y)
    raise ValueError(tag)


def exact_volume(g, t):
    """Sum of grid-cell volumes whose centre satisfies t (in y coordinates)."""
    axes = []
    for i in range(g.n):
        pts = [g.lo[i]] + g.breaks[i] + [g.hi[i]]
        axes.append([(pts[k], pts[k + 1]) for k in range(len(pts) - 1)])
    vol = Fraction(0)
    for cell in itertools.product(*axes):
        centre = [(a + b) / 2 for a, b in cell]
        if evaluate(t, centre):
            v = Fraction(1)
            for a, b in cell:
                v *= b - a
            vol += v
    return vol


def mat_inv(M):
    n = len(M)
    A = [[Fraction(M[i][j]) for j in range(n)] + [Fraction(int(i == j)) for j in range(n)]
         for i in range(n)]
    for c in range(n):
        p = next(r for r in range(c, n) if A[r][c] != 0)
        A[c], A[p] = A[p], A[c]
        pv = A[c][c]
        A[c] = [v / pv for v in A[c]]
        for r in range(n):
            if r != c and A[r][c] != 0:
                f = A[r][c]
                A[r] = [a - f * b for a, b in zip(A[r], A[c])]
    return [row[n:] for row in A]


def det(M):
    n = len(M)
    A = [[Fraction(v) for v in row] for row in M]
    d = Fraction(1)
    for c in range(n):
        p = next((r for r in range(c, n) if A[r][c] != 0), None)
        if p is None:
            return Fraction(0)
        if p != c:
            A[c], A[p] = A[p], A[c]
            d = -d
        d *= A[c][c]
        for r in range(c + 1, n):
            f = A[r][c] / A[c][c]
            A[r] = [a - f * b for a, b in zip(A[r], A[c])]
    return d


def linear(coeffs, names):
    terms = []
    for a, nm in zip(coeffs, names):
        if a == 0:
            continue
        terms.append(nm if a == 1 else "(* %s %s)" % (frac_str(a), nm))
    if not terms:
        return "0"
    return terms[0] if len(terms) == 1 else "(+ %s)" % " ".join(terms)


def to_smt(t, row_of, names):
    """row_of(i) = coefficient vector of y_i as a linear form in x."""
    tag = t[0]
    if tag == "atom":
        _, i, op, c = t
        return "(%s %s %s)" % (op, linear(row_of(i), names), frac_str(c))
    return "(%s %s)" % (tag, " ".join(to_smt(s, row_of, names) for s in t[1:]))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--count", type=int, default=25)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--dims", default="2,3,4")
    ap.add_argument("--atoms-per-dim", type=int, default=4)
    ap.add_argument("--depth", type=int, default=4)
    ap.add_argument("--shear", type=float, default=0.5, help="fraction of sheared formulas")
    ap.add_argument("--max-cells", type=int, default=200000)
    args = ap.parse_args()
    rng = random.Random(args.seed)
    dims = [int(d) for d in args.dims.split(",")]
    os.makedirs(args.out, exist_ok=True)
    made = 0
    tries = 0
    while made < args.count and tries < 50 * args.count:
        tries += 1
        n = rng.choice(dims)
        k = max(1, min(args.atoms_per_dim, int(round(args.max_cells ** (1.0 / n))) - 1))
        scale = [Fraction(10) ** rng.randint(-3, 3) * rng.randint(1, 9) for _ in range(n)]
        g = Gen(rng, n, scale, k, args.depth)
        t = g.tree(args.depth)
        vol_y = exact_volume(g, t)
        box_vol = Fraction(1)
        for i in range(n):
            box_vol *= g.hi[i] - g.lo[i]
        if vol_y == 0 or vol_y == box_vol and rng.random() < 0.8:
            continue  # empty, or (mostly skip) the whole box
        names = ["x%d" % i for i in range(n)]
        sheared = rng.random() < args.shear
        if sheared:
            while True:
                M = [[rng.randint(-2, 2) if i != j else rng.choice([1, 1, 2, 3]) for j in range(n)]
                     for i in range(n)]
                dM = det(M)
                if dM != 0:
                    break
            Minv = mat_inv(M)          # y = Minv x
            row_of = lambda i: Minv[i]
            vol = vol_y * abs(dM)
        else:
            row_of = lambda i: [Fraction(int(i == j)) for j in range(n)]
            vol = vol_y
        box = " ".join("(>= %s %s) (<= %s %s)" % (linear(row_of(i), names), frac_str(g.lo[i]),
                                                   linear(row_of(i), names), frac_str(g.hi[i]))
                       for i in range(n))
        name = "bu_s%d_%03d_n%d%s" % (args.seed, made, n, "_sh" if sheared else "")
        with open(os.path.join(args.out, name + ".smt2"), "w") as f:
            f.write("; random Boolean union of axis-aligned atoms%s (gen_bool_unions.py seed %d)\n"
                    % (" in sheared coordinates" if sheared else "", args.seed))
            f.write("; exact-volume: %r\n" % float(vol))
            f.write("; exact-volume-rational: %s\n" % vol)
            f.write("; dim: %d\n" % n)
            f.write("(set-logic QF_LRA)\n")
            for nm in names:
                f.write("(declare-fun %s () Real)\n" % nm)
            f.write("(assert (and %s))\n" % box)
            f.write("(assert %s)\n" % to_smt(t, row_of, names))
            f.write("(check-sat)\n")
        made += 1
    print("wrote %d formulas to %s" % (made, args.out))


if __name__ == "__main__":
    main()
