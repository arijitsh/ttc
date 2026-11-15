#!/usr/bin/env python3
"""Generate tests/volesti/union_smt2/*.smt2: QF_LRA unions with exactly known volume.

Every case carries an exact rational volume (Fraction) computed by an independent route
(inclusion-exclusion of boxes, exact convex clipping in 2D, combinatorics), plus a Python
membership predicate used for a Monte Carlo cross-check of the formula text.
"""
import itertools
import os
import random
import sys
from fractions import Fraction as F

OUT = sys.argv[1] if len(sys.argv) > 1 else "union_smt2"
os.makedirs(OUT, exist_ok=True)


# ---------------------------------------------------------------- SMT helpers
def num(q):
    q = F(q)
    neg = q < 0
    q = abs(q)
    d = q.denominator
    dd = d
    for p in (2, 5):
        while dd % p == 0:
            dd //= p
    if dd == 1:
        # terminating decimal
        k = 0
        while (q * 10 ** k).denominator != 1:
            k += 1
        s = str((q * 10 ** k).numerator)
        if k:
            s = s.rjust(k + 1, "0")
            s = s[:-k] + "." + s[-k:]
        else:
            s = s + ".0" if False else s
    else:
        s = "(/ %d %d)" % (q.numerator, q.denominator)
    return "(- %s)" % s if neg else s


def var(i):
    return "x%d" % i


def lin(a):
    terms = []
    for i, c in enumerate(a):
        c = F(c)
        if c == 0:
            continue
        terms.append(var(i) if c == 1 else "(* %s %s)" % (num(c), var(i)))
    if len(terms) == 1:
        return terms[0]
    return "(+ %s)" % " ".join(terms)


def hs_atom(a, b):  # a.x <= b
    return "(<= %s %s)" % (lin(a), num(b))


def box_hs(lo, hi):
    n = len(lo)
    hs = []
    for i in range(n):
        e = [0] * n
        e[i] = 1
        hs.append((list(e), F(hi[i])))
        e2 = [0] * n
        e2[i] = -1
        hs.append((e2, -F(lo[i])))
    return hs


def box_smt(lo, hi):
    parts = []
    for i in range(len(lo)):
        parts.append("(>= %s %s)" % (var(i), num(lo[i])))
        parts.append("(<= %s %s)" % (var(i), num(hi[i])))
    return "(and %s)" % " ".join(parts)


def piece_smt(hs):
    return "(and %s)" % " ".join(hs_atom(a, b) for a, b in hs)


def union_smt(pieces_smt):
    if len(pieces_smt) == 1:
        return pieces_smt[0]
    return "(or\n  %s)" % "\n  ".join(pieces_smt)


# ---------------------------------------------------------------- exact volumes
def box_vol(lo, hi):
    v = F(1)
    for l, h in zip(lo, hi):
        if h <= l:
            return F(0)
        v *= F(h) - F(l)
    return v


def box_union_vol(boxes):
    k = len(boxes)
    if k > 10:  # coordinate compression: sum over grid cells covered by some box
        n = len(boxes[0][0])
        cuts = [sorted(set(F(b[0][d]) for b in boxes) | set(F(b[1][d]) for b in boxes)) for d in range(n)]
        tot = F(0)
        for cell in itertools.product(*[range(len(c) - 1) for c in cuts]):
            mid = [(cuts[d][i] + cuts[d][i + 1]) / 2 for d, i in enumerate(cell)]
            if any(all(F(b[0][d]) <= mid[d] <= F(b[1][d]) for d in range(n)) for b in boxes):
                v = F(1)
                for d, i in enumerate(cell):
                    v *= cuts[d][i + 1] - cuts[d][i]
                tot += v
        return tot
    tot = F(0)
    for r in range(1, k + 1):
        for S in itertools.combinations(range(k), r):
            lo = [max(F(boxes[i][0][d]) for i in S) for d in range(len(boxes[0][0]))]
            hi = [min(F(boxes[i][1][d]) for i in S) for d in range(len(boxes[0][0]))]
            tot += (-1) ** (r + 1) * box_vol(lo, hi)
    return tot


def clip(poly, a, b):
    """Sutherland-Hodgman: keep {x: a.x <= b}; exact Fractions."""
    out = []
    n = len(poly)
    for i in range(n):
        P, Q = poly[i], poly[(i + 1) % n]
        fp = a[0] * P[0] + a[1] * P[1] - b
        fq = a[0] * Q[0] + a[1] * Q[1] - b
        if fp <= 0:
            out.append(P)
        if (fp < 0 < fq) or (fq < 0 < fp):
            t = fp / (fp - fq)
            out.append((P[0] + t * (Q[0] - P[0]), P[1] + t * (Q[1] - P[1])))
    return out


def hs_area2d(hs, R=10 ** 6):
    poly = [(F(-R), F(-R)), (F(R), F(-R)), (F(R), F(R)), (F(-R), F(R))]
    for a, b in hs:
        poly = clip(poly, (F(a[0]), F(a[1])), F(b))
        if not poly:
            return F(0)
    s = F(0)
    for i in range(len(poly)):
        s += poly[i][0] * poly[(i + 1) % len(poly)][1] - poly[(i + 1) % len(poly)][0] * poly[i][1]
    return abs(s) / 2


def convex_union_area2d(pieces):
    tot = F(0)
    k = len(pieces)
    for r in range(1, k + 1):
        for S in itertools.combinations(range(k), r):
            hs = [h for i in S for h in pieces[i]]
            tot += (-1) ** (r + 1) * hs_area2d(hs)
    return tot


def polygon_hs(verts):
    """CCW convex polygon -> halfplanes a.x <= b."""
    hs = []
    n = len(verts)
    for i in range(n):
        p, q = verts[i], verts[(i + 1) % n]
        dx, dy = F(q[0]) - F(p[0]), F(q[1]) - F(p[1])
        a = [dy, -dx]
        hs.append((a, a[0] * F(p[0]) + a[1] * F(p[1])))
    return hs


def inside_hs(x, hs):
    return all(sum(float(ai) * xi for ai, xi in zip(a, x)) <= float(b) + 1e-12 for a, b in hs)


# ---------------------------------------------------------------- case registry
CASES = []


def add(name, family, dim, desc, body, exact, member, bbox, timeout=None, extra=None):
    CASES.append(dict(name=name, family=family, dim=dim, desc=desc, body=body, exact=F(exact),
                      member=member, bbox=bbox, timeout=timeout, extra=extra or {}))


def add_boxes(name, family, desc, boxes, timeout=None):
    n = len(boxes[0][0])
    body = union_smt([box_smt(lo, hi) for lo, hi in boxes])
    hs = [box_hs(lo, hi) for lo, hi in boxes]
    lo = [min(F(b[0][d]) for b in boxes) for d in range(n)]
    hi = [max(F(b[1][d]) for b in boxes) for d in range(n)]
    add(name, family, n, desc, body, box_union_vol(boxes),
        lambda x, hs=hs: any(inside_hs(x, h) for h in hs), (lo, hi), timeout,
        dict(polytopes=len(boxes)))


def add_pieces2d(name, family, desc, pieces, bbox, z=None, timeout=None):
    """2D convex pieces (halfplane lists); optional extrusion z=(lo,hi) to 3D."""
    area = convex_union_area2d(pieces)
    if z is None:
        body = union_smt([piece_smt(p) for p in pieces])
        add(name, family, 2, desc, body, area,
            lambda x, ps=pieces: any(inside_hs(x, p) for p in ps), bbox, timeout,
            dict(polytopes=len(pieces)))
    else:
        p3 = [[(list(a) + [0], b) for a, b in p] for p in pieces]
        body = "(and (>= x2 %s) (<= x2 %s)\n %s)" % (num(z[0]), num(z[1]),
                                                       union_smt([piece_smt(p) for p in p3]))
        add(name, family, 3, desc, body, area * (F(z[1]) - F(z[0])),
            lambda x, ps=pieces: z[0] <= x[2] <= z[1] and any(inside_hs(x[:2], p) for p in ps),
            (bbox[0] + [F(z[0])], bbox[1] + [F(z[1])]), timeout, dict(polytopes=len(pieces)))


def rand_boxes(seed, k, n, grid=10, scale=F(1, 10)):
    rng = random.Random(seed)
    boxes = []
    for _ in range(k):
        lo, hi = [], []
        for _ in range(n):
            a = rng.randrange(0, grid - 2)
            b = rng.randrange(a + 2, grid + 1)
            lo.append(F(a) * scale)
            hi.append(F(b) * scale)
        boxes.append((lo, hi))
    return boxes


# ---- family box: unions of overlapping boxes (inclusion-exclusion)
add_boxes("box_2overlap_2d", "box", "[0,2]x[0,1] u [1,3]x[0,2]",
          [([0, 0], [2, 1]), ([1, 0], [3, 2])])
add_boxes("box_3overlap_3d", "box", "three pairwise-overlapping boxes in 3D",
          [([0, 0, 0], [2, 1, 1]), ([1, 0, 0], [3, 2, 1]), ([F(1, 2), F(1, 2), F(1, 2)], [F(5, 2), F(3, 2), 2])])
add_boxes("box_chain5_2d", "box", "chain of 5 boxes, each overlapping the next",
          [([i, 0], [i + 2, 1 + (i % 2)]) for i in range(5)])
add_boxes("box_10rand_3d", "box", "10 random boxes on a 0.1-grid in [0,1]^3",
          rand_boxes(3, 10, 3))
add_boxes("box_6rand_4d", "box", "6 random boxes on a 0.1-grid in [0,1]^4",
          rand_boxes(4, 6, 4))
add_boxes("box_4rand_6d", "box", "4 random boxes on a 0.1-grid in [0,1]^6",
          rand_boxes(6, 4, 6), timeout=120)
add_boxes("box_3overlap_8d", "box", "3 unit cubes in R^8 shifted by 1/2 along x0 and x1",
          [([0] * 8, [1] * 8), ([F(1, 2)] + [0] * 7, [F(3, 2)] + [1] * 7),
           ([0, F(1, 2)] + [0] * 6, [1, F(3, 2)] + [1] * 6)], timeout=180)

# ---- family nest: nested / identical / disjoint / touching
add_boxes("nest_identical_2d", "nest", "the same box twice", [([0, 0], [2, 1]), ([0, 0], [2, 1])])
add_boxes("nest_small_first_3d", "nest", "[0.2,0.5]^3 then [0,1]^3 (small inside big, small listed first)",
          [([F(1, 5)] * 3, [F(1, 2)] * 3), ([0] * 3, [1] * 3)])
add_boxes("nest_big_first_3d", "nest", "[0,1]^3 then [0.2,0.5]^3 (big listed first)",
          [([0] * 3, [1] * 3), ([F(1, 5)] * 3, [F(1, 2)] * 3)])
add_boxes("nest_disjoint2_2d", "nest", "two disjoint boxes", [([0, 0], [1, 1]), ([2, 0], [4, 1])])
add_boxes("nest_disjoint3_5d", "nest", "three disjoint boxes in R^5 with volumes 1, 1/32, 2",
          [([0] * 5, [1] * 5), ([2] * 5, [F(5, 2)] * 5), ([0, 0, 0, 0, 3], [2, 1, 1, 1, 4])], timeout=120)
add_boxes("nest_touching_2d", "nest", "two unit squares sharing the face x=1",
          [([0, 0], [1, 1]), ([1, 0], [2, 1])])

# ---- family shape: L-shape, frames, staircases (Boolean structure with not/or)
B3 = "(and (>= x0 0) (<= x0 3) (>= x1 0) (<= x1 3))"
add("shape_lshape_2d", "shape", 2, "L-shape: [0,3]^2 minus (1,3]x(1,3], written with not",
    "(and %s\n  (not (and (> x0 1) (> x1 1))))" % B3, 5,
    lambda x: 0 <= x[0] <= 3 and 0 <= x[1] <= 3 and not (x[0] > 1 and x[1] > 1), ([0, 0], [3, 3]))
add("shape_frame_2d", "shape", 2, "frame: [0,3]^2 minus (1,2)^2, written with not",
    "(and %s\n  (not (and (> x0 1) (< x0 2) (> x1 1) (< x1 2))))" % B3, 8,
    lambda x: 0 <= x[0] <= 3 and 0 <= x[1] <= 3 and not (1 < x[0] < 2 and 1 < x[1] < 2), ([0, 0], [3, 3]))
add("shape_frame_3d", "shape", 3, "hollow cube: [0,3]^3 minus (1,2)^3, written with not (6 overlapping slabs)",
    "(and (>= x0 0) (<= x0 3) (>= x1 0) (<= x1 3) (>= x2 0) (<= x2 3)\n"
    "  (not (and (> x0 1) (< x0 2) (> x1 1) (< x1 2) (> x2 1) (< x2 2))))", 26,
    lambda x: all(0 <= t <= 3 for t in x) and not all(1 < t < 2 for t in x), ([0] * 3, [3] * 3))
add_boxes("shape_staircase_2d", "shape", "staircase: union of [0,i]x[0,6-i], i=1..5 (area 5+4+3+2+1)",
          [([0, 0], [i, 6 - i]) for i in range(1, 6)])
add_boxes("shape_staircase_3d", "shape", "3D staircase: union of origin-anchored boxes",
          [([0, 0, 0], [3, 1, 1]), ([0, 0, 0], [2, 2, 1]), ([0, 0, 0], [1, 3, 2]), ([0, 0, 0], [1, 1, 3])])

# ---- family tri: triangles / polygons
add_pieces2d("tri_two_square_2d", "tri", "unit square as two triangles split by x0+x1=1",
             [polygon_hs([(0, 0), (1, 0), (0, 1)]), polygon_hs([(1, 0), (1, 1), (0, 1)])], ([0, 0], [1, 1]))
HEX = [(2, 0), (1, 2), (-1, 2), (-2, 0), (-1, -2), (1, -2)]
add_pieces2d("tri_hexagon6_2d", "tri", "hexagon (2,0),(1,2),(-1,2),(-2,0),(-1,-2),(1,-2) as 6 triangles from 0",
             [polygon_hs([(0, 0), HEX[i], HEX[(i + 1) % 6]]) for i in range(6)], ([-2, -2], [2, 2]))
add_pieces2d("tri_star2_2d", "tri", "two overlapping triangles (star-of-David-like)",
             [polygon_hs([(0, 0), (4, 0), (2, 3)]), polygon_hs([(0, 2), (2, -1), (4, 2)])], ([0, -1], [4, 3]))


# ---- family skinny: skinny slabs, axis-aligned and rotated
def slab(d, w, L, c=(0, 0)):
    """Parallelogram centred at c: |d.(x-c)| <= w|d|/2 (thin, width w), |dp.(x-c)| <= L|d|/2, dp = d rotated."""
    d = (F(d[0]), F(d[1]))
    dp = (-d[1], d[0])
    nd2 = d[0] ** 2 + d[1] ** 2
    # use |d| only through |d|^2: choose integer d with |d| integer (3,4),(5,12),(1,0),...
    import math
    nd = F(math.isqrt(int(nd2)))
    assert nd * nd == nd2
    hs = []
    for v, half in ((d, F(w) * nd / 2), (dp, F(L) * nd / 2)):
        off = v[0] * F(c[0]) + v[1] * F(c[1])
        hs.append(([v[0], v[1]], off + half))
        hs.append(([-v[0], -v[1]], -off + half))
    return hs


add_boxes("skinny_cross_2d", "skinny", "[0,1000]x[0,1] u [0,1]x[0,1000] (aspect 1000)",
          [([0, 0], [1000, 1]), ([0, 0], [1, 1000])])
add_pieces2d("skinny_rot2_2d", "skinny", "two crossing slabs 1x100 along (4,-3)/5 and (4,3)/5",
             [slab((3, 4), 1, 100), slab((-3, 4), 1, 100)], ([-50, -50], [50, 50]))
add_pieces2d("skinny_star3_2d", "skinny", "three crossing slabs 1x300 (aspect 300), directions (1,0),(3,4),(-3,4)",
             [slab((0, 1), 1, 300), slab((4, 3), 1, 300), slab((-4, 3), 1, 300)], ([-150, -150], [150, 150]))
add_pieces2d("skinny_rot_offset_2d", "skinny", "two parallel rotated slabs 1x1000 overlapping in half their width",
             [slab((5, 12), 1, 1000), slab((5, 12), 1, 1000, c=(F(5, 26), F(12, 26)))], ([-500, -500], [500, 500]))
add_pieces2d("skinny_rot2_3d", "skinny", "rotated crossing slabs 1x1000 extruded by x2 in [0,2]",
             [slab((3, 4), 1, 1000), slab((-3, 4), 1, 1000)], ([-500, -500], [500, 500]), z=(0, 2))
add_boxes("skinny_parallel_3d", "skinny", "[0,100]x[0,1]x[0,1] u [50,150]x[0.5,1.5]x[0,1]",
          [([0, 0, 0], [100, 1, 1]), ([50, F(1, 2), 0], [150, F(3, 2), 1])])

# ---- family tiny: many tiny disjoint boxes + one large box
tiny40 = [([F(2) + F(i % 10) / 5, F(i // 10) / 5], [F(2) + F(i % 10) / 5 + F(1, 10), F(i // 10) / 5 + F(1, 10)])
          for i in range(40)]
add_boxes("tiny40_big_last_2d", "tiny", "40 disjoint 0.1x0.1 boxes then the unit square (disjoint)",
          tiny40 + [([0, 0], [1, 1])], timeout=120)
add_boxes("tiny40_big_first_2d", "tiny", "the unit square then 40 disjoint 0.1x0.1 boxes",
          [([0, 0], [1, 1])] + tiny40, timeout=120)
tiny_in = [([F(i % 5) / 5 + F(1, 20)] * 1 + [F(i // 5) / 5 + F(1, 20)] + [F(1, 2)],
            [F(i % 5) / 5 + F(1, 10)] + [F(i // 5) / 5 + F(1, 10)] + [F(11, 20)]) for i in range(25)]
add_boxes("tiny25_inside_big_3d", "tiny", "25 boxes of side 0.05 inside [0,1]^3, then [0,1]^3",
          tiny_in + [([0] * 3, [1] * 3)], timeout=120)
add_boxes("tiny20_sum_2d", "tiny", "20 disjoint boxes with sides 1e-3..1 (volumes span 1e-6..1)",
          [([3 * i, 0], [3 * i + F(1, 2 ** (i // 2)), F(1, 2 ** ((i + 1) // 2))]) for i in range(20)], timeout=120)

# ---- family scale: dims and scales
add_boxes("scale_1e3_2d", "scale", "[0,1000]x[0,500] u [500,1500]x[0,1000] (volume 1.25e6)",
          [([0, 0], [1000, 500]), ([500, 0], [1500, 1000])])
add_boxes("scale_1em3_3d", "scale", "boxes with sides 1e-3 (volume 5e-9)",
          [([0, 0, 0], [F(2, 1000), F(1, 1000), F(1, 1000)]), ([F(1, 1000), 0, 0], [F(3, 1000), F(2, 1000), F(1, 1000)])])
add_boxes("scale_mixed_2d", "scale", "[0,1000]x[0,0.001] u [0,1]x[0,1] (sides 1e-3..1e3)",
          [([0, 0], [1000, F(1, 1000)]), ([0, 0], [1, 1])])


# ---- family parity: xor of threshold atoms in the unit cube
def parity_case(name, thr, timeout=None):
    n = len(thr)
    cube = " ".join("(>= x%d 0) (<= x%d 1)" % (i, i) for i in range(n))
    atoms = " ".join("(< x%d %s)" % (i, num(t)) for i, t in enumerate(thr))
    prod = F(1)
    for t in thr:
        prod *= 1 - 2 * F(t)
    exact = (1 - prod) / 2
    add(name, "parity", n, "unit cube, xor of x_i < t_i, t = %s (2^%d cubes)" % (
        ",".join(str(float(t)) for t in thr), n - 1),
        "(and %s\n  (xor %s))" % (cube, atoms), exact,
        lambda x, thr=thr: sum(1 for xi, t in zip(x, thr) if xi < t) % 2 == 1, ([0] * n, [1] * n), timeout)


parity_case("parity3_3d", [F(1, 2)] * 3)
parity_case("parity4_uneven_4d", [F(3, 10), F(3, 5), F(1, 4), F(4, 5)])
parity_case("parity5_5d", [F(1, 5), F(2, 5), F(1, 2), F(7, 10), F(9, 10)], timeout=120)
add("parity_shared_2d", "parity", 2, "unit square, (x0<0.3) xor (x0<0.7) xor (x1<0.2): atoms share x0",
    "(and (>= x0 0) (<= x0 1) (>= x1 0) (<= x1 1)\n  (xor (< x0 0.3) (< x0 0.7) (< x1 0.2)))",
    F(1, 2) * (1 - (1 - 2 * F(2, 5)) * (1 - 2 * F(1, 5))),
    lambda x: ((x[0] < 0.3) + (x[0] < 0.7) + (x[1] < 0.2)) % 2 == 1, ([0, 0], [1, 1]))


# ---------------------------------------------------------------- MC check + write
def mc_check(c, N=50000, seed=1):
    rng = random.Random(seed)
    lo = [float(v) for v in c["bbox"][0]]
    hi = [float(v) for v in c["bbox"][1]]
    vbox = 1.0
    for l, h in zip(lo, hi):
        vbox *= h - l
    hits = 0
    for _ in range(N):
        x = [rng.uniform(l, h) for l, h in zip(lo, hi)]
        hits += c["member"](x)
    p = hits / N
    est = p * vbox
    se = vbox * (p * (1 - p) / N) ** 0.5
    ex = float(c["exact"])
    return est, se, ex


bad = 0
for c in CASES:
    est, se, ex = mc_check(c)
    z = (est - ex) / se if se > 0 else 0.0
    flag = "" if abs(z) < 4 else "  <-- MC MISMATCH"
    bad += bool(flag)
    print("%-26s %-7s dim %d exact %-14.8g mc %-14.8g z %+5.2f%s" % (c["name"], c["family"], c["dim"], ex, est, z, flag))
    vars_ = "\n".join("(declare-const x%d Real)" % i for i in range(c["dim"]))
    hdr = ["; %s" % c["desc"],
           "; family: %s" % c["family"],
           "; dim: %d" % c["dim"],
           "; exact-volume: %s" % repr(float(c["exact"])),
           "; exact-volume-rational: %s" % c["exact"]]
    if "polytopes" in c["extra"]:
        hdr.append("; disjuncts: %d" % c["extra"]["polytopes"])
    if c["timeout"]:
        hdr.append("; timeout: %d" % c["timeout"])
    txt = "\n".join(hdr) + "\n(set-logic QF_LRA)\n" + vars_ + "\n(assert " + c["body"] + ")\n(check-sat)\n"
    with open(os.path.join(OUT, c["name"] + ".smt2"), "w") as f:
        f.write(txt)
print("%d cases, %d MC mismatches" % (len(CASES), bad))
