; 3D staircase: union of origin-anchored boxes
; family: shape
; dim: 3
; exact-volume: 10.0
; exact-volume-rational: 10
; disjuncts: 4
(set-logic QF_LRA)
(declare-const x0 Real)
(declare-const x1 Real)
(declare-const x2 Real)
(assert (or
  (and (>= x0 0) (<= x0 3) (>= x1 0) (<= x1 1) (>= x2 0) (<= x2 1))
  (and (>= x0 0) (<= x0 2) (>= x1 0) (<= x1 2) (>= x2 0) (<= x2 1))
  (and (>= x0 0) (<= x0 1) (>= x1 0) (<= x1 3) (>= x2 0) (<= x2 2))
  (and (>= x0 0) (<= x0 1) (>= x1 0) (<= x1 1) (>= x2 0) (<= x2 3))))
(check-sat)
