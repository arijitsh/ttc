; three pairwise-overlapping boxes in 3D
; family: box
; dim: 3
; exact-volume: 7.125
; exact-volume-rational: 57/8
; disjuncts: 3
(set-logic QF_LRA)
(declare-const x0 Real)
(declare-const x1 Real)
(declare-const x2 Real)
(assert (or
  (and (>= x0 0) (<= x0 2) (>= x1 0) (<= x1 1) (>= x2 0) (<= x2 1))
  (and (>= x0 1) (<= x0 3) (>= x1 0) (<= x1 2) (>= x2 0) (<= x2 1))
  (and (>= x0 0.5) (<= x0 2.5) (>= x1 0.5) (<= x1 1.5) (>= x2 0.5) (<= x2 2))))
(check-sat)
