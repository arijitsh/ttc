; three disjoint boxes in R^5 with volumes 1, 1/32, 2
; family: nest
; dim: 5
; exact-volume: 3.03125
; exact-volume-rational: 97/32
; disjuncts: 3
; timeout: 120
(set-logic QF_LRA)
(declare-const x0 Real)
(declare-const x1 Real)
(declare-const x2 Real)
(declare-const x3 Real)
(declare-const x4 Real)
(assert (or
  (and (>= x0 0) (<= x0 1) (>= x1 0) (<= x1 1) (>= x2 0) (<= x2 1) (>= x3 0) (<= x3 1) (>= x4 0) (<= x4 1))
  (and (>= x0 2) (<= x0 2.5) (>= x1 2) (<= x1 2.5) (>= x2 2) (<= x2 2.5) (>= x3 2) (<= x3 2.5) (>= x4 2) (<= x4 2.5))
  (and (>= x0 0) (<= x0 2) (>= x1 0) (<= x1 1) (>= x2 0) (<= x2 1) (>= x3 0) (<= x3 1) (>= x4 3) (<= x4 4))))
(check-sat)
