; two disjoint boxes
; family: nest
; dim: 2
; exact-volume: 3.0
; exact-volume-rational: 3
; disjuncts: 2
(set-logic QF_LRA)
(declare-const x0 Real)
(declare-const x1 Real)
(assert (or
  (and (>= x0 0) (<= x0 1) (>= x1 0) (<= x1 1))
  (and (>= x0 2) (<= x0 4) (>= x1 0) (<= x1 1))))
(check-sat)
