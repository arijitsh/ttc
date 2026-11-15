; the same box twice
; family: nest
; dim: 2
; exact-volume: 2.0
; exact-volume-rational: 2
; disjuncts: 2
(set-logic QF_LRA)
(declare-const x0 Real)
(declare-const x1 Real)
(assert (or
  (and (>= x0 0) (<= x0 2) (>= x1 0) (<= x1 1))
  (and (>= x0 0) (<= x0 2) (>= x1 0) (<= x1 1))))
(check-sat)
