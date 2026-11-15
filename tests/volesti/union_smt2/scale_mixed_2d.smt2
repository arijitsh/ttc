; [0,1000]x[0,0.001] u [0,1]x[0,1] (sides 1e-3..1e3)
; family: scale
; dim: 2
; exact-volume: 1.999
; exact-volume-rational: 1999/1000
; disjuncts: 2
(set-logic QF_LRA)
(declare-const x0 Real)
(declare-const x1 Real)
(assert (or
  (and (>= x0 0) (<= x0 1000) (>= x1 0) (<= x1 0.001))
  (and (>= x0 0) (<= x0 1) (>= x1 0) (<= x1 1))))
(check-sat)
