; [0,1000]x[0,500] u [500,1500]x[0,1000] (volume 1.25e6)
; family: scale
; dim: 2
; exact-volume: 1250000.0
; exact-volume-rational: 1250000
; disjuncts: 2
(set-logic QF_LRA)
(declare-const x0 Real)
(declare-const x1 Real)
(assert (or
  (and (>= x0 0) (<= x0 1000) (>= x1 0) (<= x1 500))
  (and (>= x0 500) (<= x0 1500) (>= x1 0) (<= x1 1000))))
(check-sat)
