; [0,1000]x[0,1] u [0,1]x[0,1000] (aspect 1000)
; family: skinny
; dim: 2
; exact-volume: 1999.0
; exact-volume-rational: 1999
; disjuncts: 2
(set-logic QF_LRA)
(declare-const x0 Real)
(declare-const x1 Real)
(assert (or
  (and (>= x0 0) (<= x0 1000) (>= x1 0) (<= x1 1))
  (and (>= x0 0) (<= x0 1) (>= x1 0) (<= x1 1000))))
(check-sat)
