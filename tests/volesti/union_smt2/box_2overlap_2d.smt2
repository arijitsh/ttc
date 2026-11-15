; [0,2]x[0,1] u [1,3]x[0,2]
; family: box
; dim: 2
; exact-volume: 5.0
; exact-volume-rational: 5
; disjuncts: 2
(set-logic QF_LRA)
(declare-const x0 Real)
(declare-const x1 Real)
(assert (or
  (and (>= x0 0) (<= x0 2) (>= x1 0) (<= x1 1))
  (and (>= x0 1) (<= x0 3) (>= x1 0) (<= x1 2))))
(check-sat)
