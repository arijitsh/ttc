; [0,100]x[0,1]x[0,1] u [50,150]x[0.5,1.5]x[0,1]
; family: skinny
; dim: 3
; exact-volume: 175.0
; exact-volume-rational: 175
; disjuncts: 2
(set-logic QF_LRA)
(declare-const x0 Real)
(declare-const x1 Real)
(declare-const x2 Real)
(assert (or
  (and (>= x0 0) (<= x0 100) (>= x1 0) (<= x1 1) (>= x2 0) (<= x2 1))
  (and (>= x0 50) (<= x0 150) (>= x1 0.5) (<= x1 1.5) (>= x2 0) (<= x2 1))))
(check-sat)
