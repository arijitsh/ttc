; staircase: union of [0,i]x[0,6-i], i=1..5 (area 5+4+3+2+1)
; family: shape
; dim: 2
; exact-volume: 15.0
; exact-volume-rational: 15
; disjuncts: 5
(set-logic QF_LRA)
(declare-const x0 Real)
(declare-const x1 Real)
(assert (or
  (and (>= x0 0) (<= x0 1) (>= x1 0) (<= x1 5))
  (and (>= x0 0) (<= x0 2) (>= x1 0) (<= x1 4))
  (and (>= x0 0) (<= x0 3) (>= x1 0) (<= x1 3))
  (and (>= x0 0) (<= x0 4) (>= x1 0) (<= x1 2))
  (and (>= x0 0) (<= x0 5) (>= x1 0) (<= x1 1))))
(check-sat)
