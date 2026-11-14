; Hexagon |y|<=1, |x+y|<=2, |x-y|<=2 with vertices (+-2,0),(+-1,+-1); area 6.
; dim: 2
; exact-volume: 6.0
(set-logic QF_LRA)
(declare-const x0 Real)
(declare-const x1 Real)
(assert (and
  (<= x1 1)
  (>= x1 (- 1))
  (<= (+ x0 x1) 2)
  (>= (+ x0 x1) (- 2))
  (<= (- x0 x1) 2)
  (>= (- x0 x1) (- 2))
))
(check-sat)
