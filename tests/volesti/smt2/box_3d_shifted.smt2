; Axis-aligned box [2,5] x [-1,1] x [10,21/2] in R^3.
; dim: 3
; exact-volume: 3.0
(set-logic QF_LRA)
(declare-const x0 Real)
(declare-const x1 Real)
(declare-const x2 Real)
(assert (and
  (>= x0 2)
  (<= x0 5)
  (>= x1 (- 1))
  (<= x1 1)
  (>= x2 10)
  (<= x2 (/ 21 2))
))
(check-sat)
