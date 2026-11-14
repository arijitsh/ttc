; Axis-aligned box [-1,1] x [-1,1] x [-1,1] x [-1,1] in R^4.
; dim: 4
; exact-volume: 16.0
(set-logic QF_LRA)
(declare-const x0 Real)
(declare-const x1 Real)
(declare-const x2 Real)
(declare-const x3 Real)
(assert (and
  (>= x0 (- 1))
  (<= x0 1)
  (>= x1 (- 1))
  (<= x1 1)
  (>= x2 (- 1))
  (<= x2 1)
  (>= x3 (- 1))
  (<= x3 1)
))
(check-sat)
