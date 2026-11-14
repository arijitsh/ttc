; Axis-aligned box [0,1] x [0,1] x [0,1] x [0,1] x [0,1] in R^5.
; dim: 5
; exact-volume: 1.0
(set-logic QF_LRA)
(declare-const x0 Real)
(declare-const x1 Real)
(declare-const x2 Real)
(declare-const x3 Real)
(declare-const x4 Real)
(assert (and
  (>= x0 0)
  (<= x0 1)
  (>= x1 0)
  (<= x1 1)
  (>= x2 0)
  (<= x2 1)
  (>= x3 0)
  (<= x3 1)
  (>= x4 0)
  (<= x4 1)
))
(check-sat)
