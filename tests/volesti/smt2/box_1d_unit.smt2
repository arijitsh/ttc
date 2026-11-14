; Axis-aligned box [0,1] in R^1.
; dim: 1
; exact-volume: 1.0
(set-logic QF_LRA)
(declare-const x0 Real)
(assert (and
  (>= x0 0)
  (<= x0 1)
))
(check-sat)
