; Standard simplex x_i >= 0, sum x_i <= 1 in R^2; vol = s^n/n!.
; dim: 2
; exact-volume: 0.5
(set-logic QF_LRA)
(declare-const x0 Real)
(declare-const x1 Real)
(assert (and
  (>= x0 0)
  (>= x1 0)
  (<= (+ x0 x1) 1)
))
(check-sat)
