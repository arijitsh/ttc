; Standard simplex x_i >= 0, sum x_i <= 2 in R^3; vol = s^n/n!.
; dim: 3
; exact-volume: 1.3333333333333333
(set-logic QF_LRA)
(declare-const x0 Real)
(declare-const x1 Real)
(declare-const x2 Real)
(assert (and
  (>= x0 0)
  (>= x1 0)
  (>= x2 0)
  (<= (+ x0 x1 x2) 2)
))
(check-sat)
