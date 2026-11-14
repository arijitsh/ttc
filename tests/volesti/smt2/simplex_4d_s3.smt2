; Standard simplex x_i >= 0, sum x_i <= 3 in R^4; vol = s^n/n!.
; dim: 4
; exact-volume: 3.375
(set-logic QF_LRA)
(declare-const x0 Real)
(declare-const x1 Real)
(declare-const x2 Real)
(declare-const x3 Real)
(assert (and
  (>= x0 0)
  (>= x1 0)
  (>= x2 0)
  (>= x3 0)
  (<= (+ x0 x1 x2 x3) 3)
))
(check-sat)
