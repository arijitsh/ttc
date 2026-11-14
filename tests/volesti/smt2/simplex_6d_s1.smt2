; Standard simplex x_i >= 0, sum x_i <= 1 in R^6; vol = s^n/n!.
; dim: 6
; exact-volume: 0.001388888888888889
(set-logic QF_LRA)
(declare-const x0 Real)
(declare-const x1 Real)
(declare-const x2 Real)
(declare-const x3 Real)
(declare-const x4 Real)
(declare-const x5 Real)
(assert (and
  (>= x0 0)
  (>= x1 0)
  (>= x2 0)
  (>= x3 0)
  (>= x4 0)
  (>= x5 0)
  (<= (+ x0 x1 x2 x3 x4 x5) 1)
))
(check-sat)
