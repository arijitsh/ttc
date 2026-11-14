; Standard simplex x_i >= 0, sum x_i <= 1 in R^8; vol = s^n/n!.
; dim: 8
; timeout: 120
; exact-volume: 2.48015873015873e-05
(set-logic QF_LRA)
(declare-const x0 Real)
(declare-const x1 Real)
(declare-const x2 Real)
(declare-const x3 Real)
(declare-const x4 Real)
(declare-const x5 Real)
(declare-const x6 Real)
(declare-const x7 Real)
(assert (and
  (>= x0 0)
  (>= x1 0)
  (>= x2 0)
  (>= x3 0)
  (>= x4 0)
  (>= x5 0)
  (>= x6 0)
  (>= x7 0)
  (<= (+ x0 x1 x2 x3 x4 x5 x6 x7) 1)
))
(check-sat)
