; Simplex x_i >= 0, sum x_i <= 0.000001 in R^5 (size sweep); vol = s^n/n!.
; dim: 5
; exact-volume: 8.333333333333334e-33
(set-logic QF_LRA)
(declare-const x0 Real)
(declare-const x1 Real)
(declare-const x2 Real)
(declare-const x3 Real)
(declare-const x4 Real)
(assert (and
  (>= x0 0)
  (>= x1 0)
  (>= x2 0)
  (>= x3 0)
  (>= x4 0)
  (<= (+ x0 x1 x2 x3 x4) 0.000001)
))
(check-sat)
