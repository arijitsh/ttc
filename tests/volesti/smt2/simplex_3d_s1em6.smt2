; Simplex x_i >= 0, sum x_i <= 0.000001 in R^3 (size sweep); vol = s^n/n!.
; dim: 3
; exact-volume: 1.6666666666666667e-19
(set-logic QF_LRA)
(declare-const x0 Real)
(declare-const x1 Real)
(declare-const x2 Real)
(assert (and
  (>= x0 0)
  (>= x1 0)
  (>= x2 0)
  (<= (+ x0 x1 x2) 0.000001)
))
(check-sat)
