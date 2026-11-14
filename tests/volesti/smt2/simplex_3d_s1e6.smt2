; Simplex x_i >= 0, sum x_i <= 1000000.0 in R^3 (size sweep); vol = s^n/n!.
; dim: 3
; exact-volume: 1.6666666666666666e+17
(set-logic QF_LRA)
(declare-const x0 Real)
(declare-const x1 Real)
(declare-const x2 Real)
(assert (and
  (>= x0 0)
  (>= x1 0)
  (>= x2 0)
  (<= (+ x0 x1 x2) 1000000.0)
))
(check-sat)
