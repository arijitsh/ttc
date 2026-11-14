; Simplex x_i >= 0, sum x_i <= 1000.0 in R^2 (size sweep); vol = s^n/n!.
; dim: 2
; exact-volume: 500000.0
(set-logic QF_LRA)
(declare-const x0 Real)
(declare-const x1 Real)
(assert (and
  (>= x0 0)
  (>= x1 0)
  (<= (+ x0 x1) 1000.0)
))
(check-sat)
