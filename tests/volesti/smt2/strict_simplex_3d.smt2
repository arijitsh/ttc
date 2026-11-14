; Open simplex x_i > 0, sum < 1 in R^3; vol 1/6.
; dim: 3
; exact-volume: 0.16666666666666666
(set-logic QF_LRA)
(declare-const x0 Real)
(declare-const x1 Real)
(declare-const x2 Real)
(assert (and
  (> x0 0)
  (> x1 0)
  (> x2 0)
  (< (+ x0 x1 x2) 1)
))
(check-sat)
