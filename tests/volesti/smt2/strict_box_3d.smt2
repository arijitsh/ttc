; Box (0,1)x(0,2)x(0,3) with strict inequalities; vol 6.
; dim: 3
; exact-volume: 6.0
(set-logic QF_LRA)
(declare-const x0 Real)
(declare-const x1 Real)
(declare-const x2 Real)
(assert (and
  (> x0 0)
  (< x0 1)
  (> x1 0)
  (< x1 2)
  (> x2 0)
  (< x2 3)
))
(check-sat)
