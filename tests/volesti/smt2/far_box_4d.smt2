; Box [1e6,1e6+1]^3 x [5e5, 5e5+3] far from the origin; vol 3.
; dim: 4
; exact-volume: 3.0
(set-logic QF_LRA)
(declare-const x0 Real)
(declare-const x1 Real)
(declare-const x2 Real)
(declare-const x3 Real)
(assert (and
  (>= x0 1000000)
  (<= x0 1000001)
  (>= x1 1000000)
  (<= x1 1000001)
  (>= x2 1000000)
  (<= x2 1000001)
  (>= x3 500000)
  (<= x3 500003)
))
(check-sat)
