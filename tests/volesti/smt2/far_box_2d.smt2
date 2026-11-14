; Box [1e6, 1e6+1] x [-1e6, -1e6+2] far from the origin; area 2.
; dim: 2
; exact-volume: 2.0
(set-logic QF_LRA)
(declare-const x0 Real)
(declare-const x1 Real)
(assert (and
  (>= x0 1000000)
  (<= x0 1000001)
  (>= x1 (- 1000000))
  (<= x1 (- 999998))
))
(check-sat)
