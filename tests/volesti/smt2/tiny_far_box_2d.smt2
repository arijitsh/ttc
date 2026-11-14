; Tiny box [1e6, 1e6+1e-3] x [1e6, 1e6+1e-3] far from the origin; area 1e-6.
; dim: 2
; exact-volume: 1e-06
(set-logic QF_LRA)
(declare-const x0 Real)
(declare-const x1 Real)
(assert (and
  (>= x0 1000000.0)
  (<= x0 1000000.001)
  (>= x1 1000000.0)
  (<= x1 1000000.001)
))
(check-sat)
