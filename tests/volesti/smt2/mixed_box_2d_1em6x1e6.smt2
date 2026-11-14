; Mixed-size box [0,0.000001] x [0,1000000.0]; area 1 (aspect 1e12).
; dim: 2
; exact-volume: 1.0
(set-logic QF_LRA)
(declare-const x0 Real)
(declare-const x1 Real)
(assert (and
  (>= x0 0)
  (<= x0 0.000001)
  (>= x1 0)
  (<= x1 1000000.0)
))
(check-sat)
