; Mixed-size box [0,0.001] x [0,1] x [0,1000.0]; vol 1 (aspect 1e6).
; dim: 3
; exact-volume: 1.0
(set-logic QF_LRA)
(declare-const x0 Real)
(declare-const x1 Real)
(declare-const x2 Real)
(assert (and
  (>= x0 0)
  (<= x0 0.001)
  (>= x1 0)
  (<= x1 1)
  (>= x2 0)
  (<= x2 1000.0)
))
(check-sat)
