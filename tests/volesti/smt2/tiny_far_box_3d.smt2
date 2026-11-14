; Box [1e6, 1e6+1e-3] x [0,1] x [0,1]; vol 1e-3 (far + thin).
; dim: 3
; exact-volume: 0.001
(set-logic QF_LRA)
(declare-const x0 Real)
(declare-const x1 Real)
(declare-const x2 Real)
(assert (and
  (>= x0 1000000.0)
  (<= x0 1000000.001)
  (>= x1 0)
  (<= x1 1)
  (>= x2 0)
  (<= x2 1)
))
(check-sat)
