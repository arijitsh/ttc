; Long box [0,1]^2 x [0,1000] (aspect ratio 1000).
; dim: 3
; exact-volume: 1000.0
(set-logic QF_LRA)
(declare-const x0 Real)
(declare-const x1 Real)
(declare-const x2 Real)
(assert (and
  (>= x0 0)
  (<= x0 1)
  (>= x1 0)
  (<= x1 1)
  (>= x2 0)
  (<= x2 1000)
))
(check-sat)
