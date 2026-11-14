; Long box [0,1]^1 x [0,10000] (aspect ratio 10000).
; dim: 2
; exact-volume: 10000.0
(set-logic QF_LRA)
(declare-const x0 Real)
(declare-const x1 Real)
(assert (and
  (>= x0 0)
  (<= x0 1)
  (>= x1 0)
  (<= x1 10000)
))
(check-sat)
