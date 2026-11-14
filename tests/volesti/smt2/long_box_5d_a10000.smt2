; Long box [0,1]^4 x [0,10000] (aspect ratio 10000).
; dim: 5
; exact-volume: 10000.0
(set-logic QF_LRA)
(declare-const x0 Real)
(declare-const x1 Real)
(declare-const x2 Real)
(declare-const x3 Real)
(declare-const x4 Real)
(assert (and
  (>= x0 0)
  (<= x0 1)
  (>= x1 0)
  (<= x1 1)
  (>= x2 0)
  (<= x2 1)
  (>= x3 0)
  (<= x3 1)
  (>= x4 0)
  (<= x4 10000)
))
(check-sat)
