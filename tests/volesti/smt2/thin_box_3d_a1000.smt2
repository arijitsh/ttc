; Thin slab [0,0.001] x [0,1]^2 (aspect 1000).
; dim: 3
; exact-volume: 0.001
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
  (<= x2 1)
))
(check-sat)
