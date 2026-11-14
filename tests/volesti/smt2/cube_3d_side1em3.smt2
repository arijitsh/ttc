; Cube [0,0.001]^3 (size sweep); vol = side^n.
; dim: 3
; exact-volume: 1e-09
(set-logic QF_LRA)
(declare-const x0 Real)
(declare-const x1 Real)
(declare-const x2 Real)
(assert (and
  (>= x0 0)
  (<= x0 0.001)
  (>= x1 0)
  (<= x1 0.001)
  (>= x2 0)
  (<= x2 0.001)
))
(check-sat)
