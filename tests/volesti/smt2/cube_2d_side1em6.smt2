; Cube [0,0.000001]^2 (size sweep); vol = side^n.
; dim: 2
; exact-volume: 1e-12
(set-logic QF_LRA)
(declare-const x0 Real)
(declare-const x1 Real)
(assert (and
  (>= x0 0)
  (<= x0 0.000001)
  (>= x1 0)
  (<= x1 0.000001)
))
(check-sat)
