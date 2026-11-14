; Cube [0,1000000.0]^2 (size sweep); vol = side^n.
; dim: 2
; exact-volume: 1000000000000.0
(set-logic QF_LRA)
(declare-const x0 Real)
(declare-const x1 Real)
(assert (and
  (>= x0 0)
  (<= x0 1000000.0)
  (>= x1 0)
  (<= x1 1000000.0)
))
(check-sat)
