; Cube [0,1000.0]^3 (size sweep); vol = side^n.
; dim: 3
; exact-volume: 1000000000.0
(set-logic QF_LRA)
(declare-const x0 Real)
(declare-const x1 Real)
(declare-const x2 Real)
(assert (and
  (>= x0 0)
  (<= x0 1000.0)
  (>= x1 0)
  (<= x1 1000.0)
  (>= x2 0)
  (<= x2 1000.0)
))
(check-sat)
