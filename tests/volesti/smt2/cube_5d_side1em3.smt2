; Cube [0,0.001]^5 (size sweep); vol = side^n.
; dim: 5
; exact-volume: 1e-15
(set-logic QF_LRA)
(declare-const x0 Real)
(declare-const x1 Real)
(declare-const x2 Real)
(declare-const x3 Real)
(declare-const x4 Real)
(assert (and
  (>= x0 0)
  (<= x0 0.001)
  (>= x1 0)
  (<= x1 0.001)
  (>= x2 0)
  (<= x2 0.001)
  (>= x3 0)
  (<= x3 0.001)
  (>= x4 0)
  (<= x4 0.001)
))
(check-sat)
