; Cube [0,1000.0]^5 (size sweep); vol = side^n.
; dim: 5
; exact-volume: 1000000000000000.0
(set-logic QF_LRA)
(declare-const x0 Real)
(declare-const x1 Real)
(declare-const x2 Real)
(declare-const x3 Real)
(declare-const x4 Real)
(assert (and
  (>= x0 0)
  (<= x0 1000.0)
  (>= x1 0)
  (<= x1 1000.0)
  (>= x2 0)
  (<= x2 1000.0)
  (>= x3 0)
  (<= x3 1000.0)
  (>= x4 0)
  (<= x4 1000.0)
))
(check-sat)
