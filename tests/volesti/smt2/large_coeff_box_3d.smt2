; Unit cube written as 0 <= 1e6*x0 <= 1e6, 0 <= 0.000001*x1 <= 0.000001, 0 <= 1000*x2 <= 1000.
; dim: 3
; exact-volume: 1.0
(set-logic QF_LRA)
(declare-const x0 Real)
(declare-const x1 Real)
(declare-const x2 Real)
(assert (and
  (>= (* 1000000 x0) 0)
  (<= (* 1000000 x0) 1000000)
  (>= (* 0.000001 x1) 0)
  (<= (* 0.000001 x1) 0.000001)
  (>= (* 1000 x2) 0)
  (<= (* 1000 x2) 1000)
))
(check-sat)
