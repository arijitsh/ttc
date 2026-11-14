; Segment x0 = 1/2, 0 <= x1 <= 1: lower-dimensional; 2D volume is 0.
; Report-only: we document what ttc prints.
; dim: 2
; exact-volume: 0.0
; report-only: true
(set-logic QF_LRA)
(declare-const x0 Real)
(declare-const x1 Real)
(assert (and
  (= x0 0.5)
  (>= x1 0)
  (<= x1 1)
))
(check-sat)
