; chain of 5 boxes, each overlapping the next
; family: box
; dim: 2
; exact-volume: 10.0
; exact-volume-rational: 10
; disjuncts: 5
(set-logic QF_LRA)
(declare-const x0 Real)
(declare-const x1 Real)
(assert (or
  (and (>= x0 0) (<= x0 2) (>= x1 0) (<= x1 1))
  (and (>= x0 1) (<= x0 3) (>= x1 0) (<= x1 2))
  (and (>= x0 2) (<= x0 4) (>= x1 0) (<= x1 1))
  (and (>= x0 3) (<= x0 5) (>= x1 0) (<= x1 2))
  (and (>= x0 4) (<= x0 6) (>= x1 0) (<= x1 1))))
(check-sat)
