; 3 unit cubes in R^8 shifted by 1/2 along x0 and x1
; family: box
; dim: 8
; exact-volume: 2.0
; exact-volume-rational: 2
; disjuncts: 3
; timeout: 180
(set-logic QF_LRA)
(declare-const x0 Real)
(declare-const x1 Real)
(declare-const x2 Real)
(declare-const x3 Real)
(declare-const x4 Real)
(declare-const x5 Real)
(declare-const x6 Real)
(declare-const x7 Real)
(assert (or
  (and (>= x0 0) (<= x0 1) (>= x1 0) (<= x1 1) (>= x2 0) (<= x2 1) (>= x3 0) (<= x3 1) (>= x4 0) (<= x4 1) (>= x5 0) (<= x5 1) (>= x6 0) (<= x6 1) (>= x7 0) (<= x7 1))
  (and (>= x0 0.5) (<= x0 1.5) (>= x1 0) (<= x1 1) (>= x2 0) (<= x2 1) (>= x3 0) (<= x3 1) (>= x4 0) (<= x4 1) (>= x5 0) (<= x5 1) (>= x6 0) (<= x6 1) (>= x7 0) (<= x7 1))
  (and (>= x0 0) (<= x0 1) (>= x1 0.5) (<= x1 1.5) (>= x2 0) (<= x2 1) (>= x3 0) (<= x3 1) (>= x4 0) (<= x4 1) (>= x5 0) (<= x5 1) (>= x6 0) (<= x6 1) (>= x7 0) (<= x7 1))))
(check-sat)
