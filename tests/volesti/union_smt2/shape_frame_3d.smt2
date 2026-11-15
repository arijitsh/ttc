; hollow cube: [0,3]^3 minus (1,2)^3, written with not (6 overlapping slabs)
; family: shape
; dim: 3
; exact-volume: 26.0
; exact-volume-rational: 26
(set-logic QF_LRA)
(declare-const x0 Real)
(declare-const x1 Real)
(declare-const x2 Real)
(assert (and (>= x0 0) (<= x0 3) (>= x1 0) (<= x1 3) (>= x2 0) (<= x2 3)
  (not (and (> x0 1) (< x0 2) (> x1 1) (< x1 2) (> x2 1) (< x2 2)))))
(check-sat)
