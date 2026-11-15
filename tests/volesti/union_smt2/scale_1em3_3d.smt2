; boxes with sides 1e-3 (volume 5e-9)
; family: scale
; dim: 3
; exact-volume: 5e-09
; exact-volume-rational: 1/200000000
; disjuncts: 2
(set-logic QF_LRA)
(declare-const x0 Real)
(declare-const x1 Real)
(declare-const x2 Real)
(assert (or
  (and (>= x0 0) (<= x0 0.002) (>= x1 0) (<= x1 0.001) (>= x2 0) (<= x2 0.001))
  (and (>= x0 0.001) (<= x0 0.003) (>= x1 0) (<= x1 0.002) (>= x2 0) (<= x2 0.001))))
(check-sat)
