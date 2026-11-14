; d03 frame translated by (1e6, 1e6): area 12.
; dim: 2
; bbox: 999997 1000003
; exact-volume: 12.0
(set-logic QF_LRA)
(declare-const x Real)
(declare-const y Real)
(assert (and (<= 999998 x) (<= x 1000002) (<= 999998 y) (<= y 1000002)
             (not (and (< 999999 x) (< x 1000001) (< 999999 y) (< y 1000001)))))
(check-sat)
