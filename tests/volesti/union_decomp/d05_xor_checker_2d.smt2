; Checkerboard: (x < 1/2) xor (y < 1/2) inside the unit square: area 1/2.
; dim: 2
; bbox: -0.5 1.5
; exact-volume: 0.5
(set-logic QF_LRA)
(declare-const x Real)
(declare-const y Real)
(assert (and (>= x 0) (<= x 1) (>= y 0) (<= y 1)
             (xor (< x (/ 1 2)) (< y 0.5))))
(check-sat)
