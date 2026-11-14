; Real-valued ite: |x| + |y| <= 1 written with (ite (>= x 0) x (- x)).  Area 2.
; dim: 2
; bbox: -2 2
; exact-volume: 2.0
(set-logic QF_LRA)
(declare-const x Real)
(declare-const y Real)
(assert (<= (+ (ite (>= x 0) x (- x)) (ite (>= y 0) y (- y))) 1))
(check-sat)
